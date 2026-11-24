function R = FABER_RAAN_Phase_Joint30d(cfg)
% FABER_RAAN_PHASE_JOINT30D
% Per-event, same-altitude/same-inclination J2 + phasing accessibility.
%
% Run directly:
%   R = FABER_RAAN_Phase_Joint30d();
%
% Or customize:
%   cfg = struct('DVcap_mps',100,'Tmax_days',30);
%   R = FABER_RAAN_Phase_Joint30d(cfg);
%
% REFERENCE: circular host/client orbits h=480.1 km, i=51.8 deg.
% TRAJECTORY ARCHITECTURE (four ideal tangential burns):
%   1. Hohmann half-ellipse to circular drift orbit aD (time tH)
%   2. Drift at aD for time tD
%   3. Hohmann half-ellipse back to reference orbit (time tH)
%   4. OPTIONAL small residual plane correction at a node, cost estimated
%      from final plane-angle difference. This correction is an IDEALIZED
%      screening allowance, not an executed, time-tagged node burn.
%
% COUPLED TERMINAL CONDITIONS (same final epoch T=2*tH+tD):
%   Δu(T) = wrap(Δu0 + 2*pi - 2*n0*tH + (nD-n0)*tD) = 0
%   ΔΩ(T) = wrap(ΔΩ0 + 2*(ΩdotH-Ωdot0)*tH
%                         + (ΩdotD-Ωdot0)*tD)
%   Phase closure is imposed analytically for each 2*pi revolution branch;
%   drift altitude and branch are searched jointly with residual plane cost.
%   All angular differences are relative HOST MINUS CLIENT at initial epoch.
%
% IMPORTANT PHYSICS LIMITATIONS:
%   * Mean circular orbit / secular J2, Hohmann arcs treated as Keplerian.
%   * Mean motion used for along-orbit phasing, rather than a high-order
%     J2-corrected mean argument-of-latitude rate.
%   * First-order averaged J2 RAAN drift integrated during Hohmann arcs
%     using ellipse aT and pT; no J2 short-period position correction.
%   * True 3D off-plane geometry and the position of the residual plane
%     burn at a node are NOT explicitly targeted (must verify in GMAT,
%     FreeFlyer or a Cartesian J2 + finite-thrust propagator).
%   * No J2 differential argument-of-perigee correction, finite burn,
%     atmospheric drag, thrust limits, stationkeeping, collision safety,
%     or schedule-dependent customer orbital propagation.
%   * 1 km drift-altitude search grid is approximate, not a global optimum.
%   * Reported ΔV is an IDEAL LOWER-FIDELITY SCREENING COST, not a tank
%     sizing-ready propellant requirement or proven rendezvous trajectory.
%
% Outputs:
%   R.minDV_mps, R.feasible, R.bestTOF_days, R.bestDriftAltitude_km,
%   R.bestDriftDays, R.bestPoweredPlaneDV_mps, R.bestResidualRAAN_deg,
%   R.bestPhaseResidual_deg, R.bestCoastDV_mps, R.coverage_percent.
% Figures and CSV/MAT exported if cfg.exportFiles=true.
%
% No Optimization Toolbox required. MATLAB R2016b+.

if nargin<1, cfg=struct(); end
cfg = applyDefaults(cfg);
validateattributes(cfg.DVcap_mps,{'numeric'},{'scalar','positive'});
validateattributes(cfg.Tmax_days,{'numeric'},{'scalar','positive'});

mu = cfg.mu_km3s2; Re = cfg.Re_km; J2 = cfg.J2;
a0 = Re + cfg.hReference_km;
i0 = deg2rad(cfg.iReference_deg);
v0_mps = 1000*sqrt(mu/a0);
n0 = sqrt(mu/a0^3);  % [rad/s], circular Keplerian mean motion
OmegaDot0 = nodeRate(a0,0,i0,mu,Re,J2);

[raanGrid,phaseGrid] = meshgrid(cfg.raanDiff_deg,cfg.phaseDiff_deg);
nPhase = size(phaseGrid,1);
nRaan = size(raanGrid,2);

% Store only candidates that satisfy BOTH event limits.
minDV = inf(nPhase,nRaan);
bestTOF = nan(nPhase,nRaan);
bestHD = nan(nPhase,nRaan);
bestDrift = nan(nPhase,nRaan);
bestRAAN = nan(nPhase,nRaan);
bestPhase = nan(nPhase,nRaan);
bestPlaneDV = nan(nPhase,nRaan);
bestOrbitDV = nan(nPhase,nRaan);
bestBranch = nan(nPhase,nRaan);

maxTime_s = cfg.Tmax_days*86400;
raan0_rad = deg2rad(cfg.raanDiff_deg(:)');

% Exact reference coincidence is feasible without an orbit change.
idxPhaseZero = find(abs(cfg.phaseDiff_deg)<1e-12);
idxRaanZero = find(abs(cfg.raanDiff_deg)<1e-12);
if ~isempty(idxPhaseZero) && ~isempty(idxRaanZero)
    minDV(idxPhaseZero,idxRaanZero) = 0;
    bestTOF(idxPhaseZero,idxRaanZero) = 0;
    bestHD(idxPhaseZero,idxRaanZero) = cfg.hReference_km;
    bestDrift(idxPhaseZero,idxRaanZero) = 0;
    bestRAAN(idxPhaseZero,idxRaanZero) = 0;
    bestPhase(idxPhaseZero,idxRaanZero) = 0;
    bestPlaneDV(idxPhaseZero,idxRaanZero) = 0;
    bestOrbitDV(idxPhaseZero,idxRaanZero) = 0;
    bestBranch(idxPhaseZero,idxRaanZero) = 0;
end

hGrid = unique(sort([cfg.hDriftCandidates_km(:); cfg.hReference_km]));
numCandidateArcs=0;

for ia = 1:numel(hGrid)
    hD = hGrid(ia);
    aD = Re + hD;
    if abs(aD-a0)<1e-9, continue; end

    % Hohmann half-ellipse from a0 to aD, repeated in reverse on return.
    aT = 0.5*(a0+aD);
    eT = abs(aD-a0)/(aD+a0);
    tH_s = pi*sqrt(aT^3/mu);  % one Hohmann half-ellipse
    tTransfer_s = 2*tH_s;
    if tTransfer_s>=maxTime_s, continue; end

    vCirc0 = sqrt(mu/a0);
    vCircD = sqrt(mu/aD);
    vOnT_at0 = sqrt(mu*(2/a0-1/aT));
    vOnT_atD = sqrt(mu*(2/aD-1/aT));
    dvOrbit_mps = 2000*(abs(vOnT_at0-vCirc0) + ...
                         abs(vOnT_atD-vCircD));
    if dvOrbit_mps > cfg.DVcap_mps, continue; end

    nD = sqrt(mu/aD^3);
    omegaDotD = nodeRate(aD,0,i0,mu,Re,J2);
    omegaDotT = nodeRate(aT,eT,i0,mu,Re,J2);
    dN = nD-n0;
    if abs(dN)<1e-14, continue; end

    % The TWO half-ellipses advance the servicer by 2*pi in the
    % ideal planar geometry, while the client advances n0*(2*tH).
    phaseOffsetArc = 2*pi-n0*tTransfer_s;
    raanOffsetArc = (omegaDotT-OmegaDot0)*tTransfer_s;
    dOmegaRate = omegaDotD-OmegaDot0;
    tMaxDrift = maxTime_s-tTransfer_s;

    for ip=1:nPhase
        du0=deg2rad(cfg.phaseDiff_deg(ip));
        % Exact phase solution on branch k:
        % du0 + phaseOffsetArc + dN*tDrift = 2*pi*k.
        kA=(du0+phaseOffsetArc+dN*cfg.minDrift_s)/(2*pi);
        kB=(du0+phaseOffsetArc+dN*tMaxDrift)/(2*pi);
        kMin=ceil(min(kA,kB)-1e-12);
        kMax=floor(max(kA,kB)+1e-12);
        if kMin>kMax, continue; end

        kvals=kMin:kMax;
        tD=(2*pi*kvals-du0-phaseOffsetArc)/dN;
        mask=tD>=cfg.minDrift_s-1e-6 & tD<=tMaxDrift+1e-6;
        tD=tD(mask);
        kvals=kvals(mask);
        if isempty(tD), continue; end

        numCandidateArcs=numCandidateArcs+numel(tD);
        for kk=1:numel(tD)
            tEvent=tD(kk)+tTransfer_s;
            dOFinal=wrapPi(raan0_rad+raanOffsetArc+dOmegaRate*tD(kk));
            residualDeg=rad2deg(dOFinal);

            % We allow at most a SMALL residual powered plane correction.
            % Its actual burn must be scheduled at an appropriate node;
            % the node constraint is NOT propagated in this screening.
            permitted=abs(residualDeg) <= cfg.maxPoweredRAAN_deg+1e-10;
            if ~any(permitted), continue; end

            gamma=acos(clamp(cos(i0)^2+sin(i0)^2*cos(dOFinal),-1,1));
            dvPlane=2*v0_mps*sin(0.5*gamma);
            dvTotal=dvOrbit_mps+dvPlane;
            permitted = permitted & dvTotal <= cfg.DVcap_mps+1e-8;
            improved=permitted & (dvTotal<minDV(ip,:));
            if ~any(improved), continue; end

            minDV(ip,improved)=dvTotal(improved);
            bestTOF(ip,improved)=tEvent/86400;
            bestHD(ip,improved)=hD;
            bestDrift(ip,improved)=tD(kk)/86400;
            bestRAAN(ip,improved)=residualDeg(improved);
            bestPlaneDV(ip,improved)=dvPlane(improved);
            bestOrbitDV(ip,improved)=dvOrbit_mps;
            bestBranch(ip,improved)=kvals(kk);
            % Phase residual must be 0 in the secular model, modulo 2*pi.
            bestPhase(ip,improved)=rad2deg(wrapPi( ...
                du0+phaseOffsetArc+dN*tD(kk)));
        end
    end
end

feasible=isfinite(minDV);
minDV(~feasible)=NaN;
coverage=100*nnz(feasible)/numel(feasible);

R=struct();
R.cfg=cfg;
R.initialRAAN_deg=raanGrid;
R.initialPhase_deg=phaseGrid;
R.minDV_mps=minDV;
R.feasible=feasible;

R.bestTOF_days=bestTOF;
R.bestDriftAltitude_km=bestHD;
R.bestDriftDays=bestDrift;
R.bestResidualRAAN_deg=bestRAAN;
R.bestPhaseResidual_deg=bestPhase;

R.bestPoweredPlaneDV_mps=bestPlaneDV;
R.bestOrbitChangeDV_mps=bestOrbitDV;
R.bestPhaseBranch=bestBranch;
R.coverage_percent=coverage;
R.screenedCandidatePhaseBranches=numCandidateArcs;

fprintf('\nFABER JOINT RAAN–PHASE ACCESSIBILITY | SCREENING ONLY\n');
fprintf('Reference h = %.2f km, i = %.2f deg\n', ...
        cfg.hReference_km,cfg.iReference_deg);
fprintf('Limits: T <= %.1f days | DV <= %.1f m/s\n', ...
        cfg.Tmax_days,cfg.DVcap_mps);
fprintf('Residual powered RAAN <= %.2f deg\n', ...
        cfg.maxPoweredRAAN_deg);
fprintf('Grid accessibility: %d / %d = %.2f%%\n', ...
        nnz(feasible),numel(feasible),coverage);
fprintf('(Hypothetical relative geometry grid; NOT commercial coverage.)\n');

sampleRaan=[0 2 -2 5 -5]; samplePhase=[20 20 20 20 20];
fprintf('\nRepresentative cases [initial RAAN, initial phase] deg:\n');
for j=1:numel(sampleRaan)
    jr=find(abs(cfg.raanDiff_deg-sampleRaan(j))<1e-9,1);
    jp=find(abs(cfg.phaseDiff_deg-samplePhase(j))<1e-9,1);
    if isempty(jr) || isempty(jp), continue; end
    if feasible(jp,jr)
        fprintf('[%+g,%+g]: %.3f m/s, %.2f d, drift h=%.1f km, plane correction %.3f m/s\n', ...
           sampleRaan(j),samplePhase(j),minDV(jp,jr),bestTOF(jp,jr), ...
           bestHD(jp,jr),bestPlaneDV(jp,jr));
    else
        fprintf('[%+g,%+g]: NOT FEASIBLE under screening limits\n', ...
            sampleRaan(j),samplePhase(j));
    end
end

if cfg.makeFigures
    makeMaps(R);
end

if cfg.exportFiles
    save('FABER_RAAN_Phase_Joint30d_results.mat','R');
    [P,O]=ndgrid(cfg.phaseDiff_deg,cfg.raanDiff_deg);
    tbl=table(O(:),P(:),feasible(:),minDV(:),bestTOF(:), ...
        bestHD(:),bestDrift(:),bestOrbitDV(:),bestPlaneDV(:), ...
        bestRAAN(:),bestPhase(:),'VariableNames', ...
        {'RAAN0_deg','Phase0_deg','Feasible','DV_mps','TOF_days', ...
         'DriftAltitude_km','DriftDays','OrbitChangeDV_mps', ...
         'PlaneCorrectionDV_mps','PreCorrectionRAAN_deg','PhaseResidual_deg'});
    writetable(tbl,'FABER_RAAN_Phase_Joint30d_grid.csv');
end
end

%% ----- helpers -----
function cfg=applyDefaults(cfg)
D=struct();
D.mu_km3s2=398600.4418;
D.Re_km=6378.1363;
D.J2=1.08262668e-3;
D.hReference_km=480.1;
D.iReference_deg=51.8;
D.Tmax_days=180;
D.DVcap_mps=100;
% Fewer drift orbits permitted at low altitude due to drag and safety.
D.hDriftCandidates_km=380:1:700;
D.raanDiff_deg=-180:1:180;
% [-180, 180) avoids duplicate phase endpoint in coverage statistic.
D.phaseDiff_deg=-180:10:170;
% 0.5 deg maximum is an assumption to control simplified plane-correction.
D.maxPoweredRAAN_deg=0.5;
D.minDrift_s=0;
D.makeFigures=true;
D.exportFiles=true;
F=fieldnames(D);
for k=1:numel(F)
    if ~isfield(cfg,F{k}),cfg.(F{k})=D.(F{k});end
end
end

function rate=nodeRate(a,e,i,mu,Re,J2)
p=a*(1-e^2);
rate=-1.5*J2*sqrt(mu/a^3)*(Re/p)^2*cos(i);
end

function x=wrapPi(x)
x=atan2(sin(x),cos(x));
end

function x=clamp(x,lo,hi)
x=max(lo,min(hi,x));
end

function makeMaps(R)
c=R.cfg;
figure('Color','w','Name','FABER joint RAAN-phase feasibility');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile;
imagesc(c.raanDiff_deg,c.phaseDiff_deg,R.minDV_mps);
axis xy; xlabel('Initial ΔRAAN [deg]');
ylabel('Initial Δphase [deg]');
title('Minimum feasible event ΔV [m/s]');
caxis([0 c.DVcap_mps]); colorbar; grid on;
nexttile;
imagesc(c.raanDiff_deg,c.phaseDiff_deg,double(R.feasible));
axis xy; xlabel('Initial ΔRAAN [deg]');
ylabel('Initial Δphase [deg]');
title(sprintf('Accessible within %.0f d and %.0f m/s', ...
    c.Tmax_days,c.DVcap_mps));
colormap(gca,[0.9 0.9 0.9;0.15 0.55 0.7]);
caxis([0 1]);colorbar;grid on;
sgtitle(sprintf('FABER screening | h %.1f km | i %.1f° | 4-burn drift orbit', ...
   c.hReference_km,c.iReference_deg));

figure('Color','w','Name','FABER feasible drift choices');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile;
imagesc(c.raanDiff_deg,c.phaseDiff_deg,R.bestDriftAltitude_km);
axis xy; xlabel('Initial ΔRAAN [deg]');ylabel('Initial Δphase [deg]');
title('Chosen drift altitude [km]'); colorbar; grid on;
nexttile;
imagesc(c.raanDiff_deg,c.phaseDiff_deg,R.bestTOF_days);
axis xy; xlabel('Initial ΔRAAN [deg]');ylabel('Initial Δphase [deg]');
title('Selected total transfer time [days]'); colorbar;grid on;
sgtitle('Feasible candidate geometry and time');
end
