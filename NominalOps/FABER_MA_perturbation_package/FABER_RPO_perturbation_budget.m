function B = FABER_RPO_perturbation_budget(R,opts)
% FABER_RPO_PERTURBATION_BUDGET  Phase-A perturbation correction estimate.
%
% Uses the ORIGINAL analytically feasible nine-phase RPO results. Does not
% modify R, retime transfers, or chain terminal errors between phases.
% For each nominal transfer:
%   1) Evaluate deterministic differential J2, drag (and optional SRP)
%      acceleration along the nominal relative trajectory.
%   2) Integrate the LINEAR HCW perturbation error for entire phase.
%   3) At alpha*T prescribe one additional correction impulse so that the
%      original terminal position can be reached at the original TOF.
%   4) Adjust terminal impulse to restore original terminal velocity.
%   5) Report nonnegative delta-v ADDITIONAL to the original trajectory.
% For short HP holds: estimate passive perturbation drift; only include a
% prescribed within-hold restoration if the position/velocity tolerance is
% exceeded. No GNC controller or instantaneous position reset.
%
% Usage:
%   R = FABER_RPO_9phase_MC(cfg);
%   opts = struct('inclination_deg',98.12,'selection','selected');
%   B = FABER_RPO_perturbation_budget(R,opts);
%
% Requires FABER_MA_dynamics.m. All distances [m], times [s], DV [m/s].
% CAUTION: differential accelerations use a circular chief ephemeris;
% this is a SHORT-ARC model, not a long-term J2-consistent propagator.
if nargin<2, opts=struct(); end
opts=setdefault(opts,'inclination_deg',98.12);
opts=setdefault(opts,'raan_deg',0);
opts=setdefault(opts,'argumentLatitude_deg',0);
opts=setdefault(opts,'correctionFraction',0.5);
opts=setdefault(opts,'selection','selected'); % selected/pareto/all
opts=setdefault(opts,'maxCandidates',50);
opts=setdefault(opts,'holdPositionTolerance_m',0.5);
opts=setdefault(opts,'holdVelocityTolerance_mps',0.001);
opts=setdefault(opts,'maxCorrectionBurnFraction',0.05);
opts=setdefault(opts,'geometryTolerance_m',0.05);
opts=setdefault(opts,'speedTolerance_mps',0.005);
opts=setdefault(opts,'thrust_N',R.cfg.thrust_N*R.cfg.nActiveThrusters*R.cfg.thrustEfficiency);
opts=setdefault(opts,'mass_kg',R.cfg.massTaxi_kg);
opts=setdefault(opts,'Isp_s',R.cfg.Isp_s);
opts=setdefault(opts,'env',FABER_MA_dynamics('defaults'));
assert(opts.correctionFraction>0 && opts.correctionFraction<1);
assert(opts.thrust_N>0 && opts.mass_kg>0 && opts.Isp_s>0);
[N,M]=size(R.valid);
mask=false(N,M);
switch lower(char(opts.selection))
    case 'all'
        mask=R.valid;
    case 'pareto'
        assert(isfield(R,'scenarioPareto'),'Missing R.scenarioPareto');
        mask=R.valid & R.scenarioPareto;
    case 'selected'
        assert(isfield(R,'selectedArchitecture'),'Missing selectedArchitecture');
        for i=1:N
            m=R.selectedArchitecture(i);
            if isfinite(m) && m>=1 && m<=M && R.valid(i,m)
                mask(i,m)=true;
            end
        end
    otherwise
        error('Selection must be selected, pareto or all.');
end
[ii,mm]=find(mask);
if isfinite(opts.maxCandidates) && numel(ii)>opts.maxCandidates
    warning('Limiting perturbation analysis to %d of %d candidates.', ...
        opts.maxCandidates,numel(ii));
    ii=ii(1:opts.maxCandidates); mm=mm(1:opts.maxCandidates);
end
B=struct(); B.options=opts; B.attempted=false(N,M);
B.validCorrection=false(N,M); B.CorrectionDV_mps=nan(N,M);
B.HoldingDV_mps=nan(N,M); B.CorrectionPropellant_kg=nan(N,M);
B.details=cell(N,M); B.failureReason=strings(N,M);
rows=cell(0,10);
for q=1:numel(ii)
    i=ii(q); m=mm(q); rec=R.records{i,m};
    if isempty(rec), continue; end
    B.attempted(i,m)=true; good=true; reason="";
    dvt=0; dvh=0; maxHoldDrift=0; diag=cell(2,9);
    for leg=1:2
        if leg==1, S=rec.Client; else, S=rec.ATLANTE; end
        for k=1:5
            idx=2*k-1;
            z=S.transferDetails{k};
            [a,fail]=transferCorrection(z,S.trace{idx},opts,R.cfg.h_km,k,R.cfg);
            diag{leg,idx}=a;
            if strlength(fail)>0
                good=false; reason=sprintf('approach%d_phase%d_%s',leg,idx,fail);
                break;
            end
            dvt=dvt+a.extraDV_mps;
            if k<=4 && S.tofPhase(2*k)>0
                [h,fail]=holdingCorrection(z.targetState(:),S.tofPhase(2*k)*60,opts,R.cfg.h_km);
                diag{leg,idx+1}=h;
                dvh=dvh+h.extraDV_mps;
                maxHoldDrift=max(maxHoldDrift,h.predictedDrift_m);
                if strlength(fail)>0
                    good=false; reason=sprintf('approach%d_hold%d_%s',leg,k,fail);
                    break;
                end
            end
        end
        if ~good, break; end
    end
    B.details{i,m}=diag;
    B.validCorrection(i,m)=good; B.failureReason(i,m)=reason;
    B.CorrectionDV_mps(i,m)=dvt;
    B.HoldingDV_mps(i,m)=dvh;
    % Incremental propellant at configured reference mass; do not add this
    % to any other perturbation allocation without checking double counting.
    extraDV=dvt+dvh;
    prop=opts.mass_kg*(1-exp(-extraDV/(opts.Isp_s*9.80665)));
    B.CorrectionPropellant_kg(i,m)=prop;
    rows(end+1,:)={i,m,good,reason,R.CycleDV_mps(i,m), ...
        R.CycleTOF_min(i,m),dvt,dvh,extraDV,maxHoldDrift}; %#ok<AGROW>
end
B.table=cell2table(rows,'VariableNames',{'Trial','Architecture','CorrectionFeasible',...
    'Reason','NominalDV_mps','NominalTOF_min','TransferCorrectionDV_mps',...
    'HoldingCorrectionDV_mps','ExtraDV_mps','MaxPredictedHoldDrift_m'});
if ~isempty(B.table)
    B.table.ExtraPropellant_kg=opts.mass_kg*(1-exp(-B.table.ExtraDV_mps/(opts.Isp_s*9.80665)));
    disp(B.table);
end
end

function [o,fail]=transferCorrection(z,trace,opts,h_km,k,cfg)
fail="";
T=z.TOF_s; n=sqrt(opts.env.mu_m3s2/(opts.env.Re_m+1000*h_km)^3);
% The stored trace is a uniformly sampled nominal trajectory. In the
% continuous case it is the nominal forced quintic reference; for impulsive
% phases it is the coast state AFTER the departure impulse.
tGrid=linspace(0,T,size(trace,2));
forced=@(t,X) errorRHS(t,X,trace,tGrid,n,opts,h_km);
o=struct('extraDV_mps',0,'predictedPosError_m',NaN, ...
    'predictedVelError_mps',NaN,'midcourseDV_mps',0,'terminalDeltaDV_mps',0, ...
    'correctionTOF_s',T,'midcourseTime_s',opts.correctionFraction*T);
try
    checkTimes=linspace(0,T,101);
    [tt,Y]=ode45(forced,checkTimes,zeros(6,1),odeset('RelTol',1e-8,'AbsTol',1e-10));
catch ME
    fail="integration_"+string(ME.identifier); return;
end
q=Y(end,:)';
o.predictedPosError_m=norm(q(1:3));
o.predictedVelError_mps=norm(q(4:6));
tr=(1-opts.correctionFraction)*T;
P=FABER_MA_dynamics('phi',n,tr);
if rcond(P(1:3,4:6))<1e-10, fail="singular_remaining_transfer";return;end
dvc=-P(1:3,4:6)\q(1:3);
dvterm=-(q(4:6)+P(4:6,4:6)*dvc);
% Nominal terminal impulse must be adjusted, not replaced.
if z.isImpulsive
    orig=z.burn2(:);
    increment=norm(dvc)+norm(orig+dvterm)-norm(orig);
else
    increment=norm(dvc)+norm(dvterm); % additional to nominal continuous control
end
o.extraDV_mps=max(0,increment); % allocation is nonnegative
o.midcourseDV_mps=norm(dvc);
o.terminalDeltaDV_mps=norm(dvterm);
% Feasibility of corrections, not just delta-v cost.
ac=opts.thrust_N/opts.mass_kg;
tBurn=(norm(dvc)+norm(dvterm))/ac;
o.estimatedCorrectionBurnTime_s=tBurn;
if tBurn/T>opts.maxCorrectionBurnFraction
    fail="correction_burn_fraction";
end
% Check linearized post-correction terminal position and velocity.
res=q+P(:,4:6)*dvc+[zeros(3,1);dvterm];
o.linearizedResidual_m=norm(res(1:3));
o.linearizedResidual_mps=norm(res(4:6));
% Verify that the FIRST-ORDER corrected path does not violate the original
% phase-specific geometry or speed constraints. This is a geometric screen,
% not a guaranteed nonlinear finite-burn path.
switch k
    case {1,2}, minAllowed=cfg.AZ_m;
    case 3, minAllowed=cfg.KOZ_m;
    otherwise, minAllowed=norm(z.endState(1:3));
end
maxSpeed=inf;
switch k
    case 1, maxSpeed=cfg.maxSpeed_homing1_mps;
    case 2, maxSpeed=cfg.maxSpeed_homing2_mps;
    case 3, maxSpeed=cfg.maxSpeed_closing1_mps;
    case 4, maxSpeed=cfg.maxSpeed_closing2_mps;
    case 5, maxSpeed=cfg.maxSpeed_final_mps;
end
minR=inf; maxV=0;
for j=1:numel(tt)
    Ynom=interp1(tGrid,trace',tt(j),'linear')';
    eState=Y(j,:)';
    if tt(j)>=o.midcourseTime_s
        eState=eState+FABER_MA_dynamics('phi',n,tt(j)-o.midcourseTime_s)*[zeros(3,1);dvc];
    end
    actual=Ynom+eState;
    minR=min(minR,norm(actual(1:3)));
    maxV=max(maxV,norm(actual(4:6)));
end
o.minCorrectedRange_m=minR;
o.maxCorrectedSpeed_mps=maxV;
if minR<minAllowed-opts.geometryTolerance_m
    fail="corrected_geometry";
elseif maxV>maxSpeed+opts.speedTolerance_mps
    fail="corrected_speed";
end
end

function [o,fail]=holdingCorrection(holdState,T,opts,h_km)
fail=""; n=sqrt(opts.env.mu_m3s2/(opts.env.Re_m+1000*h_km)^3);
o=struct('extraDV_mps',0,'predictedDrift_m',0,'predictedVelocityDrift_mps',0, ...
    'correctionRequired',false,'estimatedCorrectionBurnTime_s',0);
trace=repmat(holdState(:),1,2); tt=[0 T];
f=@(t,X) errorRHS(t,X,trace,tt,n,opts,h_km);
[~,Y]=ode45(f,[0 T],zeros(6,1),odeset('RelTol',1e-8,'AbsTol',1e-10));
q=Y(end,:)';
o.predictedDrift_m=norm(q(1:3));
o.predictedVelocityDrift_mps=norm(q(4:6));
if o.predictedDrift_m<=opts.holdPositionTolerance_m && ...
        o.predictedVelocityDrift_mps<=opts.holdVelocityTolerance_mps
    return;
end
% Prescribed correction at mid-hold, completing by its scheduled end.
% This estimates cost, it is NOT a closed-loop deadband controller.
o.correctionRequired=true;
rem=T/2;
P=FABER_MA_dynamics('phi',n,rem);
if rcond(P(1:3,4:6))<1e-10,fail="singular_hold";return;end
dv1=-P(1:3,4:6)\q(1:3);
dv2=-(q(4:6)+P(4:6,4:6)*dv1);
o.extraDV_mps=norm(dv1)+norm(dv2);
o.estimatedCorrectionBurnTime_s=o.extraDV_mps*opts.mass_kg/opts.thrust_N;
if o.estimatedCorrectionBurnTime_s/T>opts.maxCorrectionBurnFraction
    fail="hold_correction_burn_fraction";
end
end

function dx=errorRHS(t,x,trace,grid,n,opts,h_km)
% Time-interpolate 6xN stored nominal states without relying on handles
% serialized in the Monte Carlo records.
ref=interp1(grid,trace',min(grid(end),max(grid(1),t)),'linear')';
a=FABER_MA_dynamics('differential',t,ref,opts.env,h_km,...
    opts.inclination_deg,opts.raan_deg,opts.argumentLatitude_deg);
dx=[x(4:6);3*n^2*x(1)+2*n*x(5);-2*n*x(4);-n^2*x(3)];
dx(4:6)=dx(4:6)+a;
end

function s=setdefault(s,key,value)
if ~isfield(s,key)||isempty(s.(key)),s.(key)=value;end
end
