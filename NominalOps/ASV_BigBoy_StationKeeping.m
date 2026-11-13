function out = ASV_BigBoy_StationKeeping(cfg)
% ASV_BIGBOY_STATIONKEEPING  Preliminary FABER free-flyer parking trade.
%
% out = ASV_BigBoy_StationKeeping();       % illustrative defaults
% out = ASV_BigBoy_StationKeeping(cfg);    % user overrides of defaults
%
% Reference: circular host orbit; Big Boy initially parked at -V-bar
% (x,y,z) = (0,-distance,0) in host RTN/LVLH coordinates.
%
% Dynamics: nonlinear inertial 2-body + J2 + simplified corotating,
% exponential atmosphere drag for BOTH host and Big Boy.
% Controller: discrete, ideal velocity impulses at position/velocity
% deadband violations. This is an illustrative control law, NOT flight GNC.
%
% All dynamics use SI units. DO NOT interpret the estimated delta-V as a
% validated propulsion requirement. Atmosphere, Cd*A/m, parking safety,
% maneuver duration, and spacecraft design must be refined.

if nargin == 0 || isempty(cfg)
    cfg = defaultConfig();
else
    % Permit partial overrides of the illustrative baseline configuration.
    defaults = defaultConfig();
    fields = fieldnames(cfg);
    for q=1:numel(fields)
        defaults.(fields{q}) = cfg.(fields{q});
    end
    cfg = defaults;
end
validateConfig(cfg);
if cfg.saveCSV && ~exist(cfg.outputFolder,'dir')
    mkdir(cfg.outputFolder);
end

nCases = numel(cfg.distances_m)*numel(cfg.deadbands_m);
Distance_m = zeros(nCases,1);
Deadband_m = zeros(nCases,1);
DeltaV_mps = zeros(nCases,1);
Corrections = zeros(nCases,1);
MaxError_m = zeros(nCases,1);
MinRange_m = zeros(nCases,1);
MaxRange_m = zeros(nCases,1);
MaxRelSpeed_mps = zeros(nCases,1);
FiringTime_s = nan(nCases,1);
DutyCycle_pct = nan(nCases,1);
SafetyViolation = false(nCases,1);

caseResults = cell(nCases,1);
k = 0;
for i = 1:numel(cfg.distances_m)
    for j = 1:numel(cfg.deadbands_m)
        k = k+1;
        d = cfg.distances_m(i);
        deadband = cfg.deadbands_m(j);
        sim = runCase(cfg,d,deadband);
        caseResults{k} = sim;
        Distance_m(k) = d;
        Deadband_m(k) = deadband;
        DeltaV_mps(k) = sim.totalDV_mps;
        Corrections(k) = numel(sim.burnDV_mps);
        MaxError_m(k) = sim.maxPositionError_m;
        MinRange_m(k) = sim.minRange_m;
        MaxRange_m(k) = sim.maxRange_m;
        MaxRelSpeed_mps(k) = sim.maxRelativeSpeed_mps;
        FiringTime_s(k) = sim.equivalentFiringTime_s;
        DutyCycle_pct(k) = 100*sim.equivalentFiringTime_s/(cfg.missionDays*86400);
        SafetyViolation(k) = sim.safetyViolation;
    end
end

out.cfg = cfg;
out.summaryTable = table(Distance_m,Deadband_m,DeltaV_mps, ...
    Corrections,MaxError_m,MinRange_m,MaxRange_m,MaxRelSpeed_mps, ...
    FiringTime_s,DutyCycle_pct,SafetyViolation);
out.caseResults = caseResults;

disp(out.summaryTable);
if any(SafetyViolation)
    warning('At least one case violates the illustrative minimum range.');
end
if cfg.saveCSV
    writetable(out.summaryTable,fullfile(cfg.outputFolder,'BigBoy_SK_summary.csv'));
end
if cfg.makePlots
    figure('Color','w','Name','Big Boy parking trade: station-keeping Delta V');
    hold on; grid on;
    for j=1:numel(cfg.deadbands_m)
        mask = Deadband_m == cfg.deadbands_m(j);
        plot(Distance_m(mask),DeltaV_mps(mask),'-o','LineWidth',1.4, ...
            'DisplayName',sprintf('Deadband %.0f m',cfg.deadbands_m(j)));
    end
    xlabel('Big Boy nominal -V-bar distance [m]');
    ylabel('Estimated station-keeping Delta V [m/s]');
    title(sprintf('FABER Big Boy / %.1f days / preliminary drag + J2',cfg.missionDays));
    legend('Location','best');
    if cfg.saveCSV
        saveas(gcf,fullfile(cfg.outputFolder,'BigBoy_SK_trade.png'));
    end
end
end

function cfg = defaultConfig()
% Example only. These are not validated FABER specifications.
cfg.mu_m3s2 = 3.986004418e14;
cfg.Re_m = 6378136.3;
cfg.J2 = 1.08262668e-3;
cfg.omegaEarth_rads = 7.2921159e-5;
cfg.hostAltitude_m = 457250;   % change for other FABER orbit cases
cfg.inclination_deg = 40;
cfg.distances_m = [100 250 500 1000];
cfg.deadbands_m = [5 10 25];
cfg.missionDays = 7;            % short shakedown; test 30/90 later
cfg.step_s = 60;                % fixed-step RK4; test convergence
cfg.recordEvery_s = 3600;
cfg.hostMass_kg = 1000;        % PLACEHOLDER
cfg.hostCd = 2.2;              % PLACEHOLDER
cfg.hostArea_m2 = 8;           % PLACEHOLDER
cfg.bigBoyMass_kg = 250;       % PLACEHOLDER
cfg.bigBoyCd = 2.2;            % PLACEHOLDER
cfg.bigBoyArea_m2 = 2;         % PLACEHOLDER
cfg.rhoRef_kgm3 = 1e-12;       % PLACEHOLDER; not a forecast
cfg.rhoRefAltitude_m = 457250;
cfg.scaleHeight_m = 60000;     % PLACEHOLDER, local approximation
cfg.useJ2 = true;
cfg.useDrag = true;
cfg.velocityDeadband_mps = 0.02; % PLACEHOLDER
cfg.correctionTime_s = 1800;   % nominal time to reduce position error
cfg.maxImpulse_mps = 0.10;     % per command; PLACEHOLDER
cfg.minSafeRange_m = 50;       % PLACEHOLDER, not a certified KOZ
cfg.thrust_N = 0.1;            % PLACEHOLDER, for firing-time estimate
cfg.outputFolder = 'ASV_BigBoy_SK_results';
cfg.saveCSV = true;
cfg.makePlots = true;
end

function validateConfig(c)
positiveFields = {'mu_m3s2','Re_m','hostAltitude_m','step_s', ...
    'recordEvery_s','hostMass_kg','hostCd','hostArea_m2', ...
    'bigBoyMass_kg','bigBoyCd','bigBoyArea_m2', ...
    'rhoRef_kgm3','scaleHeight_m','correctionTime_s', ...
    'maxImpulse_mps','thrust_N','missionDays'};
for i=1:numel(positiveFields)
    q=c.(positiveFields{i});
    assert(isscalar(q) && isfinite(q) && q>0, ...
        '%s must be finite and positive.',positiveFields{i});
end
assert(all(isfinite(c.distances_m)) && all(c.distances_m>0));
assert(all(isfinite(c.deadbands_m)) && all(c.deadbands_m>0));
assert(all(c.distances_m>c.minSafeRange_m), ...
    'Nominal parking distance must exceed minimum safe range.');
assert(c.velocityDeadband_mps>0 && isfinite(c.velocityDeadband_mps));
assert(c.recordEvery_s>=c.step_s);
assert(islogical(c.useJ2) && islogical(c.useDrag));
end

function sim = runCase(cfg,distance_m,deadband_m)
% Initial host: circular equatorial longitude at ascending node, then
% rotated into the desired inclination. Does not model argument of latitude
% or local solar time variations in atmospheric density.
mu = cfg.mu_m3s2;
rOrbit = cfg.Re_m + cfg.hostAltitude_m;
vCircular = sqrt(mu/rOrbit);
inc = deg2rad(cfg.inclination_deg);
rH = [rOrbit;0;0];
vH = [0;vCircular*cos(inc);vCircular*sin(inc)];

[C,omegaLVLH] = frameRTN(rH,vH);
rRel = [0;-distance_m;0];
vRel = zeros(3,1);
rB = rH + C*rRel;
vB = vH + C*vRel + cross(omegaLVLH,C*rRel);
Y = [rH;vH;rB;vB];

T = cfg.missionDays*86400;
numSteps = ceil(T/cfg.step_s);
t = 0;
nextRecord = 0;
logs = zeros(0,6); % [time, x,y,z, position error, relative speed]
burnTimes = zeros(0,1);
burnDV = zeros(0,1);
totalDV = 0;
maxError = 0;
minRange = inf;
maxRange = 0;
maxSpeed = 0;
safetyViolation = false;

for it=0:numSteps
    [r,v] = relativeRTN(Y);
    e = r - [0;-distance_m;0];
    rnorm = norm(r);
    enorm = norm(e);
    vnorm = norm(v);
    maxError = max(maxError,enorm);
    minRange = min(minRange,rnorm);
    maxRange = max(maxRange,rnorm);
    maxSpeed = max(maxSpeed,vnorm);
    if rnorm < cfg.minSafeRange_m
        safetyViolation = true;
    end

    if t >= nextRecord-1e-8 || it==numSteps
        logs(end+1,:) = [t r' enorm vnorm]; %#ok<AGROW>
        nextRecord = nextRecord + cfg.recordEvery_s;
    end

    if it==numSteps || t >= T
        break;
    end

    % State-dependent correction, evaluated once at every integration step.
    % The command is an IDEAL impulsive relative velocity correction.
    if enorm > deadband_m || vnorm > cfg.velocityDeadband_mps
        vDesired = -e/cfg.correctionTime_s;
        dvRTN = vDesired - v;
        if norm(dvRTN) > cfg.maxImpulse_mps
            dvRTN = dvRTN*cfg.maxImpulse_mps/norm(dvRTN);
        end
        [C,~] = frameRTN(Y(1:3),Y(4:6));
        Y(10:12) = Y(10:12) + C*dvRTN;
        thisDV = norm(dvRTN);
        burnTimes(end+1,1) = t; %#ok<AGROW>
        burnDV(end+1,1) = thisDV; %#ok<AGROW>
        totalDV = totalDV + thisDV;
    end

    dt = min(cfg.step_s,T-t);
    k1 = dynamics(Y,cfg);
    k2 = dynamics(Y+dt*k1/2,cfg);
    k3 = dynamics(Y+dt*k2/2,cfg);
    k4 = dynamics(Y+dt*k3,cfg);
    Y = Y + dt*(k1+2*k2+2*k3+k4)/6;
    t = t+dt;
end

sim.distance_m = distance_m;
sim.deadband_m = deadband_m;
sim.totalDV_mps = totalDV;
sim.maxPositionError_m = maxError;
sim.minRange_m = minRange;
sim.maxRange_m = maxRange;
sim.maxRelativeSpeed_mps = maxSpeed;
sim.safetyViolation = safetyViolation;
sim.burnTime_s = burnTimes;
sim.burnDV_mps = burnDV;
sim.history = array2table(logs,'VariableNames', ...
    {'Time_s','X_m','Y_m','Z_m','PositionError_m','RelativeSpeed_mps'});
% Equivalent cumulative thrust on-time: m*sum(dV)/F, using constant
% Big Boy mass and ideal thrust direction. NOT a physical pulse schedule.
sim.equivalentFiringTime_s = cfg.bigBoyMass_kg*totalDV/cfg.thrust_N;
end

function dy = dynamics(Y,cfg)
rH=Y(1:3); vH=Y(4:6);
rB=Y(7:9); vB=Y(10:12);
aH=acceleration(rH,vH,cfg, ...
    cfg.hostCd*cfg.hostArea_m2/cfg.hostMass_kg);
aB=acceleration(rB,vB,cfg, ...
    cfg.bigBoyCd*cfg.bigBoyArea_m2/cfg.bigBoyMass_kg);
dy=[vH;aH;vB;aB];
end

function a = acceleration(r,v,cfg,CdAoverM)
rmag=norm(r);
mu=cfg.mu_m3s2;
a=-mu*r/rmag^3;
if cfg.useJ2
    z2 = (r(3)/rmag)^2;
    fac = 1.5*cfg.J2*mu*cfg.Re_m^2/rmag^5;
    a = a + fac*[r(1)*(5*z2-1); r(2)*(5*z2-1); r(3)*(5*z2-3)];
end
if cfg.useDrag
    alt=rmag-cfg.Re_m;
    rho=cfg.rhoRef_kgm3*exp(-(alt-cfg.rhoRefAltitude_m)/cfg.scaleHeight_m);
    vAtmos=cross([0;0;cfg.omegaEarth_rads],r);
    vAir=v-vAtmos;
    a=a-0.5*rho*CdAoverM*norm(vAir)*vAir;
end
end

function [C,omega] = frameRTN(r,v)
eR=r/norm(r);
h=cross(r,v);
eN=h/norm(h);
eT=cross(eN,eR);
C=[eR eT eN];
omega=h/norm(r)^2;
end

function [rRel,vRel] = relativeRTN(Y)
rH=Y(1:3); vH=Y(4:6);
rB=Y(7:9); vB=Y(10:12);
[C,omega]=frameRTN(rH,vH);
dr=rB-rH;
rRel=C'*dr;
vRel=C'*(vB-vH-cross(omega,dr));
end
