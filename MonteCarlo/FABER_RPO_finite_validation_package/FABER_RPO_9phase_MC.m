function R = FABER_RPO_9phase_MC(cfg)
% FABER_RPO_9PHASE_MC
% Nine-phase, four-holding-point RPO Monte Carlo.
%
% Compares 36 combinations with default gate modes:
% Homing 2: GeneralHCW, Vbar, Rbar
% Closing 1: GeneralHCW, Vbar, Rbar
% Closing 2: GeneralHCW, Continuous
%
% All GO decisions assumed true.
% No lateral corridor-width restriction.
%
% Requires:
%   analytical_manos.m
%   rpo_continuous.m

%checkConfig(cfg);

rng(cfg.seed,'twister');
validateattributes(cfg.thrust_N,{'numeric'},{'scalar','positive','finite'});
validateattributes(cfg.nActiveThrusters,{'numeric'},{'scalar','integer','positive'});
validateattributes(cfg.thrustEfficiency,{'numeric'},{'scalar','>',0,'<=',1});
validateattributes(cfg.maxBurnFraction,{'numeric'},{'scalar','>',0,'<=',1});
if ~isempty(cfg.thrusterDirections_LVLH)
    assert(size(cfg.thrusterDirections_LVLH,1)==3, ...
        'thrusterDirections_LVLH must be 3xN');
    assert(numel(cfg.thrusterForces_N)==size(cfg.thrusterDirections_LVLH,2), ...
        'One force is required for every directional thruster');
    assert(all(cfg.thrusterForces_N>0), 'Thruster forces must be positive');
    assert(exist('linprog','file')==2, ...
        'Directional thruster allocation requires Optimization Toolbox linprog');
end

validateattributes(cfg.massTaxi_kg,{'numeric'},{'scalar','positive','finite'});
validateattributes(cfg.Isp_s,{'numeric'},{'scalar','positive','finite'});
assert(cfg.shellRadius_m(1)>0 && cfg.shellRadius_m(2)>=cfg.shellRadius_m(1));
assert(ismember(lower(string(cfg.initialPositionMode)), ...
    ["alongtrack","shell","volume"]));
assert(cfg.maxBurnFraction>0 && cfg.maxBurnFraction<=1);
assert(cfg.maxPropellantCycle_kg>=0);

% mean motion
n = sqrt(cfg.mu_km3s2 / ...
    (cfg.Re_km+cfg.h_km)^3);

%% 1. GENERATE ALL METHOD COMBINATIONS
% homing 2, closing 1, closing 2, "gate mode" for HP4
% homing 1, final approach are fixed based on dynamics constraints
[A,B,C,G] = ndgrid( ...
    1:numel(cfg.methodsH2), ...
    1:numel(cfg.methodsC1), ...
    1:numel(cfg.methodsC2), ...
    1:numel(cfg.gateModes));

choices = [A(:),B(:),C(:),G(:)];

M = size(choices,1);
N = cfg.N;

methodH2 = strings(M,1);
methodC1 = strings(M,1);
methodC2 = strings(M,1);
gateMode = strings(M,1);

for m = 1:M
    methodH2(m) = string( ...
        cfg.methodsH2{choices(m,1)});

    methodC1(m) = string( ...
        cfg.methodsC1{choices(m,2)});

    methodC2(m) = string( ...
        cfg.methodsC2{choices(m,3)});

    gateMode(m) = string( ...
        cfg.gateModes{choices(m,4)});
end

%% 2. PREALLOCATE MONTE CARLO ARRAYS

% Dimensions:
%   N attempts x M architectures x 2 approaches

DV = nan(N,M,2);
TOF = nan(N,M,2);
FuelCycle_kg = nan(N,M);
BurnRatioCycle = nan(N,M);
BurnTimeCycle_s = nan(N,M);
BurnTimePhase_s = nan(N,M,2,9);
BurnFractionPhase = nan(N,M,2,9);
FuelPhase_kg = nan(N,M,2,9);

% Nine phases for each approach
DVphase = nan(N,M,2,9);
TOFphase = nan(N,M,2,9);

reason = strings(N,M);
records = cell(N,M);

%% 3. PAIRED MONTE CARLO

for i = 1:N

    % Same initial states used across all method combinations
    xc = sampleState(cfg.x0_client,cfg.dx0_client,cfg);
    xh = sampleState(cfg.x0_host,cfg.dx0_host,cfg);

    % Same sampled design used across all methods
    dc = sampleDesign(cfg,xc);
    dh = sampleDesign(cfg,xh);

    for m = 1:M

        which = { ...
            char(methodH2(m)), ...
            char(methodC1(m)), ...
            char(methodC2(m))};

        %% TAXI -> CLIENT (UNLOADED)

        cfgClient = cfg;
        cfgClient.massTaxi_kg = cfg.massTaxi_kg;

        [client,why] = simulateApproach( ...
            cfgClient,xc,dc,n,which,char(gateMode(m)));

        if ~isempty(why)
            reason(i,m) = "Client: " + string(why);
            continue;
        end

        %% PAYLOAD CAPTURE

        % Taxi has consumed fuel during client approach.
        massAfterClient_kg = ...
            cfgClient.massTaxi_kg - client.propellant_kg;

        if cfg.enablePayloadMassChange
            capturedMass_kg = cfg.payloadMass_kg;
        else
            capturedMass_kg = 0;
        end

        %% TAXI -> ATLANTE (LOADED)

        cfgHost = cfg;

        cfgHost.massTaxi_kg = ...
            massAfterClient_kg + capturedMass_kg;

        [host,why] = simulateApproach( ...
            cfgHost,xh,dh,n,which,char(gateMode(m)));

        if ~isempty(why)
            reason(i,m) = "ATLANTE: " + string(why);
            continue;
        end

        %% COMPLETE RETRIEVAL CYCLE

        % Propellant from both approaches, with the mass change
        % represented between the two sequences.
        cycleFuel = ...
            client.propellant_kg + host.propellant_kg;

        if cfg.enablePropulsionConstraints && ...
                cycleFuel > cfg.maxPropellantCycle_kg

            reason(i,m) = "Cycle propellant allocation";
            continue;
        end

        if cfg.enablePropulsionConstraints && ...
                client.dv + host.dv > cfg.maxDV_perCycle_mps

            reason(i,m) = "Cycle DV limit";
            continue;
        end

        FuelCycle_kg(i,m) = cycleFuel;
        BurnTimeCycle_s(i,m) = client.totalBurnTime_s + host.totalBurnTime_s;
        % Ratio is the largest phase-level impulsive-burn fraction,
        % not the duration of the entire retrieval cycle.
        BurnRatioCycle(i,m) = max([client.burnFractionPhase,host.burnFractionPhase]);
        BurnTimePhase_s(i,m,1,:) = reshape(client.burnTimePhase_s,1,1,1,9);
        BurnTimePhase_s(i,m,2,:) = reshape(host.burnTimePhase_s,1,1,1,9);
        BurnFractionPhase(i,m,1,:) = reshape(client.burnFractionPhase,1,1,1,9);
        BurnFractionPhase(i,m,2,:) = reshape(host.burnFractionPhase,1,1,1,9);
        FuelPhase_kg(i,m,1,:) = reshape(client.fuelPhase_kg,1,1,1,9);
        FuelPhase_kg(i,m,2,:) = reshape(host.fuelPhase_kg,1,1,1,9);
        DV(i,m,1) = client.dv;
        DV(i,m,2) = host.dv;

        TOF(i,m,1) = client.tof;
        TOF(i,m,2) = host.tof;

        DVphase(i,m,1,:) = reshape( ...
            client.dvPhase,1,1,1,9);
        DVphase(i,m,2,:) = reshape( ...
            host.dvPhase,1,1,1,9);

        TOFphase(i,m,1,:) = reshape( ...
            client.tofPhase,1,1,1,9);
        TOFphase(i,m,2,:) = reshape( ...
            host.tofPhase,1,1,1,9);

        records{i,m} = struct( ...
            'Client',client, ...
            'ATLANTE',host, ...
            'designClient',dc, ...
            'designATLANTE',dh);
    end
end

%% 4. COMPLETE-CYCLE RESULTS

cycleDV = sum(DV,3);
cycleTOF = sum(TOF,3);

valid = isfinite(cycleDV);
common = all(valid,2);

Accepted = sum(valid,1)';
Acceptance = Accepted/N;

P95ClientDV_mps = nan(M,1);
P95ATLANTEDV_mps = nan(M,1);

P95CycleDV_mps = nan(M,1);
MedianCycleDV_mps = nan(M,1);
P95CycleTOF_min = nan(M,1);
P95FuelCycle_kg = nan(M,1);
P95BurnTimeCycle_s = nan(M,1);
P95MaxPhaseBurnFraction = nan(M,1);

CommonP95CycleDV_mps = nan(M,1);

for m = 1:M

    ok = valid(:,m);

    if any(ok)

        P95ClientDV_mps(m) = p95(DV(ok,m,1));

        P95ATLANTEDV_mps(m) = ...
            p95(DV(ok,m,2));

        P95CycleDV_mps(m) = ...
            p95(cycleDV(ok,m));

        MedianCycleDV_mps(m) = ...
            median(cycleDV(ok,m));

        P95CycleTOF_min(m) = ...
            p95(cycleTOF(ok,m));
        P95FuelCycle_kg(m) = p95(FuelCycle_kg(ok,m));
        P95BurnTimeCycle_s(m) = p95(BurnTimeCycle_s(ok,m));
        P95MaxPhaseBurnFraction(m) = p95(BurnRatioCycle(ok,m));
    end

    if any(common)
        CommonP95CycleDV_mps(m) = ...
            p95(cycleDV(common,m));
    end
end

%% 5. OUTPUT STRUCTURE

R.summary = table( ...
    methodH2,methodC1,methodC2,gateMode, ...
    Accepted,Acceptance, ...
    P95ClientDV_mps,P95ATLANTEDV_mps, ...
    MedianCycleDV_mps,P95CycleDV_mps, ...
    P95CycleTOF_min,CommonP95CycleDV_mps, ...
    P95FuelCycle_kg,P95BurnTimeCycle_s,P95MaxPhaseBurnFraction);

% Discrete objective selection is performed within each scenario:
% all candidate architectures share that scenario's sampled initial states.
[selectedIdx,scenarioPareto] = selectCandidates( ...
    cfg,cycleDV,cycleTOF,valid);
R.selectedArchitecture = selectedIdx;
R.scenarioPareto = scenarioPareto;
R.selectedDV_mps = nan(N,1);
R.selectedTOF_min = nan(N,1);
for i = 1:N
    if isfinite(selectedIdx(i))
        R.selectedDV_mps(i) = cycleDV(i,selectedIdx(i));
        R.selectedTOF_min(i) = cycleTOF(i,selectedIdx(i));
    end
end

R.cfg = cfg;
R.methods = [ ...
    methodH2 methodC1 methodC2 gateMode];

R.DV_mps = DV;
R.TOF_min = TOF;

R.DV_phase_mps = DVphase;
R.TOF_phase_min = TOFphase;

R.CycleDV_mps = cycleDV;
R.CycleTOF_min = cycleTOF;
R.FuelCycle_kg = FuelCycle_kg;
R.BurnTimeCycle_s = BurnTimeCycle_s;
R.MaxPhaseBurnFraction = BurnRatioCycle;
R.BurnTimePhase_s = BurnTimePhase_s;
R.BurnFractionPhase = BurnFractionPhase;
R.FuelPhase_kg = FuelPhase_kg;

R.valid = valid;
R.commonFeasible = common;

R.rejectionReason = reason;
R.records = records;

R.phaseNames = { ...
    'Homing 1','Holding 1', ...
    'Homing 2','Holding 2', ...
    'Closing 1','Holding 3', ...
    'Closing 2','HP4 gate', ...
    'Final approach'};

fprintf('\nFABER NINE-PHASE MONTE CARLO\n');
fprintf('Orbit %.2f km | Period %.2f min\n', ...
    cfg.h_km,2*pi/n/60);

fprintf(['Attempts: %d | Architectures: %d | ' ...
         'Common feasible: %d\n'], ...
    N,M,sum(common));

disp(R.summary);

%% 6. SAVE RESULTS

if cfg.saveCSV

    if ~exist(cfg.outputFolder,'dir')
        mkdir(cfg.outputFolder);
    end

    writetable(R.summary, ...
        fullfile(cfg.outputFolder, ...
        'nine_phase_architectures.csv'));

    T = table();

    for m = 1:M

rows = table( ...
    (1:N)', ...
    repmat(m,N,1), ...
    repmat(methodH2(m),N,1), ...
    repmat(methodC1(m),N,1), ...
    repmat(methodC2(m),N,1), ...
    repmat(gateMode(m),N,1), ...
    valid(:,m), ...
    cycleDV(:,m), ...
    cycleTOF(:,m), ...
    FuelCycle_kg(:,m), BurnTimeCycle_s(:,m), BurnRatioCycle(:,m), ...
    reason(:,m), ...
    'VariableNames',{ ...
    'Attempt','Architecture', ...
    'Homing2','Closing1','Closing2', ...
    'HP4Gate', ...
    'Feasible','CycleDV_mps', ...
    'CycleTOF_min','FuelCycle_kg','BurnTimeCycle_s', ...
    'MaxPhaseBurnFraction','Reason'});

        T = [T;rows]; %#ok<AGROW>
    end

    writetable(T, ...
        fullfile(cfg.outputFolder, ...
        'nine_phase_trials.csv'));
end

%% 7. PLOT DELTA V VS TOF

if cfg.makePlots

    figure('Color','w');
    hold on;
    grid on;

    for m = 1:M
        ok = valid(:,m);
        scatter(cycleTOF(ok,m), ...
            cycleDV(ok,m),10,'filled');
    end

    xlabel('Complete cycle TOF [min]');
    ylabel('Complete cycle Delta-V [m/s]');

    title('FABER nine-phase maneuver comparison');
end

end

%% ============================================================
% SAMPLE DESIGN VARIABLES
% =============================================================

function d = sampleDesign(cfg,x0)

d.hp = [ ...
    draw(cfg.hp1_m), ...
    draw(cfg.hp2_m), ...
    draw(cfg.hp3_m), ...
    draw(cfg.hp4_m)];

d.final = draw(cfg.final_m);

d.hold = [ ...
    draw(cfg.hold1_min), ...
    draw(cfg.hold2_min), ...
    draw(cfg.hold3_min), ...
    draw(cfg.hold4_min)];

% Transfer durations [H1 H2 C1 C2 Final], seconds
d.T = zeros(1,5);

D = [ ...
    norm(x0(1:3)-[0;-d.hp(1);0]), ...
    diff(-d.hp), ...
    d.hp(4)-d.final];

D = abs(D);

bounds = { ...
    cfg.homing1_min, ...
    cfg.homing2_min, ...
    cfg.closing1_min, ...
    cfg.closing2_min, ...
    cfg.final_min};

speeds = [ ...
    cfg.maxSpeed_homing1_mps, ...
    cfg.maxSpeed_homing2_mps, ...
    cfg.maxSpeed_closing1_mps, ...
    cfg.maxSpeed_closing2_mps, ...
    cfg.maxSpeed_final_mps];

for k = 1:5

    % Necessary average-speed time bound.
    % This is valid regardless of endpoint velocities.
    lo = max([ ...
        60*bounds{k}(1), ...
        cfg.minTransferTime_s, ...
        D(k)/speeds(k)]);

    hi = 60*bounds{k}(2);

    if lo>hi
        d.valid = false;
        return;
    end

    d.T(k) = draw([lo hi]);
end

d.valid = true;

end

%% ============================================================
% SIMULATE ONE COMPLETE NINE-PHASE APPROACH
% =============================================================

function [s,reason] = simulateApproach( ...
    cfg,x0,d,n,methods,gateMode)

s = [];
reason = '';

if ~d.valid
    reason = 'time';
    return;
end

hp = d.hp;

%% 1. PRELIMINARY WAYPOINT GEOMETRY

if ~(norm(x0(1:3))>hp(1) && ...
     hp(1)>hp(2) && ...
     hp(2)>cfg.AZ_m && ...
     cfg.AZ_m>hp(3) && ...
     hp(3)>cfg.KOZ_m && ...
     cfg.KOZ_m>hp(4) && ...
     hp(4)>d.final && d.final>0)

    reason = 'geometry';
    return;
end

%% 2. STATIONARY NEGATIVE V-BAR ENDPOINTS
X = cell(1,5);

% HP1, HP2, HP3: stationary V-bar endpoints
for k = 1:3
    X{k} = [0;-hp(k);0;0;0;0];
end

% HP4: selectable full-stop or nonzero-velocity gate
switch lower(gateMode)

    case 'fullstop'
        vyGate = 0;

    case 'flythrough'
        vyGate = cfg.gateSpeed_mps;

    otherwise
        error('Unknown HP4 gate mode: %s',gateMode);
end

% Positive along-track velocity moves toward the target.
X{4} = [0;-hp(4);0;0;vyGate;0];

% Capture standoff: arrive at zero relative velocity.
X{5} = [0;-d.final;0;0;0;0];

N = max(101,cfg.nTrajectorySamples);

%% 3. DEFINE TRANSFERS

starts = {x0,X{1},X{2},X{3},X{4}};
ends = {X{1},X{2},X{3},X{4},X{5}};

family = { ...
    'GeneralHCW', ...
    methods{1}, ...
    methods{2}, ...
    methods{3}, ...
    'Continuous'};

maxSpeed = [ ...
    cfg.maxSpeed_homing1_mps, ...
    cfg.maxSpeed_homing2_mps, ...
    cfg.maxSpeed_closing1_mps, ...
    cfg.maxSpeed_closing2_mps, ...
    cfg.maxSpeed_final_mps];

windows = { ...
    cfg.homing1_min, ...
    cfg.homing2_min, ...
    cfg.closing1_min, ...
    cfg.closing2_min, ...
    cfg.final_min};

% Transfer phases in nine-phase array:
phaseIndices = [1 3 5 7 9];

traces = cell(1,9);
transferDetails = cell(1,5);
dvPhase = zeros(1,9);
tofPhase = zeros(1,9);
zoneEvents = cell(1,5);
burnTimePhase = zeros(1,9);
burnFractionPhase = zeros(1,9);
fuelDVphase = zeros(1,9);
%% 4. PROPAGATE AND VALIDATE EACH TRANSFER

for k = 1:5

    [z,why] = transfer( ...
        starts{k},ends{k},d.T(k), ...
        n,N,family{k});

    if ~isempty(why)
        reason = sprintf('phase%d_%s', ...
            phaseIndices(k),why);
        return;
    end

    % HARD TIME CONSTRAINT

    if z.TOF_s/60 < windows{k}(1)-1e-9 || ...
       z.TOF_s/60 > windows{k}(2)+1e-9

        reason = sprintf('phase%d_time', ...
            phaseIndices(k));
        return;
    end

    % For closing legs the user's speed bound may refer to
    % rho-dot (inward closing speed), not total LVLH speed.
    if k==3 || k==4
        mode = lower(string(cfg.closingSpeedMode));
        if mode=="radial"
            speedForLimit = z.maxClosingSpeed;
        elseif mode=="both"
            speedForLimit = max(z.maxSpeed,z.maxClosingSpeed);
        elseif mode=="total"
            speedForLimit = z.maxSpeed;
        else
            error('Unknown closingSpeedMode: %s',cfg.closingSpeedMode);
        end
    else
        speedForLimit = z.maxSpeed;
    end

    if speedForLimit > maxSpeed(k)+1e-9

        reason = sprintf('phase%d_speed', ...
            phaseIndices(k));
        return;
    end

    % OPTIONAL ACCELERATION LIMIT

    if cfg.enablePropulsionConstraints && ...
        z.peakAcceleration > cfg.maxPeakAccel_mps2

        reason = sprintf('phase%d_acceleration', ...
            phaseIndices(k));
        return;
    end

    %% PHASE-SPECIFIC PROPULSION SCREENING

    % Five transfers correspond to nine-phase indices:
    % H1=1, H2=3, C1=5, C2=7, Final=9
    phaseNumber = phaseIndices(k);

    phaseCfg = FABER_RPO_phase_propulsion( ...
        cfg,phaseNumber);

    if cfg.enablePropulsionConstraints

        % Thrust magnitude and direction
        % Burn duration and 5% rule
        % Continuous-control thrust feasibility
        [burnOK,burnWhy,burnInfo] = ...
            screenTransferBurns(z,phaseCfg);

    else

        % Original unconstrained propulsion case
        burnOK = true;
        burnWhy = '';

        burnInfo = struct( ...
            'totalBurn_s',NaN, ...
            'fraction',NaN);

    end

    if ~burnOK
        reason = sprintf('phase%d_%s', ...
            phaseNumber,burnWhy);
        return;
    end

    %% BRAKING DISTANCE

    if cfg.enablePropulsionConstraints && ...
            cfg.checkBrakingDistance

        [brakeOK,brakeReason] = ...
            checkBraking(z,ends{k},phaseCfg);

        if ~brakeOK
            reason = sprintf('phase%d_%s', ...
                phaseNumber,brakeReason);
            return;
        end
    end

    %% EVENT-BASED AZ / KOZ VALIDATION - NEW

    [why,events] = verifyZoneEvents( ...
        z,k,cfg,d.final);

    zoneEvents{k} = events;

    if ~isempty(why)
        reason = sprintf('phase%d_%s', ...
            phaseIndices(k),why);
        return;
    end

    j = phaseIndices(k);

    traces{j} = z.trajectory;
    % Essential metadata for optional finite-thrust verification.
    transferDetails{k} = struct( ...
        'startState',z.startState, ...
        'targetState',z.endState, ...
        'TOF_s',z.TOF_s, ...
        'method',string(family{k}), ...
        'isImpulsive',z.isImpulsive, ...
        'burn1',z.burn1, ...
        'burn2',z.burn2, ...
        'idealDV_mps',z.dvTotal);
    dvPhase(j) = z.dvTotal;
    tofPhase(j) = z.TOF_s/60;
    burnTimePhase(j) = burnInfo.totalBurn_s;
    burnFractionPhase(j) = burnInfo.fraction;
    fuelDVphase(j) = z.dvTotal;
end

%% 5. HOLDING POINTS AND GO DECISIONS

% HP1, HP2, HP3 and HP4.
% All decisions assumed GO.
%
% HP4 defaults to no dwell, but still has
% zero relative velocity at arrival.

for k = 1:4

    j = 2*k;

    location = X{k};

    traces{j} = repmat(location,1,2);

    tofPhase(j) = d.hold(k);

    dvPhase(j) = d.hold(k)* ...
        cfg.holdDVrate_mps_per_min(k);
    fuelDVphase(j) = dvPhase(j);
end

%% 6. RETURN RESULTS

s = struct( ...
    'dv',sum(dvPhase), ...
    'tof',sum(tofPhase), ...
    'dvPhase',dvPhase, ...
    'tofPhase',tofPhase, ...
    'trace',{traces}, ...
    'transferDetails',{transferDetails}, ...
    'GO',true, ...
    'methods',{methods}, ...
    'holdingPoints_m',hp);

%% PHASE-SPECIFIC PROPELLANT ESTIMATE

remainingMass = cfg.massTaxi_kg;
fuelPhase = zeros(1,9);

for jj = 1:9

    phaseCfg = FABER_RPO_phase_propulsion(cfg,jj);

    nextMass = remainingMass * exp( ...
        -fuelDVphase(jj) / ...
        (9.80665*phaseCfg.Isp_s));

    fuelPhase(jj) = remainingMass-nextMass;

    remainingMass = nextMass;

end

s.propellant_kg = sum(fuelPhase);
s.fuelPhase_kg = fuelPhase;
s.burnTimePhase_s = burnTimePhase;
s.burnFractionPhase = burnFractionPhase;
s.totalBurnTime_s = sum(burnTimePhase);

s.gateMode = string(gateMode);
s.HP4_velocity_mps = X{4}(4:6);
s.zoneEvents = zoneEvents;

end

%% ============================================================
% TRANSFER METHOD DISPATCHER
% =============================================================

function [z,reason] = transfer(x0,xf,T,n,N,method)
% TRANSFER
% General HCW, analytical V-bar/R-bar, or continuous
% quintic HCW with arbitrary endpoint velocities.
%
% NEW: stateAt(t) allows accurate crossing-time
% refinement rather than using linear interpolation.

z = [];
reason = '';
burn1 = zeros(3,1);
burn2 = zeros(3,1);
isImpulsive = false;

if ~isfinite(T) || T<=0
    reason = 'time';
    return;
end

switch lower(method)

    case 'generalhcw'
        isImpulsive = true;

        Phi = hcwSTM(n,T);
        B = Phi(1:3,4:6);

        b = xf(1:3) - ...
            Phi(1:3,1:3)*x0(1:3);

        if rcond(B)<1e-10
            reason = 'singular';
            return;
        end

        vPlus = B\b;

        burn1 = vPlus-x0(4:6);

        vArrival = ...
            Phi(4:6,1:3)*x0(1:3) + ...
            Phi(4:6,4:6)*vPlus;

        burn2 = xf(4:6)-vArrival;

        dvTotal = norm(burn1)+norm(burn2);

        Xpost = [x0(1:3);vPlus];

        stateAt = @(tau) hcwPath(n,tau,Xpost);

        peakAccel = 0;  % Impulses are idealized.

    case {'vbar','rbar'}
        isImpulsive = true;

        a = analytical_manos( ...
            x0(1:3),xf(1:3), ...
            n,method,3,false);

        if ~a.feasible
            reason = 'terminal';
            return;
        end

        T = a.TOF_s;
        burn1 = a.burn1(:);
        % Derive the arrival correction from the actual propagated
        % HCW velocity, without assuming an analytical_manos burn2 field.
        Xanalytical = [x0(1:3);x0(4:6)+burn1];
        Xarrive = hcwPath(n,T,Xanalytical);
        burn2 = xf(4:6)-Xarrive(4:6);
        dvTotal = norm(burn1)+norm(burn2);

        Xpost = [x0(1:3); ...
                 x0(4:6)+a.burn1];

        stateAt = @(tau) hcwPath(n,tau,Xpost);

        peakAccel = 0;

        % Analytical maneuver assumes both endpoints
        % stationary on V-bar.
        if norm(x0(4:6))>1e-9 || ...
           norm(xf(4:6))>1e-9

            reason = 'analytical_boundary';
            return;
        end

    case 'continuous'

        % CHANGED:
        % Generalized quintic with endpoint velocities.
        % Same HCW feedforward equations as old model.

        stateAt = @(tau) ...
            quinticBoundaryState(x0,xf,T,tau);

        tControl = linspace(0,T,max(501,N));

        [~,u] = quinticBoundaryState( ...
            x0,xf,T,tControl,n);

        uMag = vecnorm(u,2,1);

        dvTotal = trapz(tControl,uMag);
        peakAccel = max(uMag);

    otherwise
        error('Unknown maneuver family: %s',method);
end

%% Propagate complete coast/reference

t = linspace(0,T,N);
trajectory = stateAt(t);

%% Terminal position validation

if norm(trajectory(1:3,end)-xf(1:3))>1e-5
    reason = 'terminal';
    return;
end

% For impulsive transfers the arrival velocity BEFORE
% the terminal burn may differ from xf(4:6).
% Its correction is included in dvTotal.

if ~isfinite(dvTotal) || ...
   ~all(isfinite(trajectory(:)))

    reason = 'numeric';
    return;
end

z = struct();

z.dvTotal = dvTotal;
z.TOF_s = T;

z.trajectory = trajectory;
z.stateAt = stateAt;
z.startState=x0; z.endState=xf; z.meanMotion=n;
z.method = string(method);

z.maxSpeed = max( ...
    vecnorm(trajectory(4:6,:),2,1));
rr = vecnorm(trajectory(1:3,:),2,1);
rhoDot = sum(trajectory(1:3,:).*trajectory(4:6,:),1)./max(rr,eps);
z.maxClosingSpeed = max([0,-rhoDot]);
z.rhoDot_mps = rhoDot;
z.isImpulsive = isImpulsive;
z.burn1 = burn1;
z.burn2 = burn2;
z.totalBurn_s = NaN;
z.maxBurn_s = NaN;
z.peakAcceleration = peakAccel;

end

%% ============================================================
% GENERALIZED QUINTIC HCW CONTINUOUS CONTROL
% =============================================================

function [Y,u] = quinticBoundaryState(x0,xf,T,t,n)
% Position/velocity quintic interpolation.
%
% Conditions:
% r(0)=r0, v(0)=v0, a(0)=0
% r(T)=rf, v(T)=vf, a(T)=0
%
% Returns 6xN LVLH state.
% If requested, also returns 3xN HCW control accel.

t = reshape(t,1,[]);
s = t/T;

r0 = x0(1:3);
v0 = x0(4:6);

rf = xf(1:3);
vf = xf(4:6);

D = rf-r0-v0*T;
V = (vf-v0)*T;

A = 10*D-4*V;
B = -15*D+7*V;
C = 6*D-3*V;

r = r0 + v0*t + ...
    A*(s.^3) + B*(s.^4) + C*(s.^5);

v = v0 + ...
    (3*A*(s.^2) + ...
     4*B*(s.^3) + ...
     5*C*(s.^4))/T;

Y = [r;v];

if nargout>1

    a = ( ...
        6*A*s + ...
        12*B*(s.^2) + ...
        20*C*(s.^3))/T^2;

    u = [ ...
        a(1,:)-3*n^2*r(1,:)-2*n*v(2,:); ...
        a(2,:)+2*n*v(1,:); ...
        a(3,:)+n^2*r(3,:)];
end

end

%% ============================================================
% HCW STATE TRANSITION MATRIX - ORIGINAL EQUATIONS
% =============================================================

function Phi = hcwSTM(n,t)

s = sin(n*t);
c = cos(n*t);
nt = n*t;

Phi = [ ...
    4-3*c,0,0,s/n,2*(1-c)/n,0; ...
    6*(s-nt),1,0,-2*(1-c)/n,(4*s-3*nt)/n,0; ...
    0,0,c,0,0,s/n; ...
    3*n*s,0,0,c,2*s,0; ...
    6*n*(c-1),0,0,-2*s,4*c-3,0; ...
    0,0,-n*s,0,0,c];

end

%% ============================================================
% VECTORIZED HCW PROPAGATION - ORIGINAL EQUATIONS
% =============================================================

function Y = hcwPath(n,t,X)

s = sin(n*t);
c = cos(n*t);
nt = n*t;

x = X(1);
y = X(2);
z = X(3);

vx = X(4);
vy = X(5);
vz = X(6);

Y = [ ...
    (4-3*c)*x + s/n*vx + 2*(1-c)/n*vy; ...
    6*(s-nt)*x + y - 2*(1-c)/n*vx + ...
        (4*s-3*nt)/n*vy; ...
    c*z + s/n*vz; ...
    3*n*s*x + c*vx + 2*s*vy; ...
    6*n*(c-1)*x - 2*s*vx + (4*c-3)*vy; ...
    -n*s*z + c*vz];

end

%% ============================================================
% RANDOM SAMPLING
% =============================================================

function v = draw(b)

v = b(1)+(b(2)-b(1))*rand;

end

function x = sampleState(nom,dev,cfg)
% Nominal rectangular perturbations or directionally uniform 2 km shell.
x = nom+(2*rand(6,1)-1).*dev;
mode = lower(string(cfg.initialPositionMode));
if mode=="alongtrack"
    return;
end

u = randn(3,1);
u = u/norm(u);
rMin = cfg.shellRadius_m(1);
rMax = cfg.shellRadius_m(2);
if mode=="shell"
    radius = rMin+(rMax-rMin)*rand;
else
    radius = (rMin^3+(rMax^3-rMin^3)*rand)^(1/3);
end

x(1:3) = radius*u;
if strcmpi(cfg.shellVelocityMode,'ball')
    w = randn(3,1); w = w/norm(w);
    x(4:6) = cfg.shellMaxRelSpeed_mps * rand^(1/3)*w;
end
end

%% ============================================================
% EMPIRICAL P95 (NEAREST RANK)
% =============================================================

function [selected,pareto] = selectCandidates(cfg,dv,tof,valid)
[N,M] = size(valid);
selected = nan(N,1);
pareto = false(N,M);
for i = 1:N
    ids = find(valid(i,:));
    if isempty(ids), continue; end
    for m = ids
        dominates = valid(i,:) & dv(i,:)<=dv(i,m)+1e-10 & ...
            tof(i,:)<=tof(i,m)+1e-10 & ...
            (dv(i,:)<dv(i,m)-1e-10 | tof(i,:)<tof(i,m)-1e-10);
        pareto(i,m) = ~any(dominates);
    end

    frontier = find(pareto(i,:));
    switch lower(string(cfg.objective))
        case 'dv'
            [~,k] = min(dv(i,ids)); selected(i) = ids(k);
        case 'tof'
            [~,k] = min(tof(i,ids)); selected(i) = ids(k);
        case {'weighted','pareto'}
            % 'pareto' leaves all non-dominated options in R.scenarioPareto;
            % chooses an illustrative normalized compromise separately.
            scaleDV = max(dv(i,frontier))-min(dv(i,frontier));
            scaleTOF = max(tof(i,frontier))-min(tof(i,frontier));
            score = cfg.weightDV*(dv(i,frontier)-min(dv(i,frontier)))/max(scaleDV,eps) + ...
                    cfg.weightTOF*(tof(i,frontier)-min(tof(i,frontier)))/max(scaleTOF,eps);
            [~,k] = min(score); selected(i) = frontier(k);
        otherwise
            error('objective must be dv, tof, weighted, or pareto');
    end
end
end

function [ok,reason] = checkBraking(z,xf,cfg)
% Necessary radial stopping-distance screen at start and AZ/KOZ entries.
% Evaluate stopping ability in the outward radial direction. It does not
% certify finite-thrust braking or guarantee terminal-position accuracy.
ok = true; reason = '';
T = z.TOF_s;
checkT = 0;
scanT = linspace(0,T,max(101,cfg.eventScanN));
X = z.stateAt(scanT);
rhoGrid = vecnorm(X(1:3,:),2,1);
for radius = [cfg.AZ_m,cfg.KOZ_m]
    idx = find(rhoGrid(1:end-1)>radius & ...
        rhoGrid(2:end)<=radius,1,'first');
    if ~isempty(idx), checkT(end+1)=scanT(idx); end %#ok<AGROW>
end

for t0 = unique(checkT)
    Y = z.stateAt(t0);
    r = Y(1:3); v = Y(4:6);
    rho = norm(r);
    inward = max(0,-dot(r,v)/max(rho,eps));
    if inward<=eps, continue; end
    outward = r/max(rho,eps);
    Fbrake = availableThrust(outward,cfg);
    if Fbrake<=0
        ok=false; reason='brake_direction'; return;
    end

    a = Fbrake/cfg.massTaxi_kg;
    availableDistance = max(0,rho-norm(xf(1:3)) ...
                                  -cfg.brakingMargin_m);
    stopDistance = cfg.brakingSafetyFactor*inward^2/(2*a);
    if stopDistance>availableDistance+1e-6
        ok=false; reason='braking_distance'; return;
    end
end
end

function [ok,reason,B] = screenTransferBurns(z,cfg)
% Impulses: screen EACH burn vector with usable thrust in its direction.
% Continuous: screen sampled control acceleration vector against usable
% directional thrust. No finite-burn repropagation in this function.
ok = true; reason = '';
B = struct('totalBurn_s',0,'fraction',0,'propellant_kg',0, ...
    'individualBurnTime_s',[0 0],'individualBurnDV_mps',[0 0]);
g0 = 9.80665;
if z.isImpulsive
    deltaV = [norm(z.burn1),norm(z.burn2)];
    burnVec = {z.burn1,z.burn2};
    mass = cfg.massTaxi_kg;
    for ib=1:2
        if deltaV(ib)<1e-12, continue; end
        F = availableThrust(burnVec{ib}/deltaV(ib),cfg);
        if F<=0
            ok=false; reason='burn_direction'; return;
        end

        mFinal = mass*exp(-deltaV(ib)/(g0*cfg.Isp_s));
        prop = mass-mFinal;
        B.individualBurnTime_s(ib) = prop*g0*cfg.Isp_s/F;
        B.individualBurnDV_mps(ib)=deltaV(ib);
        B.propellant_kg = B.propellant_kg+prop;
        mass=mFinal;
    end

    B.totalBurn_s = sum(B.individualBurnTime_s);
    B.fraction = B.totalBurn_s/z.TOF_s;
    if any(B.individualBurnTime_s>cfg.maxBurnDuration_s) ...
            || B.totalBurn_s>z.TOF_s
        ok=false; reason='burn_duration'; return;
    end

    if B.fraction>cfg.maxBurnFraction+1e-12
        ok=false; reason='impulsive_fraction'; return;
    end
else
    % Actual continuous control acceleration in LVLH, sampled densely.
    tt=linspace(0,z.TOF_s,max(101,cfg.nThrustSamples));
    [~,controlAccel]=quinticBoundaryState( ...
        z.startState,z.endState,z.TOF_s,tt,z.meanMotion);
    for j=1:numel(tt)
        req=norm(controlAccel(:,j));
        if req<1e-12, continue; end
        F=availableThrust(controlAccel(:,j)/req,cfg);
        if F<=0 || req>F/cfg.massTaxi_kg+1e-12
            ok=false; reason='continuous_thrust'; return;
        end
    end
end
end

function F = availableThrust(direction,cfg)
% If a thruster direction matrix is given, maximize commanded force along
% the requested LVLH direction while cancelling lateral force exactly.
% Thruster directions must be expressed in the SAME LVLH frame as burns.
% No attitude rotation or torques are modeled here.
d=direction(:)/norm(direction);
if isempty(cfg.thrusterDirections_LVLH)
    F=cfg.nActiveThrusters*cfg.thrust_N*cfg.thrustEfficiency;
    return;
end

U=cfg.thrusterDirections_LVLH;
U=U./vecnorm(U,2,1);
A=U.*reshape(cfg.thrusterForces_N,1,[]);
P=eye(3)-d*d';
options=optimoptions('linprog','Display','none');
[w,~,exitflag]=linprog(-(d'*A)',[],[],P*A,zeros(3,1), ...
    zeros(size(A,2),1),ones(size(A,2),1),options);
if exitflag<=0 || isempty(w)
    F=0;
else
    F=max(0,d'*A*w)*cfg.thrustEfficiency;
end
end

function p = p95(a)

a = sort(a(isfinite(a)));

if isempty(a)
    p = NaN;
    return;
end

p = a(max(1,ceil(.95*numel(a))));

end

function [failure,events] = verifyZoneEvents( ...
    z,phaseIndex,cfg,final_m)
% VERIFYZONEEVENTS
%
% Transfer index:
% 1 = Homing 1
% 2 = Homing 2
% 3 = Closing 1 (authorized AZ crossing)
% 4 = Closing 2 (authorized KOZ crossing)
% 5 = Final approach
%
% No lateral corridor assumption is used.
%
% Event times are refined using fzero.

failure = '';

Nscan = max(101,cfg.eventScanN);

tt = linspace(0,z.TOF_s,Nscan);

S = z.stateAt(tt);

rho = vecnorm(S(1:3,:),2,1);

AZ = detectSphereCrossings( ...
    z,cfg.AZ_m,tt,rho,cfg);

KOZ = detectSphereCrossings( ...
    z,cfg.KOZ_m,tt,rho,cfg);

events = struct();
events.AZ = AZ;
events.KOZ = KOZ;
events.minRange_m = min(rho);
events.maxRange_m = max(rho);

tol = cfg.zoneTolerance_m;

switch phaseIndex

    case {1,2}
        % Before AZ authorization:
        % Entire phase must remain outside AZ.

        if any(rho<cfg.AZ_m-tol)
            failure = 'premature_AZ';
            return;
        end

    case 3
        % Closing 1: AZ entry allowed.
        % KOZ entry NOT yet authorized.

        if any(rho<cfg.KOZ_m-tol)
            failure = 'premature_KOZ';
            return;
        end

        % Exactly one inward AZ crossing.
        if numel(AZ.time_s)~=1 || ...
           AZ.direction(1)~="inward"

            failure = 'AZ_crossing';
            return;
        end

        % Must finish within AZ.
        if rho(end)>=cfg.AZ_m
            failure = 'AZ_endpoint';
            return;
        end

    case 4
        % Closing 2: KOZ entry authorized.

        % Do not re-exit AZ.
        if any(rho>cfg.AZ_m+tol)
            failure = 'AZ_reexit';
            return;
        end

        % Exactly one inward KOZ crossing.
        if numel(KOZ.time_s)~=1 || ...
           KOZ.direction(1)~="inward"

            failure = 'KOZ_crossing';
            return;
        end

        % No passage closer than final standoff.
        if any(rho<final_m-tol)
            failure = 'target_passage';
            return;
        end

    case 5
        % Final approach is already inside KOZ.
        % Remain inside KOZ, outside target standoff.

        if any(rho>cfg.KOZ_m+tol)
            failure = 'KOZ_reexit';
            return;
        end

        if any(rho<final_m-tol)
            failure = 'target_passage';
            return;
        end
end

end

%% ============================================================
% SPHERICAL BOUNDARY CROSSING EVENTS
% =============================================================

function E = detectSphereCrossings( ...
    z,radius,tt,rho,cfg)

g = rho-radius;

% Ignore exact-zero grid nodes while searching for
% intervals with opposite signs.
idx = find(abs(g)>1e-10);

tEvents = [];
rateEvents = [];
directions = strings(0,1);

for k = 1:numel(idx)-1

    i1 = idx(k);
    i2 = idx(k+1);

    if sign(g(i1))==sign(g(i2))
        continue;
    end

    % Refine event on the ACTUAL propagated trajectory.
    tr = fzero( ...
        @(tau) sphereResidual(z,tau,radius), ...
        [tt(i1),tt(i2)]);

    % Avoid duplicate root detections.
    if ~isempty(tEvents) && ...
       abs(tr-tEvents(end))<1e-5
        continue;
    end

    Y = z.stateAt(tr);

    r = Y(1:3);
    v = Y(4:6);

    rr = norm(r);

    rangeRate = dot(r,v)/rr;

    if rangeRate < -cfg.eventRangeRateTol_mps
        direction = "inward";

    elseif rangeRate > cfg.eventRangeRateTol_mps
        direction = "outward";

    else
        direction = "tangent";
    end

    tEvents(end+1,1) = tr; %#ok<AGROW>
    rateEvents(end+1,1) = rangeRate; %#ok<AGROW>

    directions(end+1,1) = direction; %#ok<AGROW>
end

E = struct( ...
    'radius_m',radius, ...
    'time_s',tEvents, ...
    'time_min',tEvents/60, ...
    'rangeRate_mps',rateEvents, ...
    'direction',directions);

end

function g = sphereResidual(z,t,radius)

Y = z.stateAt(t);
g = norm(Y(1:3))-radius;

end

function fuel = cyclePropellant(clientDV,hostDV,cfg)
% CYCLEPROPELLANT
% Compute ideal cycle fuel consumption using the Isp assigned
% to each of the nine phases.
%
% Mass decreases continuously across both approaches.
% Payload capture/release mass changes are not modeled.

g0 = 9.80665;

initialMass = cfg.massTaxi_kg;
mass = initialMass;

for approach = 1:2

    if approach == 1
        phaseDV = clientDV;
    else
        phaseDV = hostDV;
    end

    for j = 1:9

        phaseCfg = FABER_RPO_phase_propulsion(cfg,j);

        mass = mass*exp( ...
            -phaseDV(j)/(g0*phaseCfg.Isp_s));

    end
end

fuel = initialMass-mass;

end
