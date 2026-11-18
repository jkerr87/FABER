function R = FABER_RPO_9phase_MC(cfg)
% FABER_RPO_9PHASE_MC
% Nine-phase, four-holding-point RPO Monte Carlo.
%
% Compares 18 combinations:
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

% Nine phases for each approach
DVphase = nan(N,M,2,9);
TOFphase = nan(N,M,2,9);

reason = strings(N,M);
records = cell(N,M);

%% 3. PAIRED MONTE CARLO

for i = 1:N

    % Same initial states used across all method combinations
    xc = sampleState(cfg.x0_client,cfg.dx0_client);
    xh = sampleState(cfg.x0_host,cfg.dx0_host);

    % Same sampled design used across all methods
    dc = sampleDesign(cfg,xc);
    dh = sampleDesign(cfg,xh);

    for m = 1:M

        which = { ...
            char(methodH2(m)), ...
            char(methodC1(m)), ...
            char(methodC2(m))};


        % Taxi -> Client
        [client,why] = simulateApproach( ...
            cfg,xc,dc,n,which,char(gateMode(m)));
        
        if ~isempty(why)
            reason(i,m) = "Client: " + string(why);
            continue;
        end
        
        % Taxi -> ATLANTE
        [host,why] = simulateApproach( ...
            cfg,xh,dh,n,which,char(gateMode(m)));
        
        if ~isempty(why)
            reason(i,m) = "ATLANTE: " + string(why);
            continue;
        end


        if client.dv+host.dv > ...
                cfg.maxDV_perCycle_mps
            reason(i,m) = "Cycle DV limit";
            continue;
        end

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
    P95CycleTOF_min,CommonP95CycleDV_mps);

R.cfg = cfg;
R.methods = [ ...
    methodH2 methodC1 methodC2 gateMode];


R.DV_mps = DV;
R.TOF_min = TOF;

R.DV_phase_mps = DVphase;
R.TOF_phase_min = TOFphase;

R.CycleDV_mps = cycleDV;
R.CycleTOF_min = cycleTOF;

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
    reason(:,m), ...
    'VariableNames',{ ...
    'Attempt','Architecture', ...
    'Homing2','Closing1','Closing2', ...
    'HP4Gate', ...
    'Feasible','CycleDV_mps', ...
    'CycleTOF_min','Reason'});


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
dvPhase = zeros(1,9);
tofPhase = zeros(1,9);
zoneEvents = cell(1,5);
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

    % HARD TOTAL RELATIVE-SPEED CONSTRAINT

    if z.maxSpeed > maxSpeed(k)+1e-9

        reason = sprintf('phase%d_speed', ...
            phaseIndices(k));
        return;
    end

    % OPTIONAL ACCELERATION LIMIT

    if z.peakAcceleration > ...
            cfg.maxPeakAccel_mps2

        reason = sprintf('phase%d_acceleration', ...
            phaseIndices(k));
        return;
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
    dvPhase(j) = z.dvTotal;
    tofPhase(j) = z.TOF_s/60;
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
end

%% 6. RETURN RESULTS

s = struct( ...
    'dv',sum(dvPhase), ...
    'tof',sum(tofPhase), ...
    'dvPhase',dvPhase, ...
    'tofPhase',tofPhase, ...
    'trace',{traces}, ...
    'GO',true, ...
    'methods',{methods}, ...
    'holdingPoints_m',hp);

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

if ~isfinite(T) || T<=0
    reason = 'time';
    return;
end

switch lower(method)

    case 'generalhcw'

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

        a = analytical_manos( ...
            x0(1:3),xf(1:3), ...
            n,method,3,false);

        if ~a.feasible
            reason = 'terminal';
            return;
        end

        T = a.TOF_s;
        dvTotal = a.DV_mps;

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

z.maxSpeed = max( ...
    vecnorm(trajectory(4:6,:),2,1));

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

function x = sampleState(nom,dev)

x = nom+(2*rand(6,1)-1).*dev;

end

%% ============================================================
% EMPIRICAL P95 (NEAREST RANK)
% =============================================================

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
