function result = assembly_relocations_v4(mu,Re,h,Npts,timesMin,hostPoint,modulePoints,opts)
% ASSEMBLY_RELOCATIONS_V4
%
% FABER assembly module relocation analysis.
%
% Each direction:
%   1. Controlled departure: source -> exit gate
%   2. Two-impulse HCW: exit gate -> entrance gate
%   3. Controlled arrival: entrance gate -> target
%
% The return journey follows the same sequence in reverse.
%
% METHOD OPTIONS:
%   GeneralHCW - gates on source-target line
%   Vbar       - gates along local along-track direction
%   Rbar       - gates along local radial direction
%
% INPUTS:
%   mu           Earth gravitational parameter [km^3/s^2]
%   Re           Earth radius [km]
%   h            Orbit altitude [km]
%   Npts         Number of HCW trajectory samples
%   timesMin     Candidate HCW transfer durations [min]
%   hostPoint    Host assembly interface [3x1] [m]
%   modulePoints Module interfaces [Nmodules x 3] [m]
%   opts         Options struct
%
% Requires HCW_solver.m.
%
% Limitations:
%   - Ideal impulsive HCW transfer
%   - Ideal continuous quintic departure/arrival
%   - No finite-thrust simulation for impulsive burns
%   - No complete spacecraft/truss collision geometry
%   - No abort or missed-burn analysis
%   - Safety checked at discrete trajectory samples
%
% SELECTION OBJECTIVE (opts.objective):
%   dv       - minimum DV (legacy result)
%   tof      - minimum TOF
%   weighted - equal normalized DV and TOF importance by default
% J = wDV*(DV/opts.dvReference_mps) + wTOF*(TOF/opts.tofReference_min)
% The same additive score is used for each segment and method.

%% ================================================================
% 1. INITIALIZATION
% ================================================================

if nargin < 8 || isempty(opts)
    opts = struct();
end

opts = defaultOpts(opts);

assert(ismember(lower(string(opts.objective)), ...
    ["dv","tof","weighted"]), ...
    'opts.objective must be dv, tof, or weighted.');
assert(opts.dvReference_mps>0 && isfinite(opts.dvReference_mps));
assert(opts.tofReference_min>0 && isfinite(opts.tofReference_min));
assert(all(isfinite([opts.weightDV opts.weightTOF])) && ...
    opts.weightDV>=0 && opts.weightTOF>=0 && ...
    opts.weightDV+opts.weightTOF>0, ...
    'Weights must be nonnegative and at least one positive.');

assert(exist('HCW_solver','file')==2, ...
    'HCW_solver.m not found.');

assert(isfinite(mu) && mu>0);
assert(isfinite(Re) && Re>0);
assert(isfinite(h));

assert(isscalar(Npts) && Npts>=3 && Npts==floor(Npts));

assert(isvector(timesMin) && ~isempty(timesMin) && ...
    all(isfinite(timesMin)) && all(timesMin>0));

assert(isvector(opts.exitTimesMin) && ...
    ~isempty(opts.exitTimesMin) && ...
    all(isfinite(opts.exitTimesMin)) && ...
    all(opts.exitTimesMin>0));

assert(isvector(opts.finalTimesMin) && ...
    ~isempty(opts.finalTimesMin) && ...
    all(isfinite(opts.finalTimesMin)) && ...
    all(opts.finalTimesMin>0));

assert(size(modulePoints,2)==3, ...
    'modulePoints must be Nmodules-by-3.');

hostPoint = hostPoint(:);

assert(numel(hostPoint)==3 && all(isfinite(hostPoint)));
assert(all(isfinite(modulePoints(:))));

assert(opts.gateDistance_m>0);
assert(opts.minTargetRange_m>=0);
validateattributes(opts.maxTransferDV_mps,{'numeric'},{'scalar','real','nonnegative'});

assert(all([opts.maxTransferSpeed_mps, ...
    opts.maxExitSpeed_mps, ...
    opts.maxFinalSpeed_mps]>0));

assert(all([opts.maxExitAccel_mps2, ...
    opts.maxFinalAccel_mps2]>0));

% Mean orbital motion [rad/s]
n = sqrt(mu/(Re+h)^3);

methods = string(opts.methods(:));

assert(~isempty(methods) && ...
    all(ismember(lower(methods), ...
    ["generalhcw";"vbar";"rbar"])));

Nmodules = size(modulePoints,1);
M = numel(methods);

%% ================================================================
% 2. PREALLOCATE RESULTS
% ================================================================

DV = NaN(Nmodules,M);
TOF = NaN(Nmodules,M);
Scores = NaN(Nmodules,M);

DVout = NaN(Nmodules,M);
DVback = NaN(Nmodules,M);

TOFout = NaN(Nmodules,M);
TOFback = NaN(Nmodules,M);

feasible = false(Nmodules,M);
records = cell(Nmodules,M);
rejectionByDirection = strings(Nmodules,M,2);

%% ================================================================
% 3. COMPUTE ROUND TRIPS
% ================================================================

for i = 1:Nmodules

    modulePos = modulePoints(i,:)';

    for j = 1:M

        method = char(methods(j));

        % Host -> Module
        outbound = optimizeDirection( ...
            hostPoint,modulePos,method,n, ...
            Npts,timesMin,opts);

        % Module -> Host
        inbound = optimizeDirection( ...
            modulePos,hostPoint,method,n, ...
            Npts,timesMin,opts);

        rejectionByDirection(i,j,:) = reshape([outbound.reason,inbound.reason],1,1,2);

        % Both directions must be feasible
        if ~(outbound.feasible && inbound.feasible)
            continue;
        end

        feasible(i,j) = true;

        DVout(i,j) = outbound.dv_mps;
        DVback(i,j) = inbound.dv_mps;

        DV(i,j) = outbound.dv_mps + inbound.dv_mps;

        TOFout(i,j) = outbound.tof_min;
        TOFback(i,j) = inbound.tof_min;

        TOF(i,j) = outbound.tof_min + inbound.tof_min;
        Scores(i,j) = objectiveCost(DV(i,j),TOF(i,j),opts);

        records{i,j} = struct( ...
            'outbound',outbound, ...
            'returnTrip',inbound);

    end
end

%% ================================================================
% 4. SUMMARY BY METHOD
% ================================================================

% A method is globally feasible only if it works for every module.
methodFeasible = all(feasible,1)';

methodDV = NaN(M,1);
methodTOF = NaN(M,1);

for j = 1:M

    if methodFeasible(j)

        methodDV(j) = sum(DV(:,j));
        methodTOF(j) = sum(TOF(:,j));

    end
end

summary = table( ...
    methods,methodFeasible,methodDV,methodTOF, ...
    'VariableNames', ...
    {'Method','Feasible','TotalDV_mps','TotalTOF_min'});

%% ================================================================
% 5. SELECT BEST METHOD PER MODULE
% ================================================================

bestMethod = strings(Nmodules,1);

DV_perModule = NaN(Nmodules,1);
TOF_perModule = NaN(Nmodules,1);

DV_outbound = NaN(Nmodules,1);
DV_return = NaN(Nmodules,1);

% Columns: [outbound, return]
DV_exit = NaN(Nmodules,2);
DV_transfer = NaN(Nmodules,2);
DV_final = NaN(Nmodules,2);

maxExitSpeed = NaN(Nmodules,2);
maxTransferSpeed = NaN(Nmodules,2);
maxFinalSpeed = NaN(Nmodules,2);

pathsOut = cell(Nmodules,1);
pathsReturn = cell(Nmodules,1);

timePerLeg = NaN(Nmodules,2);
bestRecords = cell(Nmodules,1);

totalTimePerDirection_min = nan(Nmodules,2);
for i = 1:Nmodules

    [bestScore,j] = min(Scores(i,:),[],'omitnan');

    if isnan(bestScore) || isinf(bestScore)
        continue;
    end

    r = records{i,j};

    a = r.outbound;
    b = r.returnTrip;

    DV_perModule(i) = DV(i,j);
    TOF_perModule(i) = TOF(i,j);

    DV_outbound(i) = DVout(i,j);
    DV_return(i) = DVback(i,j);

    DV_exit(i,:) = [ ...
        a.exitDV_mps,b.exitDV_mps];

    DV_transfer(i,:) = [ ...
        a.transferDV_mps,b.transferDV_mps];

    DV_final(i,:) = [ ...
        a.finalDV_mps,b.finalDV_mps];

    maxExitSpeed(i,:) = [ ...
        a.exitMaxSpeed_mps,b.exitMaxSpeed_mps];

    maxTransferSpeed(i,:) = [ ...
        a.transferMaxSpeed_mps,b.transferMaxSpeed_mps];

    maxFinalSpeed(i,:) = [ ...
        a.finalMaxSpeed_mps,b.finalMaxSpeed_mps];

    bestMethod(i) = methods(j);

    % Combine all three trajectory segments.
    % Do not duplicate the shared gate endpoints.
    pathsOut{i} = [ ...
        a.exitPath(1:3,:), ...
        a.transferPath(1:3,2:end), ...
        a.finalPath(1:3,2:end)];

    pathsReturn{i} = [ ...
        b.exitPath(1:3,:), ...
        b.transferPath(1:3,2:end), ...
        b.finalPath(1:3,2:end)];

    % HCW transfer times [outbound,return]
    timePerLeg(i,:) = [ ...
        a.transferTOF_min,b.transferTOF_min];

    totalTimePerDirection_min(i,:) = [a.tof_min,b.tof_min];
    bestRecords{i} = r;

end

%% ================================================================
% 6. OUTPUT STRUCTURE
% ================================================================

result = struct();

result.feasible = all(isfinite(DV_perModule));
result.Nmodules = Nmodules;

result.methods = methods;
result.summary = summary;

result.feasibleByModuleMethod = feasible;

result.DV_byModuleMethod_mps = DV;
result.ObjectiveScore_byModuleMethod = Scores;
result.TOF_byModuleMethod_min = TOF;

result.records = records;
result.rejectionByDirection = rejectionByDirection;

result.bestMethodPerModule = bestMethod;
result.bestRecords = bestRecords;

result.DV_perModule = DV_perModule;

result.DV_outbound = DV_outbound;
result.DV_return = DV_return;

result.DV_exit_mps = DV_exit;
result.DV_transfer_mps = DV_transfer;
result.DV_final_mps = DV_final;

result.maxExitSpeed_mps = maxExitSpeed;
result.maxTransferSpeed_mps = maxTransferSpeed;
result.maxFinalSpeed_mps = maxFinalSpeed;

result.TOF_perModule_min = TOF_perModule;
result.timePerLeg_min = timePerLeg;
result.totalTimePerDirection_min = totalTimePerDirection_min;

result.pathsOut = pathsOut;
result.pathsReturn = pathsReturn;

result.orbitalMeanMotion_radps = n;
result.options = opts;
result.objective = string(opts.objective);

if result.feasible

    result.DV = sum(DV_perModule);
    result.TOF_min = sum(TOF_perModule);

else

    result.DV = NaN;
    result.TOF_min = NaN;

end

%% ================================================================
% 7. PRINT SUMMARY
% ================================================================

if opts.printSummary

    fprintf('\n========================================\n');
    fprintf('FABER ASSEMBLY RELOCATIONS V4\n');
    fprintf('========================================\n');

    fprintf('Altitude: %.2f km\n',h);
    fprintf('Modules: %d\n',Nmodules);
    fprintf('Gate radius: %.1f m\n', ...
        opts.gateDistance_m);

    disp(summary);

    disp(table( ...
        (1:Nmodules)', ...
        bestMethod, ...
        DV_outbound, ...
        DV_return, ...
        DV_perModule, ...
        TOF_perModule, ...
        'VariableNames', ...
        {'Module','Method','DV_Out_mps', ...
         'DV_Return_mps','DV_Total_mps', ...
         'TOF_Total_min'}));

    if result.feasible

        fprintf('Best-per-module total DV: %.5f m/s\n', ...
            result.DV);

        fprintf('Best-per-module total TOF: %.2f min\n', ...
            result.TOF_min);

    else

        fprintf('Some modules have no feasible round-trip candidate.\n');

    end
end

% Detailed trajectory diagnostics are available in result.bestRecords.
% No unconditional display during sensitivity sweeps.

end

%% ================================================================
% LOCAL FUNCTION: OPTIMIZE A SINGLE DIRECTION
% ================================================================

function out = optimizeDirection( ...
    source,target,method,n,Npts,timesMin,opts)

out = struct('feasible',false,'reason',"gate separation");

g = opts.gateDistance_m;

% Require enough separation for two distinct gates
if norm(target-source)<=2*g
    return;
end

%% 1. DEFINE EXIT AND ENTRANCE GATES

[exitGate,entryGate] = chooseGates( ...
    source,target,method,g);

% Gates must lie gate distance from their own interfaces.
if norm(exitGate-source)<g-1e-7 || ...
   norm(entryGate-target)<g-1e-7
    return;
end

% Neither gate may lie within the opposite interface sphere.
if norm(exitGate-target)<=g || ...
   norm(entryGate-source)<=g
    return;
end

if norm(exitGate-entryGate)<=1e-8
    return;
end

%% ================================================================
% SEGMENT 1: CONTROLLED EXIT
% ================================================================

exitSeg = bestControlledSegment( ...
    source,exitGate,n,max(Npts,201), ...
    opts.exitTimesMin, ...
    opts.maxExitSpeed_mps, ...
    opts.maxExitAccel_mps2,target,g,opts);

if ~exitSeg.feasible
    out.reason = "no feasible departure duration";
    return;
end

%% ================================================================
% SEGMENT 2: TWO-IMPULSE HCW TRANSFER
% ================================================================

transfer = struct('feasible',false,'dv',Inf,'cost',Inf);

for tMin = timesMin(:)'

    % Gate-to-gate transfer, stationary before/after the impulses.
    sol = HCW_solver(exitGate,zeros(3,1),entryGate,zeros(3,1), ...
        60*tMin,n,Npts);
    X = sol.trajectory;
    if ~sol.ok || sol.positionError > 1e-3 || sol.velocityError > 1e-6
        continue;
    end

    speed = max(vecnorm(X(4:6,:),2,1));
    rhoSource = vecnorm(X(1:3,:)-source,2,1);
    rhoTarget = vecnorm(X(1:3,:)-target,2,1);
    if speed > opts.maxTransferSpeed_mps || min(rhoSource) < g-1e-5 || ...
            min(rhoTarget) < max(g-1e-5,opts.minTargetRange_m)
        continue;
    end

    if sol.dvTotal > opts.maxTransferDV_mps || ...
            FABER_path_clearance(X(1:3,:),opts.pathGeometry) <= 0
        continue;
    end

    % Select minimum objective over feasible transfer times.
    cost = objectiveCost(sol.dvTotal,tMin,opts);
    if cost<transfer.cost

        transfer = struct( ...
            'feasible',true, ...
            'dv',sol.dvTotal, ...
            'cost',cost, ...
            'T_min',tMin, ...
            'path',X, ...
            'maxSpeed_mps',speed, ...
            'minSourceRange_m',min(rhoSource), ...
            'minTargetRange_m',min(rhoTarget));

    end

end

if ~transfer.feasible
    out.reason = "no feasible coast duration";
    return;
end

%% ================================================================
% SEGMENT 3: CONTROLLED ARRIVAL
% ================================================================

finalSeg = bestControlledSegment( ...
    entryGate,target,n,max(Npts,201), ...
    opts.finalTimesMin, ...
    opts.maxFinalSpeed_mps, ...
    opts.maxFinalAccel_mps2,source,g,opts);

if ~finalSeg.feasible
    out.reason = "no feasible arrival duration";
    return;
end

%% ================================================================
% COMBINE THREE SEGMENTS
% ================================================================

out.feasible = true;
out.reason = "";

out.method = string(method);
out.source_m = source;
out.target_m = target;

out.exitGate_m = exitGate;
out.entranceGate_m = entryGate;

% Compatibility with earlier arrival-gate output
out.gate_m = entryGate;

%% DELTA V

out.exitDV_mps = exitSeg.dv;
out.transferDV_mps = transfer.dv;
out.finalDV_mps = finalSeg.dv;

out.dv_mps = ...
    exitSeg.dv + transfer.dv + finalSeg.dv;

%% TIME OF FLIGHT

out.exitTOF_min = exitSeg.T_min;
out.transferTOF_min = transfer.T_min;
out.finalTOF_min = finalSeg.T_min;

out.tof_min = ...
    exitSeg.T_min + transfer.T_min + finalSeg.T_min;

%% TRAJECTORIES

out.exitPath = exitSeg.path;
out.transferPath = transfer.path;
out.finalPath = finalSeg.path;

%% MAXIMUM SPEEDS

out.exitMaxSpeed_mps = exitSeg.maxSpeed_mps;
out.transferMaxSpeed_mps = transfer.maxSpeed_mps;
out.finalMaxSpeed_mps = finalSeg.maxSpeed_mps;

%% MAXIMUM CONTINUOUS ACCELERATIONS

out.exitMaxAccel_mps2 = exitSeg.maxAccel_mps2;
out.finalMaxAccel_mps2 = finalSeg.maxAccel_mps2;

%% MINIMUM HCW CLEARANCE

out.transferMinSourceRange_m = ...
    transfer.minSourceRange_m;

out.transferMinTargetRange_m = ...
    transfer.minTargetRange_m;

end

%% ================================================================
% LOCAL FUNCTION: FIND BEST CONTROLLED SEGMENT
% ================================================================

function seg = bestControlledSegment( ...
    r0,rf,n,N,times,limSpeed,limAccel,opposite,gateRadius,opts)

seg = struct('feasible',false,'dv',Inf,'cost',Inf);

for tMin = times(:)'

    T = 60*tMin;

    t = linspace(0,T,N);

    % Ideal continuously controlled HCW trajectory
    [X,u] = quinticHCW(r0,rf,T,n,t);

    if any(~isfinite(X(:))) || ...
       any(~isfinite(u(:)))
        continue;
    end

    % Actual relative speed along reference
    speeds = vecnorm(X(4:6,:),2,1);

    % Required continuous feedforward acceleration
    acc = vecnorm(u,2,1);

    vMax = max(speeds);
    aMax = max(acc);

    if vMax>limSpeed || aMax>limAccel
        continue;
    end

    if min(vecnorm(X(1:3,:)-opposite,2,1)) < gateRadius-1e-5 || ...
            FABER_path_clearance(X(1:3,:),opts.pathGeometry) <= 0
        continue;
    end

    % Ideal continuous-control Delta V
    dv = trapz(t,acc);

    cost = objectiveCost(dv,tMin,opts);

    if cost<seg.cost

        seg = struct( ...
            'feasible',true, ...
            'dv',dv, ...
            'cost',cost, ...
            'T_min',tMin, ...
            'path',X, ...
            'maxSpeed_mps',vMax, ...
            'maxAccel_mps2',aMax);

    end

end

end

%% ================================================================
% LOCAL FUNCTION: GATE GEOMETRY
% ================================================================

function [exitGate,entryGate] = chooseGates( ...
    source,target,method,g)

delta = target-source;

switch lower(method)

    case 'generalhcw'

        % Along line connecting source and target
        u = delta/norm(delta);

        exitGate = source+g*u;
        entryGate = target-g*u;

    case 'vbar'

        % Along-track gate geometry
        u = [0;signOrPlus(delta(2));0];

        exitGate = source+g*u;
        entryGate = target-g*u;

    case 'rbar'

        % Radial gate geometry
        u = [signOrPlus(delta(1));0;0];

        exitGate = source+g*u;
        entryGate = target-g*u;

    otherwise

        error('Unknown gate geometry: %s',method);

end

end

function s = signOrPlus(a)

s = sign(a);

if s==0
    s = 1;
end

end

%% ================================================================
% LOCAL FUNCTION: QUINTIC HCW TRAJECTORY
% ================================================================

function [X,u] = quinticHCW(r0,rf,T,n,t)
% Stationary-to-stationary controlled relative translation.
%
% Enforces zero relative velocity and acceleration at endpoints
% for the prescribed position profile.
%
% Includes HCW feedforward acceleration.

s = reshape(t,1,[])/T;

D = rf(:)-r0(:);

% Quintic position profile
q = 10*s.^3-15*s.^4+6*s.^5;

% First derivative
qd = (30*s.^2-60*s.^3+30*s.^4)/T;

% Second derivative
qdd = (60*s-180*s.^2+120*s.^3)/T^2;

% Desired relative motion
r = r0(:)+D*q;
v = D*qd;
a = D*qdd;

X = [r;v];

% Required control acceleration from HCW
u = [ ...
    a(1,:)-3*n^2*r(1,:)-2*n*v(2,:); ...
    a(2,:)+2*n*v(1,:); ...
    a(3,:)+n^2*r(3,:)];

end

%% ================================================================
% LOCAL FUNCTION: OBJECTIVE FUNCTION
% ================================================================

function J = objectiveCost(dv_mps,tof_min,opts)
% The reference values are FIXED across all candidates and modules.
% This ensures consistent comparison and additive phase costs.

switch lower(char(opts.objective))
    case 'dv'
        J = dv_mps;
    case 'tof'
        J = tof_min;
    case 'weighted'
        wsum = opts.weightDV+opts.weightTOF;
        J = (opts.weightDV/wsum) * ...
            (dv_mps/opts.dvReference_mps) + ...
            (opts.weightTOF/wsum) * ...
            (tof_min/opts.tofReference_min);
    otherwise
        error('Unknown objective: %s',string(opts.objective));
end
end

%% ================================================================
% LOCAL FUNCTION: DEFAULT OPTIONS
% ================================================================

function opts = defaultOpts(opts)

defs = struct( ...
    'methods',{{'GeneralHCW','Vbar','Rbar'}}, ...
    'gateDistance_m',30, ...
    'exitTimesMin',5:5:30, ...
    'finalTimesMin',5:5:30, ...
    'maxTransferSpeed_mps',Inf, ...
    'maxTransferDV_mps',Inf, ...
    'pathGeometry',struct(), ...
    'maxExitSpeed_mps',Inf, ...
    'maxFinalSpeed_mps',Inf, ...
    'maxExitAccel_mps2',Inf, ...
    'maxFinalAccel_mps2',Inf, ...
    'minTargetRange_m',0, ...
    'printSummary',true, ...
    'objective','weighted', ...
    'weightDV',0.5, ...
    'weightTOF',0.5, ...
    'dvReference_mps',0.1, ...
    'tofReference_min',60);

f = fieldnames(defs);

for k = 1:numel(f)

    if ~isfield(opts,f{k})
        opts.(f{k}) = defs.(f{k});
    end

end

end
