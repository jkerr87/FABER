function result = assembly_relocations_v3(mu,Re,h,Npts,timesMin,hostPoint,modulePoints,opts)
% ASSEMBLY_RELOCATIONS_V3
% Round trip for each module:
%   Host -> 30-m module gate (HCW) -> module interface (continuous)
%   Module -> 30-m host gate   (HCW) -> host interface   (continuous)
%
% METHODS (these are GATE GEOMETRIES, not analytic V-bar/R-bar solutions):
%   'GeneralHCW' : gate on line from departure to target
%   'Vbar'       : gate 30 m from target on +/- along-track (y) axis
%   'Rbar'       : gate 30 m from target on +/- radial (x) axis
% All transfer-to-gate legs use existing two-impulse HCW_solver.
% Both final approaches use a controlled quintic reference with HCW
% feedforward acceleration and integral of |u| dt as ideal thrust DV.
%
% INPUTS
%   mu,Re,h       [km^3/s^2], [km], [km]
%   Npts          trajectory samples per leg
%   timesMin      candidate HCW transfer times [min], vector
%   hostPoint     [3x1] LVLH location [m]
%   modulePoints  [Nmodules x 3] LVLH locations [m]
%   opts          optional struct, fields below
%
% OPTIONS
%   methods               default {'GeneralHCW','Vbar','Rbar'}
%   gateDistance_m        default 30
%   finalTimesMin         default 5:5:30 [min]
%   maxTransferSpeed_mps  default Inf
%   maxFinalSpeed_mps     default Inf
%   maxFinalAccel_mps2    default Inf
%   minTargetRange_m      default 0 [m] (0 disables exclusion)
%   printSummary          default true
%
% CAUTIONS
%  - The input target points MUST be assembly/capture interfaces or
%    appropriate offset locations, not necessarily body centres.
%  - No 3D body/arm collision clearance, pointing, plume, or abort model.
%  - Transfer burns are ideal impulses; final is ideal continuous control.
%  - No propellant mass, finite thrust, or attitude dynamics.
%  - This is a minimum-DV discrete timing study, not a TOF Pareto solver.
%
% REQUIRES existing HCW_solver.m with sol.dvTotal and sol.trajectory.

% Maximum relative speeds [m/s]
opts.maxTransferSpeed_mps = 0.3;
opts.maxFinalSpeed_mps = 0.1;

if nargin < 8 || isempty(opts), opts = struct(); end
opts = defaultOpts(opts);
assert(exist('HCW_solver','file')==2,'HCW_solver.m not found');
assert(isfinite(mu) && mu>0 && isfinite(Re) && Re>0 && isfinite(h));
assert(isscalar(Npts) && Npts>=3 && Npts==floor(Npts));
assert(isvector(timesMin) && all(isfinite(timesMin)) && all(timesMin>0));
assert(isvector(opts.finalTimesMin) && all(isfinite(opts.finalTimesMin)) && all(opts.finalTimesMin>0));
assert(size(modulePoints,2)==3,'modulePoints must be Nmodules-by-3');
hostPoint = hostPoint(:);
assert(numel(hostPoint)==3 && all(isfinite(hostPoint)));
assert(all(isfinite(modulePoints(:))));
assert(opts.gateDistance_m>0 && opts.minTargetRange_m>=0);

n = sqrt(mu/(Re+h)^3);                    % rad/s
methods = string(opts.methods(:));
assert(all(ismember(lower(methods),["generalhcw";"vbar";"rbar"])));
Nmodules = size(modulePoints,1);
M = numel(methods);

DV = NaN(Nmodules,M);
TOF = NaN(Nmodules,M);
DVout = NaN(Nmodules,M);
DVback = NaN(Nmodules,M);
TOFout = NaN(Nmodules,M);
TOFback = NaN(Nmodules,M);
feasible = false(Nmodules,M);
records = cell(Nmodules,M);

for i = 1:Nmodules
    modulePos = modulePoints(i,:)';
    for j = 1:M
        method = char(methods(j));

        out = optimizeDirection(hostPoint,modulePos,method,n,Npts,timesMin,opts);
        back = optimizeDirection(modulePos,hostPoint,method,n,Npts,timesMin,opts);
        if ~(out.feasible && back.feasible), continue; end

        feasible(i,j) = true;
        DVout(i,j) = out.dv_mps;
        DVback(i,j) = back.dv_mps;
        DV(i,j) = out.dv_mps + back.dv_mps;
        TOFout(i,j) = out.tof_min;
        TOFback(i,j) = back.tof_min;
        TOF(i,j) = out.tof_min + back.tof_min;
        records{i,j} = struct('outbound',out,'returnTrip',back);
    end
end

% A method used for ALL modules must have every module feasible.
methodFeasible = all(feasible,1)';
methodDV = NaN(M,1); methodTOF = NaN(M,1);
for j=1:M
    if methodFeasible(j)
        methodDV(j) = sum(DV(:,j));
        methodTOF(j) = sum(TOF(:,j));
    end
end
summary = table(methods,methodFeasible,methodDV,methodTOF, ...
    'VariableNames',{'Method','Feasible','TotalDV_mps','TotalTOF_min'});

% Choose the best method PER MODULE (could vary between modules).
bestMethod = strings(Nmodules,1);
DV_perModule = NaN(Nmodules,1);
TOF_perModule = NaN(Nmodules,1);
DV_outbound = NaN(Nmodules,1);
DV_return = NaN(Nmodules,1);
pathsOut = cell(Nmodules,1);
pathsReturn = cell(Nmodules,1);
timePerLeg = NaN(Nmodules,2); % HCW portion only: outbound, return
bestRecords = cell(Nmodules,1);
for i = 1:Nmodules
    [v,j] = min(DV(i,:),[],'omitnan');
    if isnan(v) || isinf(v), continue; end
    DV_perModule(i) = v;
    TOF_perModule(i) = TOF(i,j);
    DV_outbound(i) = DVout(i,j);
    DV_return(i) = DVback(i,j);
    bestMethod(i) = methods(j);
    r = records{i,j};
    pathsOut{i} = [r.outbound.transferPath(1:3,:),r.outbound.finalPath(1:3,2:end)];
    pathsReturn{i} = [r.returnTrip.transferPath(1:3,:),r.returnTrip.finalPath(1:3,2:end)];
    timePerLeg(i,:) = [r.outbound.transferTOF_min,r.returnTrip.transferTOF_min];
    bestRecords{i} = r;
end

result = struct();
result.feasible = all(isfinite(DV_perModule));
result.Nmodules = Nmodules;
result.methods = methods;
result.summary = summary;
result.feasibleByModuleMethod = feasible;
result.DV_byModuleMethod_mps = DV;
result.TOF_byModuleMethod_min = TOF;
result.records = records;
result.bestMethodPerModule = bestMethod;
result.bestRecords = bestRecords;
result.DV_perModule = DV_perModule;
result.DV_outbound = DV_outbound;
result.DV_return = DV_return;
result.TOF_perModule_min = TOF_perModule;
result.timePerLeg_min = timePerLeg;
result.pathsOut = pathsOut;
result.pathsReturn = pathsReturn;
result.orbitalMeanMotion_radps = n;
result.options = opts;
if result.feasible
    result.DV = sum(DV_perModule);
    result.TOF_min = sum(TOF_perModule);
else
    result.DV = NaN;
    result.TOF_min = NaN;
end

if opts.printSummary
    fprintf('\nFABER ASSEMBLY RELOCATIONS (30-m gated)\n');
    fprintf('Altitude %.2f km | Modules %d\n',h,Nmodules);
    disp(summary);
    disp(table((1:Nmodules)',bestMethod,DV_outbound,DV_return, ...
        DV_perModule,TOF_perModule,'VariableNames', ...
        {'Module','Method','DV_Out_mps','DV_Return_mps', ...
         'DV_Total_mps','TOF_Total_min'}));
    fprintf('Best-per-module total DV: %.5f m/s\n',result.DV);
    fprintf('Best-per-module total TOF: %.2f min\n',result.TOF_min);
end
end

function out = optimizeDirection(source,target,method,n,Npts,timesMin,opts)
out = struct('feasible',false);
delta = target-source;
d = norm(delta);
if d <= opts.gateDistance_m
    % Not enough separation to introduce a 30-m gate.
    return;
end

gate = chooseGate(source,target,method,opts.gateDistance_m);
if norm(gate-source)<=1e-8, return; end

% HCW transfer -> gate. Endpoint relative velocities are zero.
bestTransfer = struct( ...
    'feasible',true, ...
    'dv',sol.dvTotal, ...
    'T_min',tMin, ...
    'path',X, ...
    'maxSpeed_mps',max(vecnorm(X(4:6,:),2,1)));

for tMin = timesMin(:)'
    try
        sol = HCW_solver(source,zeros(3,1),gate,zeros(3,1), ...
            60*tMin,n,Npts);
        if ~isfield(sol,'dvTotal') || ~isfield(sol,'trajectory'), continue; end
        X = sol.trajectory;
        if size(X,1)<6 || any(~isfinite(X(:))) || ~isfinite(sol.dvTotal),continue;end
        if norm(X(1:3,end)-gate)>1e-3,continue;end
        if max(vecnorm(X(4:6,:),2,1))>opts.maxTransferSpeed_mps,continue;end
        % Never enter the 30-m gate sphere during the transfer.
        rho = vecnorm(X(1:3,:)-target,2,1);
        if min(rho)<opts.gateDistance_m-1e-5,continue;end
        if opts.minTargetRange_m>0 && min(rho)<opts.minTargetRange_m,continue;end
        if sol.dvTotal < bestTransfer.dv
            bestTransfer = struct('feasible',true,'dv',sol.dvTotal, ...
                'T_min',tMin,'path',X);
        end
    catch
        % Individual candidate is infeasible/numerically singular.
        continue;
    end
end
if ~bestTransfer.feasible,return;end

% Continuous final (same on all method families for a given gate).
bestFinal = struct( ...
    'feasible',true, ...
    'dv',dv, ...
    'T_min',tMin, ...
    'path',X, ...
    'maxSpeed_mps',max(vecnorm(X(4:6,:),2,1)), ...
    'maxAccel_mps2',max(vecnorm(u,2,1)));

for tMin = opts.finalTimesMin(:)'
    t = linspace(0,60*tMin,max(Npts,201));
    [X,u] = quinticHCW(gate,target,60*tMin,n,t);
    if any(~isfinite(X(:))) || any(~isfinite(u(:))),continue;end
    if max(vecnorm(X(4:6,:),2,1))>opts.maxFinalSpeed_mps,continue;end
    if max(vecnorm(u,2,1))>opts.maxFinalAccel_mps2,continue;end
    rho = vecnorm(X(1:3,:)-target,2,1);
    if opts.minTargetRange_m>0 && min(rho(1:end-1))<opts.minTargetRange_m,continue;end
    dv = trapz(t,vecnorm(u,2,1));
    if dv < bestFinal.dv
        bestFinal = struct('feasible',true,'dv',dv, ...
            'T_min',tMin,'path',X, ...
            'maxAccel_mps2',max(vecnorm(u,2,1)));
    end
end
if ~bestFinal.feasible,return;end
out.feasible = true;
out.method = string(method);
out.source_m = source;
out.target_m = target;
out.gate_m = gate;
out.transferDV_mps = bestTransfer.dv;
out.finalDV_mps = bestFinal.dv;
out.dv_mps = bestTransfer.dv+bestFinal.dv;
out.transferTOF_min = bestTransfer.T_min;
out.finalTOF_min = bestFinal.T_min;
out.tof_min = bestTransfer.T_min+bestFinal.T_min;
out.transferPath = bestTransfer.path;
out.finalPath = bestFinal.path;
out.finalMaxAccel_mps2 = bestFinal.maxAccel_mps2;

out.transferMaxSpeed_mps = bestTransfer.maxSpeed_mps;
out.finalMaxSpeed_mps = bestFinal.maxSpeed_mps;

end

function gate = chooseGate(source,target,method,gateDistance)
delta = source-target;
switch lower(method)
    case 'generalhcw'
        direction = delta/norm(delta);
    case 'vbar'
        direction = [0;signOrPlus(delta(2));0];
    case 'rbar'
        direction = [signOrPlus(delta(1));0;0];
    otherwise
        error('Unknown approach geometry: %s',method);
end
gate = target+gateDistance*direction;
end

function s = signOrPlus(a)
s = sign(a);
if s==0,s=1;end
end

function [X,u] = quinticHCW(r0,rf,T,n,t)
% Rest-to-rest quintic, with accelerations satisfying controlled HCW.
s = reshape(t,1,[])/T;
D = rf(:)-r0(:);
q = 10*s.^3-15*s.^4+6*s.^5;
qd = (30*s.^2-60*s.^3+30*s.^4)/T;
qdd = (60*s-180*s.^2+120*s.^3)/T^2;
r = r0(:)+D*q;
v = D*qd;
a = D*qdd;
X = [r;v];
u = [a(1,:)-3*n^2*r(1,:)-2*n*v(2,:); ...
     a(2,:)+2*n*v(1,:); ...
     a(3,:)+n^2*r(3,:)];
end

function opts = defaultOpts(opts)
defs = struct('methods',{{'GeneralHCW','Vbar','Rbar'}}, ...
    'gateDistance_m',30,'finalTimesMin',5:5:30, ...
    'maxTransferSpeed_mps',Inf,'maxFinalSpeed_mps',Inf, ...
    'maxFinalAccel_mps2',Inf,'minTargetRange_m',0, ...
    'printSummary',true);
f = fieldnames(defs);
for k=1:numel(f)
    if ~isfield(opts,f{k}),opts.(f{k})=defs.(f{k});end
end
end
