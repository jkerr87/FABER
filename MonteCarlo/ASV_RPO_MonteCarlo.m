
function out = ASV_RPO_MonteCarlo(cfg)
% ASV_RPO_MONTECARLO
% FABER RPO design-space Monte Carlo (client + host).
%
% Usage:
%   cfg = ASV_RPO_MC_config();
%   out = ASV_RPO_MonteCarlo(cfg);
%
% Phase sequence:
%   Homing -> Holding 1 -> Closing -> Holding 2 -> Final
%
% Homing and closing:
%   Randomized multiple two-impulse HCW transfers.
%
% Final approach:
%   Continuous quintic reference trajectory with HCW feedforward.
%
% Monte Carlo design variables:
%   - HP1 and HP2 locations
%   - Homing, closing and final approach TOF
%   - Holding durations
%   - Number of homing and closing transfers
%   - Intermediate waypoint locations
%   - Time allocation among transfers
%
% IMPORTANT:
%   This is Monte Carlo sampling, NOT optimization.
%   Each attempt uses one randomly sampled design per approach.
%   TOF is sampled independently of the orbital period.
%
% P95 describes accepted feasible sampled designs.
% It does NOT represent mission reliability.

%% 1. INITIALIZATION

validateCfg(cfg);

rng(cfg.seed, 'twister');

% Reference orbital mean motion [rad/s]
n = sqrt(cfg.mu_km3s2 / ...
    (cfg.Re_km + cfg.h_km)^3);

phaseNames = {'Homing','Holding 1','Closing', ...
              'Holding 2','Final approach'};

labels = [strcat('Client - ',phaseNames), ...
          strcat('Host - ',phaseNames)];

maxAttempts = max(cfg.N, ...
    ceil(cfg.N * cfg.maxAttemptsFactor));

DV_phase = nan(cfg.N,10);
TOF_phase = nan(cfg.N,10);

% Nine design variables stored per approach:
% HP1, HP2, Homing TOF, Hold 1 TOF, Closing TOF,
% Hold 2 TOF, Final TOF,
% Number of homing transfers, Number of closing transfe
settings = nan(cfg.N,18);

dv_client = nan(cfg.N,1);
dv_host = nan(cfg.N,1);

tof_client = nan(cfg.N,1);
tof_host = nan(cfg.N,1);

methodClient = strings(cfg.N,1);
methodHost = strings(cfg.N,1);

peakSpeed = nan(cfg.N,1);

X0_client = nan(cfg.N,6);
X0_host   = nan(cfg.N,6);

fail = struct( ...
    'geometry',0, ...
    'time',0, ...
    'terminal',0, ...
    'singular',0, ...
    'KOZ',0, ...
    'excursion',0, ...
    'speed_homing',0, ...
    'speed_closing',0, ...
    'speed_final',0, ...
    'acceleration',0, ...
    'dv_limit',0, ...
    'numeric',0);

%% 2. MONTE CARLO LOOP

accepted = 0;
attempts = 0;

while accepted < cfg.N && attempts < maxAttempts

    attempts = attempts + 1;

    % Sample initial conditions (NEW)

    x0_client = sampleInitialState( ...
        cfg.x0_client, cfg.dx0_client);

    x0_host = sampleInitialState( ...
        cfg.x0_host, cfg.dx0_host);

    % One randomly sampled client approach
    [c,reasonC] = sampleLeg(cfg,x0_client,n);  % MODIFIED

    if ~isempty(reasonC)
        fail.(reasonC) = fail.(reasonC) + 1;
        continue;
    end

    % One randomly sampled host approach
    [h,reasonH] = sampleLeg(cfg,x0_host,n);    % MODIFIED

    if ~isempty(reasonH)
        fail.(reasonH) = fail.(reasonH) + 1;
        continue;
    end

    % Total cycle Delta V constraint
    if c.dvTotal + h.dvTotal > cfg.maxDV_perCycle_mps
        fail.dv_limit = fail.dv_limit + 1;
        continue;
    end

    % Accept feasible realization
    accepted = accepted + 1;

    % Store initial conditions for accepted trial
    X0_client(accepted,:) = x0_client';
    X0_host(accepted,:)   = x0_host';

    DV_phase(accepted,:) = [c.dvPhase h.dvPhase];
    TOF_phase(accepted,:) = [c.tofPhase h.tofPhase];

    settings(accepted,:) = [c.draw h.draw];

    dv_client(accepted) = c.dvTotal;
    dv_host(accepted) = h.dvTotal;

    tof_client(accepted) = c.tofTotal;
    tof_host(accepted) = h.tofTotal;

    peakSpeed(accepted) = ...
        max(c.peakSpeed,h.peakSpeed);

    methodClient(accepted) = string(c.method);
    methodHost(accepted) = string(h.method);

end

%% 3. FEASIBILITY DIAGNOSTICS

if accepted == 0

    fprintf('\n========================================\n');
    fprintf('MONTE CARLO FEASIBILITY DIAGNOSTICS\n');
    fprintf('========================================\n');

    fprintf('Attempted cycles: %d\n',attempts);
    fprintf('Accepted cycles: %d\n',accepted);

    fprintf('\nRejection counts:\n');
    disp(fail);

    error(['No feasible trials after %d attempts. ', ...
        'Check rejection counts above.'],attempts);

end

if accepted < cfg.N
    warning( ...
        'Only %d/%d samples accepted after %d attempts.', ...
        accepted,cfg.N,attempts);
end

%% 4. TRIM ARRAYS

DV_phase = DV_phase(1:accepted,:);
TOF_phase = TOF_phase(1:accepted,:);
settings = settings(1:accepted,:);

dv_client = dv_client(1:accepted);
dv_host = dv_host(1:accepted);

tof_client = tof_client(1:accepted);
tof_host = tof_host(1:accepted);

peakSpeed = peakSpeed(1:accepted);

methodClient = methodClient(1:accepted);
methodHost = methodHost(1:accepted);

dv_cycle = dv_client + dv_host;
tof_cycle = tof_client + tof_host;

X0_client = X0_client(1:accepted,:);
X0_host   = X0_host(1:accepted,:);

%% PARETO FRONTIER: MINIMIZE TOF AND DELTA V

% Sort feasible samples by increasing TOF, then Delta V
[~,sortIdx] = sortrows([tof_cycle,dv_cycle],[1 2]);

paretoMask = false(accepted,1);

bestDV = inf;

for j = 1:accepted

    i = sortIdx(j);

    if dv_cycle(i) < bestDV
        paretoMask(i) = true;
        bestDV = dv_cycle(i);
    end

end

% Extract Pareto designs, ordered by TOF
paretoIdx = sortIdx(paretoMask(sortIdx));

paretoTOF = tof_cycle(paretoIdx);
paretoDV = dv_cycle(paretoIdx);

%% 5. PHASE STATISTICS

phaseMeanDV = mean(DV_phase,1)';
phaseP95DV = colPercentile(DV_phase,95)';

phaseMeanTOF = mean(TOF_phase,1)';
phaseP95TOF = colPercentile(TOF_phase,95)';

phaseTable = table( ...
    string(labels(:)), ...
    phaseMeanDV,phaseP95DV, ...
    phaseMeanTOF,phaseP95TOF, ...
    'VariableNames', ...
    {'Phase','MeanDV_mps','P95DV_mps', ...
     'MeanTOF_min','P95TOF_min'});

%% 6. TOTAL CYCLE STATISTICS

summaryLabels = [ ...
    "Client RPO"; ...
    "Host RPO"; ...
    "Complete retrieval RPO"];

DVs = [dv_client dv_host dv_cycle];
Ts = [tof_client tof_host tof_cycle];

summaryTable = table( ...
    summaryLabels, ...
    mean(DVs,1)', ...
    colPercentile(DVs,50)', ...
    colPercentile(DVs,95)', ...
    max(DVs,[],1)', ...
    mean(Ts,1)', ...
    colPercentile(Ts,50)', ...
    colPercentile(Ts,95)', ...
    max(Ts,[],1)', ...
    'VariableNames', ...
    {'Scope','MeanDV_mps','P50DV_mps', ...
     'P95DV_mps','MaxDV_mps', ...
     'MeanTOF_min','P50TOF_min', ...
     'P95TOF_min','MaxTOF_min'});

%% 7. MONTE CARLO TRIAL TABLE

trialTable = table( ...
    (1:accepted)', ...
    dv_client,dv_host,dv_cycle, ...
    tof_client,tof_host,tof_cycle, ...
    peakSpeed,methodClient,methodHost, ...
    'VariableNames', ...
    {'Trial','ClientDV_mps','HostDV_mps', ...
     'CycleDV_mps','ClientTOF_min', ...
     'HostTOF_min','CycleTOF_min', ...
     'PeakSpeed_mps','ClientMethod','HostMethod'});

drawNames = { ...
    'ClientHP1_m', ...
    'ClientHP2_m', ...
    'ClientHoming_min', ...
    'ClientHold1_min', ...
    'ClientClosing_min', ...
    'ClientHold2_min', ...
    'ClientFinal_min', ...
    'ClientN_Homing', ...
    'ClientN_Closing', ...
    'HostHP1_m', ...
    'HostHP2_m', ...
    'HostHoming_min', ...
    'HostHold1_min', ...
    'HostClosing_min', ...
    'HostHold2_min', ...
    'HostFinal_min', ...
    'HostN_Homing', ...
    'HostN_Closing'};

trialTable = [trialTable ...
    array2table(settings,'VariableNames',drawNames)];

initialNames = { ...
    'ClientX0_m', 'ClientY0_m', 'ClientZ0_m', ...
    'ClientVx0_mps', 'ClientVy0_mps', 'ClientVz0_mps', ...
    'HostX0_m', 'HostY0_m', 'HostZ0_m', ...
    'HostVx0_mps', 'HostVy0_mps', 'HostVz0_mps'};

trialTable = [trialTable ...
    array2table([X0_client X0_host], ...
    'VariableNames',initialNames)];

%% SAVE PARETO DESIGN PARAMETERS

paretoTable = trialTable(paretoIdx,:);


%% 8. PRINT RESULTS

fprintf('\n============================================\n');
fprintf('FABER RPO MONTE CARLO\n');
fprintf('%d/%d accepted (%.1f%%)\n', ...
    accepted,attempts,100*accepted/attempts);
fprintf('============================================\n');

disp(summaryTable);
disp(phaseTable);

fprintf('\nFeasibility rejection counts:\n');
disp(fail);

[~,idxMinDV] = min(dv_cycle);
[~,idxMinTOF] = min(tof_cycle);

fprintf( ...
    'Minimum sampled DV: %.5f m/s, TOF %.2f min\n', ...
    dv_cycle(idxMinDV),tof_cycle(idxMinDV));

fprintf( ...
    'Minimum sampled TOF: %.2f min, DV %.5f m/s\n', ...
    tof_cycle(idxMinTOF),dv_cycle(idxMinTOF));

fprintf('One random design per approach.\n');
fprintf('No objective-based optimization performed.\n');

fprintf('NOTE: P95 describes feasible sampled designs.\n');

%% 9. OUTPUT STRUCTURE

out = struct();

out.paretoTable = paretoTable;
out.paretoIndices = paretoIdx;
out.paretoDV_mps = paretoDV;
out.paretoTOF_min = paretoTOF;

out.cfg = cfg;
out.accepted = accepted;
out.attempts = attempts;

out.failCounts = fail;

out.summaryTable = summaryTable;
out.phaseTable = phaseTable;
out.trialTable = trialTable;

out.DV_phase_mps = DV_phase;
out.TOF_phase_min = TOF_phase;

out.bestDV_trial = idxMinDV;
out.bestTOF_trial = idxMinTOF;

out.P95_DV_cycle_mps = percentile(dv_cycle,95);
out.P95_TOF_cycle_min = percentile(tof_cycle,95);

out.selectedClientMethod = methodClient;
out.selectedHostMethod = methodHost;

out.jointP95Exceedance = mean( ...
    dv_cycle <= out.P95_DV_cycle_mps & ...
    tof_cycle <= out.P95_TOF_cycle_min);

%% 10. SAVE RESULTS

if cfg.saveCSV || cfg.makePlots
    if ~exist(cfg.outputFolder,'dir')
        mkdir(cfg.outputFolder);
    end
end

if cfg.saveCSV

    writetable(trialTable, ...
        fullfile(cfg.outputFolder,'RPO_MC_trials.csv'));

    writetable(summaryTable, ...
        fullfile(cfg.outputFolder,'RPO_MC_summary.csv'));

    writetable(phaseTable, ...
        fullfile(cfg.outputFolder,'RPO_MC_phases.csv'));

    writetable(paretoTable, ...
        fullfile(cfg.outputFolder,'RPO_MC_Pareto.csv'));

end

if cfg.makePlots

    makePlots(cfg,dv_cycle,tof_cycle, ...
        DV_phase,TOF_phase,labels, ...
        out.P95_DV_cycle_mps, ...
        out.P95_TOF_cycle_min, ...
        idxMinDV,idxMinTOF);

end

end

%% ================================================================
% LOCAL FUNCTION: RANDOM APPROACH
% ================================================================

function [trial,reason] = sampleLeg(cfg,x0,n)

% Generate exactly ONE random design.
% No best-of-K selection or objective minimization.

d = sampleDesign(cfg);

% Randomly choose an allowed method
method = char(cfg.methods{randi(numel(cfg.methods))});

[trial,reason] = simulateLeg(cfg,x0,n,d,method);

end

%% ================================================================
% LOCAL FUNCTION: SAMPLE DESIGN VARIABLES
% ================================================================

function d = sampleDesign(cfg)
% SAMPLE DESIGN
% Randomized architecture and TOF.
% No objective-based optimization.

%% Geometry [m]
d.hp1 = draw(cfg.hp1_m);
d.hp2 = draw(cfg.hp2_m);
d.final_m = draw(cfg.final_m);

%% Phase TOFs [s]
% Independent of orbital period
d.Th = 60*draw(cfg.homing_min);
d.Tc = 60*draw(cfg.closing_min);
d.Tfinal = 60*draw(cfg.final_min);

%% Hold durations [min]
d.T_hold1 = draw(cfg.hold1_min);
d.T_hold2 = draw(cfg.hold2_min);

%% Number of HCW arcs
d.Nh = randi(cfg.N_homing);
d.Nc = randi(cfg.N_closing);

%% Minimum time per HCW arc [s]
% Configurable, default 120 seconds
if isfield(cfg,'minTransferTime_s')
    minLeg = cfg.minTransferTime_s;
else
    minLeg = 120;
end

if ~isscalar(minLeg) || ...
        ~isfinite(minLeg) || minLeg <= 0
    error('minTransferTime_s must be positive.');
end

%% Feasibility of sampled phase durations
d.validTime = ...
    d.Th >= d.Nh*minLeg && ...
    d.Tc >= d.Nc*minLeg;

if ~d.validTime
    d.homingFrac = [];
    d.closingFrac = [];
    d.homingLegTimes_s = [];
    d.closingLegTimes_s = [];
    return;
end

%% Intermediate waypoint fractions
% Strictly increasing with final fraction = 1
d.homingFrac = sort([rand(1,d.Nh-1),1]);
d.closingFrac = sort([rand(1,d.Nc-1),1]);

%% Segment fractions
hSeg = diff([0,d.homingFrac]);
cSeg = diff([0,d.closingFrac]);

%% Randomized time allocation
% Longer spatial segments receive more time,
% with +/-20% random weighting.
hVariation = 0.8 + 0.4*rand(1,d.Nh);
cVariation = 0.8 + 0.4*rand(1,d.Nc);

hWeights = hSeg.*hVariation;
cWeights = cSeg.*cVariation;

hWeights = hWeights/sum(hWeights);
cWeights = cWeights/sum(cWeights);

%% Allocate available time
% Every leg receives at least minLeg seconds.
d.homingLegTimes_s = minLeg + ...
    (d.Th-d.Nh*minLeg)*hWeights;

d.closingLegTimes_s = minLeg + ...
    (d.Tc-d.Nc*minLeg)*cWeights;

%% Numerical verification
tol = 1e-8;

d.validTime = ...
    all(d.homingLegTimes_s >= minLeg-tol) && ...
    all(d.closingLegTimes_s >= minLeg-tol) && ...
    abs(sum(d.homingLegTimes_s)-d.Th) < tol && ...
    abs(sum(d.closingLegTimes_s)-d.Tc) < tol;

end


%% ================================================================
% LOCAL FUNCTION: SIMULATE ONE APPROACH
% ================================================================


function [res,fail,nSVD] = simulateLeg(cfg,x0,n,d,method)

res = [];
fail = '';
nSVD = 0;

if ~d.validTime
    fail = 'time';
    return;
end

hp1 = d.hp1;
hp2 = d.hp2;
final_m = d.final_m;

Th = d.Th;
Tc = d.Tc;
Tfinal = d.Tfinal;

T_hold1 = d.T_hold1;
T_hold2 = d.T_hold2;

%% GEOMETRY

if ~(hp1>hp2 && ...
     hp2>cfg.KOZ_m && ...
     hp2<norm(x0(1:3)) && ...
     final_m>0 && ...
     final_m<cfg.KOZ_m && ...
     cfg.AZ_m>cfg.KOZ_m && ...
     cfg.AZ_m<hp1 && ...
     cfg.AZ_m>hp2 && ...
     cfg.corridorHalfWidth_m>0)

    fail = 'geometry';
    return;
end

%% ================================================================
% HOMING: INTERMEDIATE FLY-THROUGH WAYPOINTS
% ================================================================

x_start = x0;
r_initial = x0(1:3);
r_final = [0;-hp1;0];

dh = 0;
spH = 0;
minH = inf;

for k = 1:d.Nh

    alpha = d.homingFrac(k);

    r_target = r_initial + ...
        alpha*(r_final-r_initial);

    % Final homing waypoint = HP1 (full stop)
    stopAtArrival = (k == d.Nh);

    % Desired position is always prescribed.
    % Intermediate arrival velocity is free.
    % HP1 terminal velocity is zero.
    x_target = [r_target;zeros(3,1)];

    T_leg = d.homingLegTimes_s(k);

    [sol,reason,usedSVD] = safeHCWTransfer( ...
        x_start,x_target,T_leg,n, ...
        cfg.nTrajectorySamples,stopAtArrival);

    nSVD = nSVD+usedSVD;

    if ~isempty(reason)
        fail = reason;
        return;
    end

    % First burn for each arc.
    % Second burn included only at HP1.
    dh = dh+sol.dvTotal;

    spH = max(spH,sol.maxSpeed);
    minH = min(minH,sol.minRange);

    % Carry actual arrival velocity forward.
    % At intermediate waypoints it is nonzero.
    % At HP1 it is zero after braking.
    x_start = sol.finalState;

end

%% ================================================================
% APPROACH GEOMETRY
% ================================================================

if strcmpi(method,'Vbar')

    targetHP2 = [0;-hp2;0;0;0;0];

    rStart = [0;-hp2;0];
    rEnd = [0;-final_m;0];

    dvRadialHold2 = 0;

elseif strcmpi(method,'Rbar')

    targetHP2 = [hp2;0;0;0;0;0];

    rStart = [hp2;0;0];
    rEnd = [final_m;0;0];

    % Ideal radial stationkeeping
    dvRadialHold2 = ...
        3*n^2*hp2*(T_hold2*60);

else
    error('Unknown approach method.');
end

%% ================================================================
% CLOSING: INTERMEDIATE FLY-THROUGH WAYPOINTS
% ================================================================

% Start from HP1 at rest.
x_start = [0;-hp1;0;0;0;0];

r_initial = x_start(1:3);
r_final = targetHP2(1:3);

dcl = 0;
spC = 0;
minC = inf;
maxC = 0;

for k = 1:d.Nc

    alpha = d.closingFrac(k);

    r_target = r_initial + ...
        alpha*(r_final-r_initial);

    % Final closing waypoint = HP2 (full stop)
    stopAtArrival = (k == d.Nc);

    x_target = [r_target;zeros(3,1)];

    T_leg = d.closingLegTimes_s(k);

    [sol,reason,usedSVD] = safeHCWTransfer( ...
        x_start,x_target,T_leg,n, ...
        cfg.nTrajectorySamples,stopAtArrival);

    nSVD = nSVD+usedSVD;

    if ~isempty(reason)
        fail = reason;
        return;
    end

    dcl = dcl+sol.dvTotal;

    spC = max(spC,sol.maxSpeed);
    minC = min(minC,sol.minRange);
    maxC = max(maxC,sol.maxRange);

    % Preserve fly-through velocity between arcs.
    x_start = sol.finalState;

end

%% ================================================================
% FINAL APPROACH
% ================================================================

[dFinal,spF,accF] = finalQuinticAxis( ...
    rStart,rEnd,Tfinal,n, ...
    max(cfg.nTrajectorySamples,101));

peak = max([spH spC spF]);

if ~all(isfinite( ...
        [dh dcl dFinal peak minH minC maxC accF]))

    fail = 'numeric';
    return;
end

%% ================================================================
% SAFETY
% ================================================================

% KOZ exclusion before final approach
if min(minH,minC) < cfg.KOZ_m
    fail = 'KOZ';
    return;
end

if spH > cfg.maxSpeed_homing_mps
    fail = 'speed_homing';
    return;
end

if spC > cfg.maxSpeed_closing_mps
    fail = 'speed_closing';
    return;
end

if spF > cfg.maxSpeed_final_mps
    fail = 'speed_final';
    return;
end

if accF > cfg.maxPeakAccel_mps2
    fail = 'acceleration';
    return;
end

%% ================================================================
% OUTPUT
% ================================================================

res.dvPhase = [ ...
    dh, ...
    T_hold1*cfg.hold1_DVrate_mps_per_min, ...
    dcl, ...
    T_hold2*cfg.hold2_DVrate_mps_per_min + ...
        dvRadialHold2, ...
    dFinal];

res.tofPhase = [ ...
    Th/60, ...
    T_hold1, ...
    Tc/60, ...
    T_hold2, ...
    Tfinal/60];

res.dvTotal = sum(res.dvPhase);
res.tofTotal = sum(res.tofPhase);

res.peakSpeed = peak;
res.method = method;

res.draw = [ ...
    hp1 hp2 Th/60 T_hold1 Tc/60 T_hold2 Tfinal/60 ...
    d.Nh d.Nc];

end



%% ================================================================
% LOCAL FUNCTION: CONTINUOUS FINAL APPROACH
% ================================================================

function [dv,maxSpeed,peakAccel] = ...
    finalQuinticAxis(rStart,rEnd,T,n,N)

t = linspace(0,T,N);
s = t/T;

% Quintic position profile
p = 10*s.^3-15*s.^4+6*s.^5;

% First derivative
pd = (30*s.^2-60*s.^3+30*s.^4)/T;

% Second derivative
pdd = (60*s-180*s.^2+120*s.^3)/T^2;

dr = rEnd-rStart;

r = rStart+dr*p;
v = dr*pd;
a = dr*pdd;

% Exact HCW feedforward acceleration
u = [ ...
    a(1,:)-3*n*n*r(1,:)-2*n*v(2,:); ...
    a(2,:)+2*n*v(1,:); ...
    a(3,:)+n*n*r(3,:)];

aMag = sqrt(sum(u.^2,1));

dv = trapz(t,aMag);

maxSpeed = max(sqrt(sum(v.^2,1)));
peakAccel = max(aMag);

end

%% ================================================================
% LOCAL FUNCTION: UNIFORM RANDOM SAMPLE
% ================================================================

function x = draw(bounds)

x = bounds(1)+(bounds(2)-bounds(1))*rand();

end

%% ================================================================
% LOCAL FUNCTION: PERCENTILE
% ================================================================

function p = percentile(x,q)

x = sort(x(:));
n = numel(x);

if n == 1
    p = x;
    return;
end

r = 1+(n-1)*q/100;

lo = floor(r);
hi = ceil(r);

p = x(lo)+(r-lo)*(x(hi)-x(lo));

end

function vals = colPercentile(X,q)

vals = zeros(1,size(X,2));

for j = 1:size(X,2)
    vals(j) = percentile(X(:,j),q);
end

end

%% =====================================================================
% fn to give the MC a new initial state each run
% ======================================================================
function x0 = sampleInitialState(x0_nom, dx0)

% Uniform variations around nominal initial state
% x0_nom: nominal 6x1 relative state
% dx0: maximum absolute variation in each component

delta = (2*rand(6,1)-1).*dx0;

x0 = x0_nom + delta;

end

%% ================================================================
% LOCAL FUNCTION: CONFIGURATION VALIDATION
% ================================================================

function validateCfg(cfg)

assert(isstruct(cfg), ...
    'Input must be a configuration struct.');

assert(cfg.N>=1 && cfg.N==floor(cfg.N), ...
    'N must be a positive integer.');

assert(cfg.nTrajectorySamples>=3, ...
    'Need at least 3 trajectory samples.');

assert(cfg.maxAttemptsFactor>=1);

assert(cfg.h_km>100 && ...
       cfg.Re_km>0 && ...
       cfg.mu_km3s2>0);

assert(all(size(cfg.x0_client)==[6 1]) && ...
       all(size(cfg.x0_host)==[6 1]));

% Maneuver-count bounds
assert( ...
    numel(cfg.N_homing)==2 && ...
    all(cfg.N_homing>=1) && ...
    all(cfg.N_homing==floor(cfg.N_homing)) && ...
    cfg.N_homing(2)>=cfg.N_homing(1), ...
    'Invalid homing transfer count.');

assert( ...
    numel(cfg.N_closing)==2 && ...
    all(cfg.N_closing>=1) && ...
    all(cfg.N_closing==floor(cfg.N_closing)) && ...
    cfg.N_closing(2)>=cfg.N_closing(1), ...
    'Invalid closing transfer count.');

keys = { ...
    'hp1_m','hp2_m','final_m', ...
    'homing_min','closing_min', ...
    'final_min','hold1_min','hold2_min'};

for k = 1:numel(keys)

    b = cfg.(keys{k});

    assert( ...
        isnumeric(b) && ...
        numel(b)==2 && ...
        all(isfinite(b)) && ...
        b(1)>0 && b(2)>=b(1), ...
        'Invalid bounds for %s.',keys{k});

end

assert(cfg.KOZ_m>0 && ...
       cfg.AZ_m>cfg.KOZ_m);

assert(cfg.maxSpeed_homing_mps>0 && ...
       cfg.maxSpeed_closing_mps>0 && ...
       cfg.maxSpeed_final_mps>0);

assert(cfg.maxPeakAccel_mps2>0);

assert(iscell(cfg.methods) && ...
       ~isempty(cfg.methods));

assert(all(ismember( ...
    lower(string(cfg.methods)),["vbar","rbar"])));

assert(all(size(cfg.dx0_client)==[6 1]) && ...
    all(size(cfg.dx0_host)==[6 1]), ...
    'Initial-condition variations must be 6x1.');

assert(all(isfinite(cfg.dx0_client)) && ...
    all(isfinite(cfg.dx0_host)) && ...
    all(cfg.dx0_client>=0) && ...
    all(cfg.dx0_host>=0), ...
    'Initial-condition variation bounds must be finite and nonnegative.');

end

%% ================================================================
% LOCAL FUNCTION: PLOTS
% ================================================================

function makePlots(cfg,dv,tof,DVphase,TOFphase, ...
    labels,p95dv,p95tof,iDV,iTOF)

%% DELTA V HISTOGRAM

f1 = figure('Color','w', ...
    'Name','RPO Monte Carlo - Delta V');

histogram(dv,40);
hold on;
grid on;

xline(p95dv,'--','P95','LineWidth',1.6);

xlabel('Client + Host RPO Delta V [m/s]');
ylabel('Accepted designs');

title(sprintf( ...
    'RPO Delta V (N=%d, P95=%.3f m/s)', ...
    numel(dv),p95dv));

saveas(f1,fullfile( ...
    cfg.outputFolder,'RPO_MC_DV_histogram.png'));

%% TOF HISTOGRAM

f2 = figure('Color','w', ...
    'Name','RPO Monte Carlo - TOF');

histogram(tof/60,40);
hold on;
grid on;

xline(p95tof/60,'--','P95','LineWidth',1.6);

xlabel('Client + Host RPO TOF [hours]');
ylabel('Accepted designs');

title(sprintf( ...
    'RPO TOF (P95=%.2f h)',p95tof/60));

saveas(f2,fullfile( ...
    cfg.outputFolder,'RPO_MC_TOF_histogram.png'));

%% TOF VS DELTA V

f3 = figure('Color','w', ...
    'Name','RPO Monte Carlo - TOF vs Delta V');

scatter(tof/60,dv,11,'filled');
hold on;
grid on;

plot(tof(iDV)/60,dv(iDV),'kp', ...
    'MarkerSize',13,'MarkerFaceColor','k');

plot(tof(iTOF)/60,dv(iTOF),'ks', ...
    'MarkerSize',9,'LineWidth',2);

xline(p95tof/60,'--','P95 TOF');
yline(p95dv,'--','P95 DV');

xlabel('Complete retrieval RPO TOF [hours]');
ylabel('Complete retrieval RPO Delta V [m/s]');

legend('Feasible trials','Min DV','Min TOF', ...
    'Location','best');

title('TOF-DV trade space');

saveas(f3,fullfile( ...
    cfg.outputFolder,'RPO_MC_TOF_vs_DV.png'));

%% PHASE DELTA V

f4 = figure('Color','w', ...
    'Name','RPO Monte Carlo - Phase P95 DV');

barh(colPercentile(DVphase,95));
grid on;

set(gca, ...
    'YTick',1:numel(labels), ...
    'YTickLabel',labels, ...
    'YDir','reverse');

xlabel('Per-phase P95 Delta V [m/s]');
title('Individual phase P95 (do not sum)');

saveas(f4,fullfile( ...
    cfg.outputFolder,'RPO_MC_Phase_P95_DV.png'));

%% PHASE TOF

f5 = figure('Color','w', ...
    'Name','RPO Monte Carlo - Phase P95 TOF');

barh(colPercentile(TOFphase,95));
grid on;

set(gca, ...
    'YTick',1:numel(labels), ...
    'YTickLabel',labels, ...
    'YDir','reverse');

xlabel('Per-phase P95 TOF [min]');
title('Individual phase P95 (do not sum)');

saveas(f5,fullfile( ...
    cfg.outputFolder,'RPO_MC_Phase_P95_TOF.png'));

%% TOF VS DV WITH PARETO FRONTIER

f3 = figure('Color','w', ...
    'Name','FABER RPO - Pareto Frontier');

% Feasible Monte Carlo designs
scatter(tof/60,dv,12, ...
    [0.65 0.65 0.65], ...
    'filled','MarkerFaceAlpha',0.35);

hold on;
grid on;

% Calculate nondominated designs
[~,sortIdx] = sortrows([tof(:),dv(:)],[1 2]);

paretoMask = false(numel(dv),1);
bestDV = inf;

for j = 1:numel(sortIdx)

    i = sortIdx(j);

    if dv(i) < bestDV
        paretoMask(i) = true;
        bestDV = dv(i);
    end

end

paretoIdx = sortIdx(paretoMask(sortIdx));

% Plot Pareto frontier
plot(tof(paretoIdx)/60,dv(paretoIdx), ...
    '-r','LineWidth',2.2, ...
    'Marker','o','MarkerSize',5, ...
    'MarkerFaceColor','r');

% Minimum DV design
plot(tof(iDV)/60,dv(iDV), ...
    'kp','MarkerSize',13, ...
    'MarkerFaceColor','k');

% Minimum TOF design
plot(tof(iTOF)/60,dv(iTOF), ...
    'ks','MarkerSize',10, ...
    'LineWidth',2);

xlabel('Complete retrieval RPO TOF [hours]');
ylabel('Complete retrieval Delta V [m/s]');

title('FABER RPO Monte Carlo - Delta V / TOF Pareto Frontier');

legend( ...
    'Feasible Monte Carlo designs', ...
    'Sampled Pareto frontier', ...
    'Minimum Delta V', ...
    'Minimum TOF', ...
    'Location','best');

saveas(f3,fullfile( ...
    cfg.outputFolder,'RPO_MC_Pareto_Frontier.png'));

end


function [sol,reason,usedSVD] = safeHCWTransfer( ...
    x0,xf,T,n,N,stopAtArrival)
% SAFEHCWTRANSFER
%
% stopAtArrival = true:
%   Two-impulse transfer, final velocity prescribed.
%
% stopAtArrival = false:
%   Position-constrained fly-through.
%   No terminal braking impulse.
%   Natural coast arrival velocity is retained.
%
% Uses HCW_solver.m for well-conditioned transfers.
% Uses SVD fallback for ill-conditioned transfers.

sol = [];
reason = '';
usedSVD = false;

if ~isfinite(T) || T<=0 || N<2
    reason = 'time';
    return;
end

%% HCW state matrix

A = [ ...
    0 0 0 1 0 0; ...
    0 0 0 0 1 0; ...
    0 0 0 0 0 1; ...
    3*n^2 0 0 0 2*n 0; ...
    0 0 0 -2*n 0 0; ...
    0 0 -n^2 0 0 0];

Phi = expm(A*T);

Prr = Phi(1:3,1:3);
Prv = Phi(1:3,4:6);
Pvr = Phi(4:6,1:3);
Pvv = Phi(4:6,4:6);

b = xf(1:3)-Prr*x0(1:3);

%% Matrix conditioning

rcondLimit = 1e-10;
rc = rcond(Prv);

if ~isfinite(rc)
    reason = 'numeric';
    return;
end

%% ================================================================
% WELL-CONDITIONED TRANSFER
% ================================================================

if rc >= rcondLimit

    raw = HCW_solver( ...
        x0(1:3),x0(4:6), ...
        xf(1:3),xf(4:6),T,n,N);

    if ~raw.ok
        reason = 'numeric';
        return;
    end

    traj = raw.trajectory;

    % Coast arrival velocity, before terminal burn
    vArrival = traj(4:6,end);

    if stopAtArrival

        % Full stop at HP1 or HP2
        sol.dvTotal = raw.dvTotal;
        sol.finalState = raw.finalState;

    else

        % Fly-through: no terminal burn
        sol.dvTotal = raw.dv1;
        sol.finalState = [traj(1:3,end);vArrival];

    end

else

    %% ============================================================
    % ILL-CONDITIONED TRANSFER: SVD FALLBACK
    % ============================================================

    usedSVD = true;

    [U,S,V] = svd(Prv);

    singularValues = diag(S);

    tolSVD = max(size(Prv))* ...
        eps(max(singularValues));

    invS = zeros(size(singularValues));

    for j = 1:numel(singularValues)

        if singularValues(j)>tolSVD
            invS(j) = 1/singularValues(j);
        end

    end

    % Minimum-norm departure velocity satisfying
    % the reachable position constraints.
    vAfter1 = V*diag(invS)*U'*b;

    positionResidual = norm(Prv*vAfter1-b);

    positionTolerance = ...
        1e-6+1e-9*max(1,norm(xf(1:3)));

    if ~isfinite(positionResidual) || ...
            positionResidual>positionTolerance

        reason = 'singular';
        return;
    end

    burn1 = vAfter1-x0(4:6);

    % Arrival velocity before braking
    vArrival = Pvr*x0(1:3)+Pvv*vAfter1;

    % Propagate coast
    traj = zeros(6,N);
    Xafter1 = [x0(1:3);vAfter1];

    for k = 1:N
        tk = T*(k-1)/(N-1);
        traj(:,k) = expm(A*tk)*Xafter1;
    end

    if stopAtArrival

        burn2 = xf(4:6)-vArrival;

        sol.dvTotal = norm(burn1)+norm(burn2);

        sol.finalState = traj(:,end);
        sol.finalState(4:6) = ...
            sol.finalState(4:6)+burn2;

    else

        % Fly-through: retain incoming coast velocity
        sol.dvTotal = norm(burn1);
        sol.finalState = traj(:,end);

    end

end

%% ================================================================
% TERMINAL VERIFICATION
% ================================================================

positionError = norm( ...
    sol.finalState(1:3)-xf(1:3));

if stopAtArrival

    velocityError = norm( ...
        sol.finalState(4:6)-xf(4:6));

else

    % Fly-through velocity is unconstrained.
    % Check consistency with coast propagation instead.
    velocityError = norm( ...
        sol.finalState(4:6)-traj(4:6,end));

end

if ~isfinite(positionError) || ...
        ~isfinite(velocityError) || ...
        positionError>1e-5 || ...
        velocityError>1e-7

    reason = 'terminal';
    return;
end

if ~isfinite(sol.dvTotal) || ...
        ~all(isfinite(traj(:)))

    reason = 'numeric';
    return;
end

%% ================================================================
% TRAJECTORY METRICS
% ================================================================

r = traj(1:3,:);
v = traj(4:6,:);

ranges = vecnorm(r,2,1);
speeds = vecnorm(v,2,1);

sol.maxSpeed = max(speeds);
sol.minRange = min(ranges);
sol.maxRange = max(ranges);

sol.trajectory = traj;
sol.positionError = positionError;
sol.velocityError = velocityError;

end
