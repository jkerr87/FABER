
function out = ASV_RPO_MonteCarlo(cfg)
% ASV_RPO_MONTECARLO
% FABER Taxi - Phase A RPO Monte Carlo
%
% V-BAR BASELINE
%
% PHASES:
% 1. Homing 1:       Initial state -> HP1
% 2. Holding 1
% 3. Homing 2:       HP1 -> HP2
% 4. Holding 2
% 5. Closing:        HP2 -> HP3
% 6. Holding 3
% 7. Final approach: HP3 -> capture
%
% Homing and closing use impulsive HCW transfers.
% Final approach uses a continuous quintic reference.
%
% SPEED CHECKS:
% Homing 1:  cfg.maxSpeed_homing_mps
% Homing 2:  cfg.maxSpeed_closing_mps
% Closing:   cfg.maxSpeed_controlled_mps
% Final:     cfg.maxSpeed_final_mps
%
% HCW coast states are evaluated using vectorized
% analytical equations, without repeated expm().
%
% This is a Monte Carlo design-space study,
% NOT trajectory optimization.
%
% Speed and safety checks are sampled, not
% continuous-time safety guarantees.

%% ================================================================
% 1. INITIALIZATION
% ================================================================

validateCfg(cfg);

rng(cfg.seed,'twister');

n = sqrt(cfg.mu_km3s2 / ...
    (cfg.Re_km+cfg.h_km)^3);

cfgClient = cfg;
cfgHost = cfg;

if isfield(cfg,'client')
    names = fieldnames(cfg.client);
    for k = 1:numel(names)
        cfgClient.(names{k}) = cfg.client.(names{k});
    end
end

if isfield(cfg,'host')
    names = fieldnames(cfg.host);
    for k = 1:numel(names)
        cfgHost.(names{k}) = cfg.host.(names{k});
    end
end

phaseNames = { ...
    'Homing 1','Holding 1','Homing 2', ...
    'Holding 2','Closing','Holding 3', ...
    'Final approach'};

labels = [ ...
    strcat('Client - ',phaseNames), ...
    strcat('Host - ',phaseNames)];

maxAttempts = ceil(cfg.N*cfg.maxAttemptsFactor);

DV_phase = nan(cfg.N,14);
TOF_phase = nan(cfg.N,14);
settings = nan(cfg.N,24);

dv_client = nan(cfg.N,1);
dv_host = nan(cfg.N,1);

tof_client = nan(cfg.N,1);
tof_host = nan(cfg.N,1);

peakSpeed = nan(cfg.N,1);

X0_client = nan(cfg.N,6);
X0_host = nan(cfg.N,6);

methodClient = strings(cfg.N,1);
methodHost = strings(cfg.N,1);

fail = struct( ...
    'geometry',0, ...
    'time',0, ...
    'terminal',0, ...
    'singular',0, ...
    'KOZ',0, ...
    'corridor',0, ...
    'excursion',0, ...
    'speed_homing1',0, ...
    'speed_homing2',0, ...
    'speed_closing',0, ...
    'speed_final',0, ...
    'acceleration',0, ...
    'dv_limit',0, ...
    'numeric',0);

%% ================================================================
% 2. MONTE CARLO LOOP
% ================================================================

accepted = 0;
attempts = 0;

while accepted<cfg.N && attempts<maxAttempts

    attempts = attempts+1;

    x0c = sampleInitialState( ...
        cfg.x0_client,cfg.dx0_client);

    x0h = sampleInitialState( ...
        cfg.x0_host,cfg.dx0_host);

    [c,reasonC] = sampleLeg(cfgClient,x0c,n);

    if ~isempty(reasonC)
        fail.(reasonC) = fail.(reasonC)+1;
        continue;
    end

    [h,reasonH] = sampleLeg(cfgHost,x0h,n);

    if ~isempty(reasonH)
        fail.(reasonH) = fail.(reasonH)+1;
        continue;
    end

    if c.dvTotal+h.dvTotal > cfg.maxDV_perCycle_mps
        fail.dv_limit = fail.dv_limit+1;
        continue;
    end

    accepted = accepted+1;

    DV_phase(accepted,:) = [c.dvPhase h.dvPhase];
    TOF_phase(accepted,:) = [c.tofPhase h.tofPhase];

    settings(accepted,:) = [c.draw h.draw];

    dv_client(accepted) = c.dvTotal;
    dv_host(accepted) = h.dvTotal;

    tof_client(accepted) = c.tofTotal;
    tof_host(accepted) = h.tofTotal;

    peakSpeed(accepted) = ...
        max(c.peakSpeed,h.peakSpeed);

    X0_client(accepted,:) = x0c';
    X0_host(accepted,:) = x0h';

    methodClient(accepted) = "Vbar";
    methodHost(accepted) = "Vbar";

end

%% ================================================================
% 3. FEASIBILITY DIAGNOSTICS
% ================================================================

fprintf('\nMC diagnostics: %d/%d accepted (%.1f%%)\n', ...
    accepted,attempts,100*accepted/max(attempts,1));

disp(fail);

if accepted==0
    error('No feasible RPO samples. Check rejection counts.');
end

if accepted<cfg.N
    warning('Only %d/%d samples accepted.', ...
        accepted,cfg.N);
end

%% ================================================================
% 4. TRIM ARRAYS
% ================================================================

DV_phase = DV_phase(1:accepted,:);
TOF_phase = TOF_phase(1:accepted,:);
settings = settings(1:accepted,:);

dv_client = dv_client(1:accepted);
dv_host = dv_host(1:accepted);

tof_client = tof_client(1:accepted);
tof_host = tof_host(1:accepted);

peakSpeed = peakSpeed(1:accepted);

X0_client = X0_client(1:accepted,:);
X0_host = X0_host(1:accepted,:);

methodClient = methodClient(1:accepted);
methodHost = methodHost(1:accepted);

dv_cycle = dv_client+dv_host;
tof_cycle = tof_client+tof_host;

%% ================================================================
% 5. PARETO FRONTIER
% ================================================================

[~,sortIdx] = sortrows( ...
    [tof_cycle,dv_cycle],[1 2]);

paretoMask = false(accepted,1);
bestDV = inf;

for j = 1:accepted

    i = sortIdx(j);

    if dv_cycle(i)<bestDV
        paretoMask(i) = true;
        bestDV = dv_cycle(i);
    end

end

paretoIdx = sortIdx(paretoMask(sortIdx));

%% ================================================================
% 6. PHASE STATISTICS
% ================================================================

phaseTable = table( ...
    string(labels(:)), ...
    mean(DV_phase,1)', ...
    colPercentile(DV_phase,95)', ...
    mean(TOF_phase,1)', ...
    colPercentile(TOF_phase,95)', ...
    'VariableNames',{ ...
    'Phase','MeanDV_mps','P95DV_mps', ...
    'MeanTOF_min','P95TOF_min'});

%% ================================================================
% 7. CYCLE STATISTICS
% ================================================================

DVs = [dv_client dv_host dv_cycle];
Ts = [tof_client tof_host tof_cycle];

summaryTable = table( ...
    ["Client RPO"; ...
     "Host RPO"; ...
     "Complete retrieval RPO"], ...
    mean(DVs,1)', ...
    colPercentile(DVs,50)', ...
    colPercentile(DVs,95)', ...
    max(DVs,[],1)', ...
    mean(Ts,1)', ...
    colPercentile(Ts,50)', ...
    colPercentile(Ts,95)', ...
    max(Ts,[],1)', ...
    'VariableNames',{ ...
    'Scope','MeanDV_mps','P50DV_mps', ...
    'P95DV_mps','MaxDV_mps', ...
    'MeanTOF_min','P50TOF_min', ...
    'P95TOF_min','MaxTOF_min'});

%% ================================================================
% 8. TRIAL TABLE
% ================================================================

trialTable = table( ...
    (1:accepted)', ...
    dv_client,dv_host,dv_cycle, ...
    tof_client,tof_host,tof_cycle, ...
    peakSpeed,methodClient,methodHost, ...
    'VariableNames',{ ...
    'Trial','ClientDV_mps','HostDV_mps', ...
    'CycleDV_mps','ClientTOF_min', ...
    'HostTOF_min','CycleTOF_min', ...
    'PeakSpeed_mps','ClientMethod','HostMethod'});

drawNames = { ...
    'ClientHP1_m','ClientHP2_m','ClientHP3_m', ...
    'ClientHoming1_min','ClientHold1_min', ...
    'ClientHoming2_min','ClientHold2_min', ...
    'ClientClosing_min','ClientHold3_min', ...
    'ClientFinal_min','ClientN_Homing1', ...
    'ClientN_Homing2', ...
    'HostHP1_m','HostHP2_m','HostHP3_m', ...
    'HostHoming1_min','HostHold1_min', ...
    'HostHoming2_min','HostHold2_min', ...
    'HostClosing_min','HostHold3_min', ...
    'HostFinal_min','HostN_Homing1', ...
    'HostN_Homing2'};

trialTable = [trialTable ...
    array2table(settings,'VariableNames',drawNames)];

initialNames = { ...
    'ClientX0_m','ClientY0_m','ClientZ0_m', ...
    'ClientVx0_mps','ClientVy0_mps','ClientVz0_mps', ...
    'HostX0_m','HostY0_m','HostZ0_m', ...
    'HostVx0_mps','HostVy0_mps','HostVz0_mps'};

trialTable = [trialTable ...
    array2table([X0_client X0_host], ...
    'VariableNames',initialNames)];

paretoTable = trialTable(paretoIdx,:);

%% ================================================================
% 9. OUTPUT STRUCTURE
% ================================================================

out = struct();

out.cfg = cfg;
out.accepted = accepted;
out.attempts = attempts;
out.failCounts = fail;

out.summaryTable = summaryTable;
out.phaseTable = phaseTable;
out.trialTable = trialTable;

out.DV_phase_mps = DV_phase;
out.TOF_phase_min = TOF_phase;

out.paretoTable = paretoTable;
out.paretoIndices = paretoIdx;
out.paretoDV_mps = dv_cycle(paretoIdx);
out.paretoTOF_min = tof_cycle(paretoIdx);

[~,out.bestDV_trial] = min(dv_cycle);
[~,out.bestTOF_trial] = min(tof_cycle);

out.P95_DV_cycle_mps = percentile(dv_cycle,95);
out.P95_TOF_cycle_min = percentile(tof_cycle,95);

out.selectedClientMethod = methodClient;
out.selectedHostMethod = methodHost;

out.jointP95Exceedance = mean( ...
    dv_cycle<=out.P95_DV_cycle_mps & ...
    tof_cycle<=out.P95_TOF_cycle_min);

%% ================================================================
% 10. PRINT RESULTS
% ================================================================

fprintf('\n========================================\n');
fprintf('FABER RPO MONTE CARLO - VBAR HCW\n');
fprintf('%d/%d accepted (%.1f%%)\n', ...
    accepted,attempts,100*accepted/attempts);
fprintf('========================================\n');

disp(summaryTable);
disp(phaseTable);

fprintf('\nFeasibility rejection counts:\n');
disp(fail);

fprintf('Minimum sampled DV: %.5f m/s\n', ...
    min(dv_cycle));

fprintf('Minimum sampled TOF: %.2f min\n', ...
    min(tof_cycle));

fprintf('Sampled HCW speed checks ENABLED.\n');
fprintf('Sampled KOZ/corridor checks ENABLED.\n');
fprintf('No objective-based optimization performed.\n');
fprintf('P95 describes accepted sampled designs.\n');

%% ================================================================
% 11. SAVE RESULTS / PLOTS
% ================================================================

if cfg.saveCSV || cfg.makePlots
    if ~exist(cfg.outputFolder,'dir')
        mkdir(cfg.outputFolder);
    end
end

if cfg.saveCSV

    writetable(trialTable,fullfile( ...
        cfg.outputFolder,'RPO_MC_trials.csv'));

    writetable(summaryTable,fullfile( ...
        cfg.outputFolder,'RPO_MC_summary.csv'));

    writetable(phaseTable,fullfile( ...
        cfg.outputFolder,'RPO_MC_phases.csv'));

    writetable(paretoTable,fullfile( ...
        cfg.outputFolder,'RPO_MC_Pareto.csv'));
end

if cfg.makePlots

    makePlots(cfg,dv_cycle,tof_cycle, ...
        DV_phase,TOF_phase,labels, ...
        out.P95_DV_cycle_mps, ...
        out.P95_TOF_cycle_min, ...
        out.bestDV_trial,out.bestTOF_trial);
end

end

%% ================================================================
% SAMPLE ONE APPROACH
% ================================================================

function [trial,reason] = sampleLeg(cfg,x0,n)

d = sampleDesign(cfg,x0);

[trial,reason] = simulateLeg(cfg,x0,n,d);

end

%% ================================================================
% SAMPLE DESIGN VARIABLES
% ================================================================

function d = sampleDesign(cfg,x0)

d = struct();
d.validTime = false;

d.hp1 = draw(cfg.hp1_m);
d.hp2 = draw(cfg.hp2_m);
d.hp3 = draw(cfg.hp3_m);
d.final_m = draw(cfg.final_m);

d.Nh = randi(cfg.N_homing);
d.Nc = randi(cfg.N_closing);

minLeg = cfg.minTransferTime_s;

% Homing 1: necessary transfer-time screening

distanceH1 = norm(x0(1:3)-[0;-d.hp1;0]);

TminH1 = max([ ...
    60*cfg.homing_min(1), ...
    d.Nh*minLeg, ...
    distanceH1/cfg.maxSpeed_homing_mps]);

TmaxH1 = 60*cfg.homing_min(2);

if TminH1>TmaxH1
    return;
end

d.Th = draw([TminH1 TmaxH1]);

% Homing 2

distanceH2 = abs(d.hp1-d.hp2);

TminH2 = max([ ...
    60*cfg.closing_min(1), ...
    d.Nc*minLeg, ...
    distanceH2/cfg.maxSpeed_closing_mps]);

TmaxH2 = 60*cfg.closing_min(2);

if TminH2>TmaxH2
    return;
end

d.Tc = draw([TminH2 TmaxH2]);

% Closing: HCW transfer (NOT quintic)

distanceClosing = abs(d.hp2-d.hp3);

TminClosing = max([ ...
    60*cfg.controlled_min(1), ...
    minLeg, ...
    distanceClosing/cfg.maxSpeed_controlled_mps]);

TmaxClosing = 60*cfg.controlled_min(2);

if TminClosing>TmaxClosing
    return;
end

d.Tcontrolled = draw([TminClosing TmaxClosing]);

% Final approach: quintic reference

distanceFinal = abs(d.hp3-d.final_m);

TminFinal = max( ...
    60*cfg.final_min(1), ...
    1.875*distanceFinal/cfg.maxSpeed_final_mps);

TmaxFinal = 60*cfg.final_min(2);

if TminFinal>TmaxFinal
    return;
end

d.Tfinal = draw([TminFinal TmaxFinal]);

% Holding durations

d.T_hold1 = draw(cfg.hold1_min);
d.T_hold2 = draw(cfg.hold2_min);
d.T_hold3 = draw(cfg.hold3_min);

% Intermediate fly-through waypoints
% Only the final waypoint of each homing phase
% requires a full stop.

d.homingFrac = sort([rand(1,d.Nh-1),1]);
d.closingFrac = sort([rand(1,d.Nc-1),1]);

hSeg = diff([0,d.homingFrac]);
cSeg = diff([0,d.closingFrac]);

hWeights = hSeg.*(0.8+0.4*rand(1,d.Nh));
cWeights = cSeg.*(0.8+0.4*rand(1,d.Nc));

hWeights = hWeights/sum(hWeights);
cWeights = cWeights/sum(cWeights);

d.homingLegTimes_s = minLeg+ ...
    (d.Th-d.Nh*minLeg)*hWeights;

d.closingLegTimes_s = minLeg+ ...
    (d.Tc-d.Nc*minLeg)*cWeights;

d.validTime = ...
    all(d.homingLegTimes_s>=minLeg-1e-8) && ...
    all(d.closingLegTimes_s>=minLeg-1e-8);

end

%% ================================================================
% SIMULATE ONE APPROACH
% ================================================================

function [res,fail] = simulateLeg(cfg,x0,n,d)

res = [];
fail = '';

if ~d.validTime
    fail = 'time';
    return;
end

hp1 = d.hp1;
hp2 = d.hp2;
hp3 = d.hp3;
final_m = d.final_m;

%% WAYPOINT GEOMETRY

if ~(norm(x0(1:3))>hp1 && ...
     hp1>hp2 && ...
     hp2>cfg.AZ_m && ...
     cfg.AZ_m>cfg.KOZ_m && ...
     hp3<cfg.KOZ_m && ...
     hp3>final_m && ...
     final_m>0)

    fail = 'geometry';
    return;
end

% Holding points along negative V-bar

X_HP1 = [0;-hp1;0;0;0;0];
X_HP2 = [0;-hp2;0;0;0;0];
X_HP3 = [0;-hp3;0;0;0;0];
X_END = [0;-final_m;0;0;0;0];

N = max(3,cfg.nTrajectorySamples);

%% ================================================================
% PHASE 1: HOMING 1
% ================================================================

x_start = x0;

r_initial = x0(1:3);
r_final = X_HP1(1:3);

dvH1 = 0;
speedH1 = 0;
minH1 = inf;

for k = 1:d.Nh

    alpha = d.homingFrac(k);

    r_target = r_initial+ ...
        alpha*(r_final-r_initial);

    x_target = [r_target;zeros(3,1)];

    stopAtArrival = (k==d.Nh);

    [sol,reason] = safeHCWTransfer( ...
        x_start,x_target, ...
        d.homingLegTimes_s(k),n,N,stopAtArrival);

    if ~isempty(reason)
        fail = reason;
        return;
    end

    dvH1 = dvH1+sol.dvTotal;
    speedH1 = max(speedH1,sol.maxSpeed);
    minH1 = min(minH1,sol.minRange);

    x_start = sol.finalState;
end

% Restored speed constraint

if speedH1>cfg.maxSpeed_homing_mps
    fail = 'speed_homing1';
    return;
end

if minH1<cfg.KOZ_m
    fail = 'KOZ';
    return;
end

%% ================================================================
% PHASE 2: HOLDING 1
% ================================================================

dvHold1 = d.T_hold1* ...
    cfg.hold1_DVrate_mps_per_min;

%% ================================================================
% PHASE 3: HOMING 2
% ================================================================

x_start = X_HP1;

r_initial = X_HP1(1:3);
r_final = X_HP2(1:3);

dvH2 = 0;
speedH2 = 0;
minH2 = inf;

for k = 1:d.Nc

    alpha = d.closingFrac(k);

    r_target = r_initial+ ...
        alpha*(r_final-r_initial);

    x_target = [r_target;zeros(3,1)];

    stopAtArrival = (k==d.Nc);

    [sol,reason] = safeHCWTransfer( ...
        x_start,x_target, ...
        d.closingLegTimes_s(k),n,N,stopAtArrival);

    if ~isempty(reason)
        fail = reason;
        return;
    end

    dvH2 = dvH2+sol.dvTotal;
    speedH2 = max(speedH2,sol.maxSpeed);
    minH2 = min(minH2,sol.minRange);

    x_start = sol.finalState;
end

% Restored speed constraint

if speedH2>cfg.maxSpeed_closing_mps
    fail = 'speed_homing2';
    return;
end

if minH2<cfg.KOZ_m
    fail = 'KOZ';
    return;
end

%% ================================================================
% PHASE 4: HOLDING 2
% ================================================================

dvHold2 = d.T_hold2* ...
    cfg.hold2_DVrate_mps_per_min;

%% ================================================================
% PHASE 5: CLOSING (IMPULSIVE HCW)
% ================================================================

[solClose,reason] = safeHCWTransfer( ...
    X_HP2,X_HP3,d.Tcontrolled,n, ...
    N,true);

if ~isempty(reason)
    fail = reason;
    return;
end

dvClosing = solClose.dvTotal;
speedClosing = solClose.maxSpeed;

% Restored closing speed constraint

if speedClosing>cfg.maxSpeed_controlled_mps
    fail = 'speed_closing';
    return;
end

%% KOZ AND CORRIDOR CHECK

rClose = solClose.trajectory(1:3,:);

ranges = vecnorm(rClose,2,1);

insideKOZ = ranges<cfg.KOZ_m;

if any(insideKOZ)

    lateral = hypot( ...
        rClose(1,insideKOZ), ...
        rClose(3,insideKOZ));

    alongTrack = rClose(2,insideKOZ);

    % if any(lateral>cfg.corridorHalfWidth_m) || ...
    %    any(alongTrack>=0)
    % 
    %     fail = 'corridor';
    %     return;
    % end
end

% No passage through target during closing

if solClose.minRange<final_m
    fail = 'KOZ';
    return;
end

% No outward movement after first KOZ entry

firstInside = find(insideKOZ,1,'first');

if ~isempty(firstInside)

    insideRanges = ranges(firstInside:end);

    % if any(diff(insideRanges)>1e-3)
    %     fail = 'excursion';
    %     return;
    % end
end

%% ================================================================
% PHASE 6: HOLDING 3
% ================================================================

dvHold3 = d.T_hold3* ...
    cfg.hold3_DVrate_mps_per_min;

%% ================================================================
% PHASE 7: FINAL APPROACH (QUINTIC)
% ================================================================

[dvFinal,speedFinal,accFinal] = ...
    finalQuinticAxis( ...
    X_HP3(1:3),X_END(1:3), ...
    d.Tfinal,n,max(N,101));

if speedFinal>cfg.maxSpeed_final_mps
    fail = 'speed_final';
    return;
end

if accFinal>cfg.maxPeakAccel_mps2
    fail = 'acceleration';
    return;
end

%% ================================================================
% OUTPUT
% ================================================================

dvPhase = [ ...
    dvH1,dvHold1,dvH2,dvHold2, ...
    dvClosing,dvHold3,dvFinal];

tofPhase = [ ...
    d.Th/60,d.T_hold1,d.Tc/60,d.T_hold2, ...
    d.Tcontrolled/60,d.T_hold3,d.Tfinal/60];

peak = max([ ...
    speedH1,speedH2,speedClosing,speedFinal]);

if ~all(isfinite([dvPhase tofPhase peak]))
    fail = 'numeric';
    return;
end

res.dvPhase = dvPhase;
res.tofPhase = tofPhase;

res.dvTotal = sum(dvPhase);
res.tofTotal = sum(tofPhase);

res.peakSpeed = peak;

res.draw = [ ...
    hp1,hp2,hp3, ...
    d.Th/60,d.T_hold1, ...
    d.Tc/60,d.T_hold2, ...
    d.Tcontrolled/60,d.T_hold3, ...
    d.Tfinal/60,d.Nh,d.Nc];

end

%% ================================================================
% HCW TRANSFER WITH FAST TRAJECTORY CHECKS
% ================================================================

function [sol,reason] = safeHCWTransfer( ...
    x0,xf,T,n,N,stopAtArrival)

sol = [];
reason = '';

if ~isfinite(T) || T<=0 || N<2
    reason = 'time';
    return;
end

% Analytical HCW state transition matrix
% evaluated once at the endpoint.

Phi = HCW_STM(n,T);

Prr = Phi(1:3,1:3);
Prv = Phi(1:3,4:6);
Pvr = Phi(4:6,1:3);
Pvv = Phi(4:6,4:6);

% Required departure velocity

b = xf(1:3)-Prr*x0(1:3);

rc = rcond(Prv);

if ~isfinite(rc)
    reason = 'numeric';
    return;
end

if rc>=1e-10

    vAfter1 = Prv\b;

else

    % SVD fallback

    [U,S,V] = svd(Prv);

    sv = diag(S);

    tol = max(size(Prv))*eps(max(sv));

    invS = zeros(size(sv));

    for k = 1:numel(sv)
        if sv(k)>tol
            invS(k) = 1/sv(k);
        end
    end

    vAfter1 = V*diag(invS)*U'*b;

    if norm(Prv*vAfter1-b)>1e-5
        reason = 'singular';
        return;
    end
end

% Departure impulse

burn1 = vAfter1-x0(4:6);

% Coast arrival velocity

vArrival = Pvr*x0(1:3)+Pvv*vAfter1;

% Fast analytical trajectory evaluation.
% No loop over expm() or numerical ODE integration.

Xafter1 = [x0(1:3);vAfter1];

t = linspace(0,T,N);

traj = HCW_trajectory(n,t,Xafter1);

% Arrival braking impulse

if stopAtArrival

    burn2 = xf(4:6)-vArrival;

    dvTotal = norm(burn1)+norm(burn2);

    finalState = traj(:,end);
    finalState(4:6) = ...
        finalState(4:6)+burn2;

else

    % Fly-through intermediate waypoint

    dvTotal = norm(burn1);
    finalState = traj(:,end);
end

%% TERMINAL CHECKS

positionError = norm( ...
    finalState(1:3)-xf(1:3));

if stopAtArrival
    velocityError = norm( ...
        finalState(4:6)-xf(4:6));
else
    velocityError = norm( ...
        finalState(4:6)-traj(4:6,end));
end

if positionError>1e-5 || velocityError>1e-7
    reason = 'terminal';
    return;
end

if ~isfinite(dvTotal) || ...
   ~all(isfinite(traj(:)))

    reason = 'numeric';
    return;
end

%% TRAJECTORY METRICS

ranges = vecnorm(traj(1:3,:),2,1);
speeds = vecnorm(traj(4:6,:),2,1);

sol.dvTotal = dvTotal;
sol.finalState = finalState;

sol.maxSpeed = max(speeds);
sol.minRange = min(ranges);
sol.maxRange = max(ranges);

sol.trajectory = traj;

end

%% ================================================================
% ANALYTICAL HCW STATE TRANSITION MATRIX
% ================================================================

function Phi = HCW_STM(n,t)

nt = n*t;

c = cos(nt);
s = sin(nt);

Phi = zeros(6,6);

Phi(1:3,1:3) = [ ...
    4-3*c,0,0; ...
    6*(s-nt),1,0; ...
    0,0,c];

Phi(1:3,4:6) = [ ...
    s/n,2*(1-c)/n,0; ...
    -2*(1-c)/n,(4*s-3*nt)/n,0; ...
    0,0,s/n];

Phi(4:6,1:3) = [ ...
    3*n*s,0,0; ...
    6*n*(c-1),0,0; ...
    0,0,-n*s];

Phi(4:6,4:6) = [ ...
    c,2*s,0; ...
    -2*s,4*c-3,0; ...
    0,0,c];

end

%% ================================================================
% FAST VECTORIZED HCW TRAJECTORY
% ================================================================

function traj = HCW_trajectory(n,t,X)
% t = 1xN time array [s]
% X = 6x1 post-departure-burn state
%
% Returns:
% traj = 6xN states along unforced HCW coast
%
% Uses analytical HCW equations directly.
% No repeated matrix exponentials.

nt = n*t;

c = cos(nt);
s = sin(nt);

x = X(1);
y = X(2);
z = X(3);

vx = X(4);
vy = X(5);
vz = X(6);

traj = [ ...
    (4-3*c)*x + (s/n)*vx + (2*(1-c)/n)*vy; ...
    6*(s-nt)*x + y - 2*(1-c)/n*vx + ...
        (4*s-3*nt)/n*vy; ...
    c*z + (s/n)*vz; ...
    3*n*s*x + c*vx + 2*s*vy; ...
    6*n*(c-1)*x - 2*s*vx + (4*c-3)*vy; ...
    -n*s*z + c*vz];

end

%% ================================================================
% FINAL QUINTIC APPROACH
% ================================================================

function [dv,maxSpeed,peakAccel] = ...
    finalQuinticAxis(rStart,rEnd,T,n,N)

t = linspace(0,T,N);
s = t/T;

p = 10*s.^3-15*s.^4+6*s.^5;

pd = (30*s.^2-60*s.^3+30*s.^4)/T;

pdd = (60*s-180*s.^2+120*s.^3)/T^2;

dr = rEnd-rStart;

r = rStart+dr*p;
v = dr*pd;
a = dr*pdd;

% HCW feedforward control

u = [ ...
    a(1,:)-3*n^2*r(1,:)-2*n*v(2,:); ...
    a(2,:)+2*n*v(1,:); ...
    a(3,:)+n^2*r(3,:)];

aMag = sqrt(sum(u.^2,1));

dv = trapz(t,aMag);

maxSpeed = 1.875*norm(dr)/T;
peakAccel = max(aMag);

end

%% ================================================================
% RANDOM SAMPLING
% ================================================================

function x = draw(bounds)

x = bounds(1)+ ...
    (bounds(2)-bounds(1))*rand();

end

function x0 = sampleInitialState(x0_nom,dx0)

x0 = x0_nom+(2*rand(6,1)-1).*dx0;

end

%% ================================================================
% PERCENTILE FUNCTIONS
% ================================================================

function p = percentile(x,q)

x = sort(x(:));

N = numel(x);

if N==1
    p = x;
    return;
end

r = 1+(N-1)*q/100;

lo = floor(r);
hi = ceil(r);

p = x(lo)+(r-lo)*(x(hi)-x(lo));

end

function vals = colPercentile(X,q)

vals = zeros(1,size(X,2));

for k = 1:size(X,2)
    vals(k) = percentile(X(:,k),q);
end

end

%% ================================================================
% CONFIGURATION VALIDATION
% ================================================================

function validateCfg(cfg)

assert(isstruct(cfg));
assert(cfg.N>=1 && cfg.N==floor(cfg.N));
assert(cfg.maxAttemptsFactor>=1);
assert(cfg.nTrajectorySamples>=3);

assert(cfg.mu_km3s2>0);
assert(cfg.Re_km>0);
assert(cfg.h_km>100);

assert(isequal(size(cfg.x0_client),[6 1]));
assert(isequal(size(cfg.x0_host),[6 1]));

assert(isequal(size(cfg.dx0_client),[6 1]));
assert(isequal(size(cfg.dx0_host),[6 1]));

keys = { ...
    'hp1_m','hp2_m','hp3_m','final_m', ...
    'homing_min','closing_min', ...
    'controlled_min','final_min', ...
    'hold1_min','hold2_min','hold3_min'};

for k = 1:numel(keys)

    b = cfg.(keys{k});

    assert(isnumeric(b) && numel(b)==2 && ...
        all(isfinite(b)) && ...
        b(1)>0 && b(2)>=b(1), ...
        'Invalid bounds: %s',keys{k});
end

assert(cfg.hp1_m(1)>cfg.hp2_m(2));
assert(cfg.hp2_m(1)>cfg.AZ_m);
assert(cfg.AZ_m>cfg.KOZ_m);
assert(cfg.hp3_m(2)<cfg.KOZ_m);
assert(cfg.hp3_m(1)>cfg.final_m(2));

assert(cfg.corridorHalfWidth_m>0);

assert(numel(cfg.N_homing)==2);
assert(numel(cfg.N_closing)==2);

assert(all(cfg.N_homing>=1));
assert(all(cfg.N_closing>=1));

assert(all(cfg.N_homing==floor(cfg.N_homing)));
assert(all(cfg.N_closing==floor(cfg.N_closing)));

assert(cfg.N_homing(2)>=cfg.N_homing(1));
assert(cfg.N_closing(2)>=cfg.N_closing(1));

assert(cfg.maxSpeed_homing_mps>0);
assert(cfg.maxSpeed_closing_mps>0);
assert(cfg.maxSpeed_controlled_mps>0);
assert(cfg.maxSpeed_final_mps>0);

assert(cfg.maxPeakAccel_mps2>0);

assert(iscell(cfg.methods) && ...
    ~isempty(cfg.methods) && ...
    all(strcmpi(cfg.methods,'Vbar')), ...
    'This version supports V-bar only.');

end

%% ================================================================
% PLOTS
% ================================================================

function makePlots(cfg,dv,tof,DVphase,TOFphase, ...
    labels,p95dv,p95tof,iDV,iTOF)

% Delta V histogram

f = figure('Color','w');

histogram(dv,40);
grid on;

xline(p95dv,'--','P95');

xlabel('Complete retrieval Delta V [m/s]');
ylabel('Accepted designs');
title('FABER RPO - Delta V');

saveas(f,fullfile( ...
    cfg.outputFolder,'RPO_MC_DV_histogram.png'));

% TOF histogram

f = figure('Color','w');

histogram(tof/60,40);
grid on;

xline(p95tof/60,'--','P95');

xlabel('Complete retrieval TOF [hours]');
ylabel('Accepted designs');
title('FABER RPO - TOF');

saveas(f,fullfile( ...
    cfg.outputFolder,'RPO_MC_TOF_histogram.png'));

% TOF versus Delta V

f = figure('Color','w');

scatter(tof/60,dv,12,'filled');
hold on;
grid on;

plot(tof(iDV)/60,dv(iDV),'kp', ...
    'MarkerSize',12,'MarkerFaceColor','k');

plot(tof(iTOF)/60,dv(iTOF),'ks', ...
    'MarkerSize',9,'LineWidth',2);

xline(p95tof/60,'--','P95 TOF');
yline(p95dv,'--','P95 DV');

xlabel('TOF [hours]');
ylabel('Delta V [m/s]');
title('FABER RPO - TOF vs Delta V');

saveas(f,fullfile( ...
    cfg.outputFolder,'RPO_MC_TOF_vs_DV.png'));

% Phase Delta V

f = figure('Color','w');

barh(colPercentile(DVphase,95));
grid on;

set(gca,'YTick',1:numel(labels), ...
    'YTickLabel',labels,'YDir','reverse');

xlabel('Phase P95 Delta V [m/s]');
title('Individual phase P95 Delta V');

saveas(f,fullfile( ...
    cfg.outputFolder,'RPO_MC_Phase_P95_DV.png'));

% Phase TOF

f = figure('Color','w');

barh(colPercentile(TOFphase,95));
grid on;

set(gca,'YTick',1:numel(labels), ...
    'YTickLabel',labels,'YDir','reverse');

xlabel('Phase P95 TOF [min]');
title('Individual phase P95 TOF');

saveas(f,fullfile( ...
    cfg.outputFolder,'RPO_MC_Phase_P95_TOF.png'));

% Pareto frontier

[~,idx] = sortrows([tof(:),dv(:)],[1 2]);

pareto = false(numel(dv),1);
best = inf;

for k = 1:numel(idx)

    i = idx(k);

    if dv(i)<best
        pareto(i) = true;
        best = dv(i);
    end
end

pidx = idx(pareto(idx));

f = figure('Color','w');

scatter(tof/60,dv,12,'filled');
hold on;
grid on;

plot(tof(pidx)/60,dv(pidx),'-ro', ...
    'LineWidth',2);

xlabel('TOF [hours]');
ylabel('Delta V [m/s]');
title('FABER RPO - Sampled Pareto Frontier');

saveas(f,fullfile( ...
    cfg.outputFolder,'RPO_MC_Pareto_Frontier.png'));

end
