function [study,budget] = MA_main(outputFolder,budget,cfg)
% Compare RPO paths, departure epochs and holding strategies.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root,fullfile(root,'Environment'),fullfile(root,'ORBMECH'), ...
    fullfile(root,'ProximityOps'),fullfile(root,'MonteCarlo','FABER_safety_v3_nine_phase'));

if nargin<1 || isempty(outputFolder), outputFolder = fullfile(pwd,'MA_Duty_Cycle'); end
if nargin<2 || isempty(budget)
    budget = MA_vehicle_budget_phaseA(struct('N_RPO_MC',120,'captureWait_min',45, ...
        'handoverWait_min',10,'veicLoading_min',5,'outboundAtlanteRAAN_deg',180));
elseif ~isfield(budget.cfg,'kikiTransferPropulsion')
    phaseDefaults = FABER_phaseA_defaults();
    budget.cfg.kikiTransferPropulsion = phaseDefaults.kikiTransferPropulsion;
end

if nargin<3, cfg = struct(); end

defaults = struct('epochUTC',datetime(2038,2,1,'TimeZone','UTC'), ...
    'useGMAT',true,'useToolbox',false,'sunlightRpoTrials',120, ...
    'sunlightRpoSearchLimit',12,'taxiDepartureSearch_min',0:15:180, ...
    'sweepHoldStrategies',[0 0;1 2;2 2;3 3],'captureActive_min',5, ...
    'assemblyModuleIndices',1:5,'assemblyTimelineDetails',true);

if ~isfield(cfg,'retrievals')
    cfg.retrievals = struct('moduleId',{'A','B','C','D','E'}, ...
        'phaseAngle_deg',num2cell(linspace(20,10,5)),'cargoCount',1);
end

for field = string(fieldnames(defaults))'
    if ~isfield(cfg,field), cfg.(field) = defaults.(field); end
end

cfg.outputFolder = outputFolder;
study = run_faber_environment(budget,cfg);
trials = study.scheduleSweep;
study.convergence = table();
for strategy = unique(trials.Strategy)'
    group = trials(trials.Strategy==strategy,:);
    stages = unique([1 1;ceil([.25;.5;.75;1]* ...
        [max(group.CandidateRank) max(group.DepartureRank)])],'rows','stable');
    previous = NaN;
    for j = 1:size(stages,1)
        included = group.CandidateRank<=stages(j,1) & group.DepartureRank<=stages(j,2);
        accepted = included & group.ScheduleAccepted;
        dv = group.TaxiDV_mps(accepted)+group.KikiFleetDV_mps(accepted);
        best = NaN; p95 = NaN;
        if ~isempty(dv), best = min(dv); p95 = prctile(dv,95); end
        row = table(strategy,group.TaxiExtraHomingHolds(1),group.AssemblyExtraHolds(1), ...
            stages(j,1),stages(j,2),nnz(included),nnz(accepted), ...
            nnz(accepted)/nnz(included),best,best-previous,p95, ...
            'VariableNames',{'Strategy','TaxiExtraHomingHolds','AssemblyExtraHolds', ...
            'CandidateCount','DepartureCount','Trials','Accepted','AcceptanceFraction', ...
            'BestRawDV_mps','ChangeInBestRawDV_mps','AcceptedP95RawDV_mps'});
        study.convergence = [study.convergence;row];
        previous = best;
    end
end

% Keep the complete Phase A rollup and reconcile just the replayed components.
% Campaign totals across separate craft are bookkeeping sums.
base = budget.budgetTable;
e = study.events;
taxi = e.Vehicle=="Taxi";
kiki = startsWith(e.Vehicle,"Kiki");

oldTaxi = base.Vehicle=="Taxi" & ismember(base.Segment, ...
    ["Outbound far","Inbound far","Outbound RPO","Inbound RPO"]);
oldKiki = base.Vehicle=="Kiki" & ismember(base.Maneuver, ...
    ["Loaded relocation","Unloaded return"]);

count = [nnz(study.retrievalSummary.Status=="Delivered");height(study.assemblySorties)];
phaseA = [sum(base.RawDV_mps(oldTaxi));sum(base.RawDV_mps(oldKiki))];

if ~any(oldTaxi), phaseA(1) = NaN; end
if ~any(oldKiki), phaseA(2) = NaN; end

scheduled = [sum(e.DeltaV_mps(taxi));sum(e.DeltaV_mps(kiki))];
meanDV = nan(2,1);

if study.ready
    meanDV(count>0) = scheduled(count>0)./count(count>0);
end

study.deltaVComparison = table(["Taxi";"Kiki fleet"],["Retrieval";"Assembly round trip"], ...
    count,phaseA,scheduled,meanDV,meanDV-phaseA,repmat(study.ready,2,1), ...
    'VariableNames',{'Vehicle','EventUnit','CompletedEvents','PhaseARawPerEvent_mps', ...
    'ScheduledCampaignRaw_mps','ScheduledRawPerEvent_mps','ChangePerEvent_mps','ReplayReady'});
study.deltaVEvents = e(e.DeltaV_mps~=0, ...
    {'Vehicle','Activity','Retrieval','StartUTC','EndUTC','DeltaV_mps'});
totals = budget.vehicleTotals;
totals.ReconciledKnownRawDV_mps = totals.KnownRawDV_mps;
totals.ReplayApplied = false(height(totals),1);

for j = 1:2
    mask = oldTaxi;
    if j==2, mask = oldKiki; end
    if ~any(mask), continue; end
    unit = unique(base.BudgetUnit(mask));
    if ~isfinite(meanDV(j)) || ~isscalar(unit), continue; end
    row = totals.Vehicle==base.Vehicle(find(mask,1)) & totals.BudgetUnit==unit;
    totals.ReconciledKnownRawDV_mps(row) = totals.KnownRawDV_mps(row) ...
        -sum(base.RawDV_mps(mask),'omitnan')+meanDV(j);
    totals.ReplayApplied(row) = true;
end

study.deltaVVehicleTotals = totals;
study.selectedReplayAcceptance = table(study.scheduleReady,study.epsReady, ...
    study.turnaroundReady,study.ready,'VariableNames', ...
    {'SchedulePass','EPSReservePass','EndTargetPass','Accepted'});
study.candidateScreen = table(cfg.sunlightRpoTrials,study.rpoPlanning.Candidates, ...
    numel(unique(trials.CandidateId(trials.CandidateId>0))),study.rpoPlanning.Diagnostic, ...
    'VariableNames',{'RequestedMCTrials','FeasibleAlternates','ScheduledAlternates','Diagnostic'});

% output
writetable(study.candidateScreen,fullfile(outputFolder,'MA_Candidate_Screen.csv'));
writetable(trials,fullfile(outputFolder,'MA_Schedule_Trials.csv'));
writetable(study.convergence,fullfile(outputFolder,'MA_Search_Convergence.csv'));
writetable(study.selectedReplayAcceptance,fullfile(outputFolder,'MA_Selected_Replay_Acceptance.csv'));
writetable(base,fullfile(outputFolder,'MA_PhaseA_Budget.csv'));
writetable(totals,fullfile(outputFolder,'MA_DeltaV_Vehicle_Totals.csv'));
writetable(study.deltaVComparison,fullfile(outputFolder,'MA_Scheduled_DeltaV_Comparison.csv'));
writetable(study.deltaVEvents,fullfile(outputFolder,'MA_Scheduled_DeltaV_Events.csv'));
disp(study.convergence);
disp(study.deltaVComparison);
end
