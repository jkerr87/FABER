function [screen,planning] = FABER_environment_rpo_candidates(budget,cfg,mu,Re)
% Replan environment-only RPO candidates around stationary V-bar entry gates.
% Powered legs remain subject to the normal nine-phase feasibility screens.
screen = budget;
planning = struct('Candidates',0,'Diagnostic',"");
phase = budget.cfg;
phase.x0_client = [0;-cfg.rpoEntryHoldDistance_m;0;0;0;0];
phase.x0_atlante = phase.x0_client;
phase.N_RPO_MC = cfg.sunlightRpoTrials;
try
    outbound = FABER_outbound_lambert(mu,Re,phase);
catch err
    if strcmp(err.identifier,'FABER:OutboundLambertNoSolution')
        planning.Diagnostic = string(err.message);
        return;
    end
    rethrow(err);
end
mc = FABER_phaseA_rpo_config(phase,mu,Re,outbound);
% Target the real return Lambert arc after the capture and wait schedule
% has been determined. Screen the loaded approach from the stationary gate.
mc = rmfield(mc,'returnHandoff');
mc.arrivalStatesHost = repmat(phase.x0_atlante,1,mc.N);
mc.methodsH2 = {'GeneralHCW'};
mc.homing1_min = [30 cfg.sunlightRpoMaxTransfer_min];
mc.homing2_min = [20 cfg.sunlightRpoMaxTransfer_min];
mc.closing1_min = [10 cfg.sunlightRpoMaxTransfer_min];
mc.closing2_min = [5 cfg.sunlightRpoMaxClosing_min];
mc.final_min = [5 cfg.sunlightRpoMaxFinal_min];
% Preserve applicable mission constraints from the baseline run.
source = budget.RPO_MC.cfg;
names = string(fieldnames(source))';
for f = names
    if ~isempty(regexp(f,'^(max|AZ|KOZ|clientGeometry|hostGeometry|terminal|enablePropulsion|checkBrakingDistance|propulsion|phasePropulsion|closingSpeedMode)', ...
            'once'))
        mc.(f) = source.(f);
    end
end
if ~isempty(fieldnames(cfg.clientHoldGeometry)), mc.clientGeometry=cfg.clientHoldGeometry; end
if ~isempty(fieldnames(cfg.atlanteHoldGeometry)), mc.hostGeometry=cfg.atlanteHoldGeometry; end
rngState = rng;
restore = onCleanup(@() rng(rngState)); %#ok<NASGU>
R = [];
evalc('R = FABER_RPO_9phase_MC(mc);');
valid = find(R.valid(:,1));
terminal = arrayfun(@(i) sum(R.records{i,1}.Client.tofPhase(7:9)),valid);
[~,order] = sortrows([terminal(:) R.CycleTOF_min(valid,1)],[1 2]);
planning.CandidateIds = valid(order);
planning.Candidates = numel(valid);
if isempty(valid)
    planning.Diagnostic = "No shorter paired candidate passed the nine-phase feasibility checks.";
end
screen.cfg = phase;
screen.outboundLambert = outbound;
screen.RPO_MC = R;
end
