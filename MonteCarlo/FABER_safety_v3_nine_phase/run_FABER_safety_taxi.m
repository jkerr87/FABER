function [S,R] = run_FABER_safety_taxi(camGroups_LVLH_N,payloadMass_kg,recoveryTarget)
% RUN_FABER_SAFETY_TAXI
% Example:
% [S,R] = run_FABER_safety_taxi(verifiedNetForceVectors,actualPayloadMass,[]);
% Outputs are workspace structs. Plots/results are NOT automatically saved.
%
% Inputs:
% camGroups_LVLH_N 3xG physically verified net translation force vectors,
%   each column one independently executable simultaneous firing group [N].
%   Body-fixed thrusters require valid attitude-to-LVLH conversion BEFORE use.
% payloadMass_kg actual captured mass transported by Taxi [kg]
% recoveryTarget optional 6x1 SAFE LVLH target state (default disabled).
if nargin<3,recoveryTarget=[];end
assert(nargin>=2 && size(camGroups_LVLH_N,1)==3 && ...
    ~isempty(camGroups_LVLH_N) && ...
    all(isfinite(camGroups_LVLH_N(:))) && ...
    isscalar(payloadMass_kg) && isfinite(payloadMass_kg) && ...
    payloadMass_kg>=0, ...
    'Provide verified 3xG CAM force groups and numeric captured mass.');
if payloadMass_kg==0
    warning('Captured payload is zero: unloaded return DIAGNOSTIC only.');
end
cfg=ASV_RPO_MC_config();
cfg.N=100;cfg.seed=42;
cfg.makePlots=false;cfg.saveCSV=false;
cfg.finiteThrustValidation='none';
cfg.safetySaveBurns=true;
cfg.enablePayloadMassChange=true;
cfg.payloadMass_kg=payloadMass_kg;
cfg.enablePropulsionConstraints=true;
cfg.propulsionMode='phase';
cfg.maxBurnFraction=0.10; % match your recent MC; adjust if required
fprintf('Running paired nine-phase Taxi MC with captured mass %.3f kg.\n', ...
    payloadMass_kg);
R=FABER_RPO_9phase_MC(cfg);
% Assemble events from one Pareto-selected architecture per feasible trial.
events=struct('xPre',{},'dvNom',{},'phase',{},'vehicle',{}, ...
    'index',{},'mass_kg',{});
for i=1:cfg.N
    m=R.selectedArchitecture(i);
    if ~isfinite(m),continue;end
    rec=R.records{i,m};
    if isempty(rec),continue;end
    ec=rec.Client.safetyEvents;
    eh=rec.ATLANTE.safetyEvents;
    for k=1:numel(ec),ec(k).phase=['Client: ' ec(k).phase];end
    for k=1:numel(eh),eh(k).phase=['ATLANTE: ' eh(k).phase];end
    events=[events ec eh]; %#ok<AGROW>
end
assert(~isempty(events), ...
    'No accepted nine-phase safety events. Review RPO constraints.');
scfg=struct();
scfg.N=200;scfg.seed=4242;
scfg.vehicleMass_kg=cfg.massTaxi_kg;
scfg.camForceVectors_LVLH_N=camGroups_LVLH_N;
scfg.hazards=struct('centers_m',[0 0 0]);
scfg.recovery=struct('enabled',~isempty(recoveryTarget), ...
    'targetState',recoveryTarget);
% Add real vehicle envelope and verified recovery target/physics separately.
S=FABER_safety_MC_taxi(events,cfg.mu_km3s2, ...
    cfg.Re_km,cfg.h_km,scfg);
disp(S.summary);
fprintf('Safety events: %d | MC faults: %d | model NOT flight qualified.\n', ...
    numel(events),scfg.N);
end
