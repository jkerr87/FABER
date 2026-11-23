function R = FABER_safety_MC_taxi(events,mu,Re,h,cfg)
% Taxi failure MC operating on actual saved individual nominal burn events.
% Each event requires xPre [6x1], dvNom [3x1], phase, vehicle, index.
% Obtain events by enabling event capture in instrumented Taxi MC, NOT
% by reverse-engineering sum(DV_phase) or TOF_phase arrays.
if nargin<5 || isempty(cfg),cfg=struct();end
cfg=FABER_safety_defaults(cfg,mu,Re,h,'taxi');
if isempty(events)
    error(['Taxi MC contains no recorded burn events. ', ...
        'Run instrumented FABER_RPO_9phase_MC with cfg.safetySaveBurns=true.']);
end
R=FABER_safety_MC_core(events,cfg);
R.label='Taxi RPO impulsive-burn execution screening';
end
