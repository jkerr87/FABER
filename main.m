
% DV BUDGET
% cfg = struct('N_RPO_MC', 30, 'captureWait_min', 45, ...
%     'outboundDepartureEpoch_s', 600);
% out = MA_vehicle_budget_phaseA(cfg);
% print_FABER_budget(out);


% ENVIRO ANALYSIS
addpath('Environment');
% cfg = struct('epochUTC', datetime(2038,1,1,'TimeZone','UTC'), ...
%     'outputFolder', 'FABER_Environment_Study');
% cfg.retrievals = struct( ...
%     'moduleId', {'Module A', 'Module B', 'Module C'}, ...
%     'phaseAngle_deg', {20, 15, 10}, ...
%     'cargoCount', {1, 1, 1});
% 
% study = run_faber_environment(out, cfg);
study = run_faber_environment_campaign('FABER_Environment_Five', 5);
%study = run_faber_environment_campaign('FABER_Environment_Three', 3);