function study = run_faber_environment_campaign(outputFolder,moduleCount)
% Reproduce two to five deliveries with one active assembly location.
if nargin<1, outputFolder = ""; end
if nargin<2, moduleCount = 2; end
validateattributes(moduleCount,{'numeric'},{'scalar','integer','>=',2,'<=',5});
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);

% Illustrative June geometry. All initial client phase angles are at most 20 deg.
budget = MA_vehicle_budget_phaseA(struct('N_RPO_MC',30, ...
    'captureWait_min',45,'handoverWait_min',10,'veicLoading_min',5, ...
    'outboundAtlanteRAAN_deg',180,'runAnnualSK',false, ...
    'includeMomentumUnloading',false));

cfg = struct('epochUTC',datetime(2038,6,21,'TimeZone','UTC'), ...
    'captureActive_min',5,'kikiHandover_min',5, ...
    'assemblyPair',[1 2],'assemblyModuleIndices',1:moduleCount, ...
    'outputFolder',outputFolder);
ids = {'Module A','Module B','Module C','Module D','Module E'};
angles = num2cell(linspace(20,10,moduleCount));
cfg.retrievals = struct('moduleId',ids(1:moduleCount), ...
    'phaseAngle_deg',angles,'cargoCount',1);
study = run_faber_environment(budget,cfg);
end
