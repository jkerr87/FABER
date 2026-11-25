function analysis = run_FABER_sunlight_eclipse_study(budget,outputFolder,moduleCount,makePlots)
% Compare active-RPO sunlight gating and eclipse-tolerant far-field coast.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root,fullfile(root,'Environment'),fullfile(root,'ORBMECH'));
if nargin<1 || isempty(budget)
    budget = MA_vehicle_budget_phaseA(struct('N_RPO_MC',30,'captureWait_min',45, ...
        'handoverWait_min',10,'veicLoading_min',5,'outboundAtlanteRAAN_deg',180, ...
        'runAnnualSK',false,'includeMomentumUnloading',false));
end
if nargin<2 || isempty(outputFolder), outputFolder = fullfile(pwd,'SunlightEclipseResults'); end
if nargin<3, moduleCount = 5; end
if nargin<4, makePlots = false; end
validateattributes(moduleCount,{'numeric'},{'scalar','integer','>=',1,'<=',5});
if ~isfolder(outputFolder), mkdir(outputFolder); end

epochs = [datetime(2038,2,1,'TimeZone','UTC'),datetime(2038,6,21,'TimeZone','UTC')];
seasons = ["Winter","Summer"];
extra = 0:2;
ids = "Module "+string(1:moduleCount);
angles = 20;
if moduleCount>1, angles = linspace(20,10,moduleCount); end
cfg = struct('useToolbox',false,'captureActive_min',min(5,budget.cfg.captureWait_min), ...
    'assemblyModuleIndices',1:moduleCount,'makePlots',makePlots,'captureRequiresSun',true);
cfg.retrievals = struct('moduleId',cellstr(ids),'phaseAngle_deg',num2cell(angles),'cargoCount',1);
summary = table();
farField = table();
runs = cell(numel(epochs),numel(extra));
Re = 1000*astroConstants(23);
for season = 1:numel(epochs)
    cfg.epochUTC = epochs(season);
    for j = 1:numel(extra)
        cfg.rpoExtraHomingHolds = extra(j);
        caseName = seasons(season)+"_Holds"+extra(j);
        cfg.outputFolder = fullfile(outputFolder,caseName);
        study = run_faber_environment(budget,cfg);
        runs{season,j} = study;
        far = auditFarField(study,caseName,Re);
        farField = [farField;far]; %#ok<AGROW>
        e = study.events;
        staging = contains(e.Activity,"Homing 1 stage") & (startsWith(e.Activity,"Client") | startsWith(e.Activity,"Atlante"));
        waits = startsWith(e.Activity,"Eclipse wait before ");
        violation = sum(study.lightingSummary.SunRequiredViolation_min);
        summary = [summary;table(caseName,cfg.epochUTC,extra(j), ...
            nnz(staging & ~endsWith(e.Activity,["2 of 2","3 of 3","4 of 4"])), ...
            nnz(waits & e.Vehicle=="Taxi"),nnz(waits & contains(e.Activity,"Homing 1 stage")), ...
            nnz(study.retrievalSummary.Status=="Delivered"),numel(study.completedAssemblyModules), ...
            sum(far.CoastEclipse_min),nnz(far.FullSunThroughout), ...
            nnz(far.DepartSunFraction<cfgSunFraction(study))+nnz(far.ArriveSunFraction<cfgSunFraction(study)), ...
            violation,sum(e.DeltaV_mps(e.Vehicle=="Taxi"),'omitnan'),study.epsReady, ...
            'VariableNames',{'Case','EpochUTC','RequestedExtraHoldsPerHoming', ...
            'ExtraGatesUsed','TaxiScheduledWaits','StagedHomingWaits','Deliveries','Installations', ...
            'FarFieldEclipse_min','LegsPassingContinuousSunAudit','BurnEpochsBelowFullSun', ...
            'ActiveSunViolation_min','TaxiScheduledDV_mps','EPSPass'})]; %#ok<AGROW>
    end
end
analysis = struct('summary',summary,'farField',farField,'runs',{runs});
writetable(summary,fullfile(outputFolder,'Sunlight_Eclipse_Comparison.csv'));
writetable(farField,fullfile(outputFolder,'Far_Field_Lighting.csv'));
disp(summary);
end

function far = auditFarField(study,caseName,Re)
e = study.events(study.events.Vehicle=="Taxi" & ...
    ismember(study.events.Activity,["Outbound Lambert","Return Lambert"]),:);
sun = FABER_environment_solar(study.options.epochUTC,study.options.maxDays,false);
fraction = nan(height(e),2);
for i = 1:height(e)
    direction = 1+(e.Activity(i)=="Return Lambert");
    arc = study.retrievalArcs{e.Retrieval(i),direction};
    t = [e.Start_s(i),e.End_s(i)];
    x = [arc.departurePreBurnECI,arc.arrivalPostBurnECI];
    fraction(i,:) = FABER_environment_light(sun,t,x,study.options.aspectAxisLVLH,Re);
end
full = e.Eclipse_min<=1e-8 & all(fraction>=cfgSunFraction(study),2);
far = table(repmat(caseName,height(e),1),e.Retrieval,e.Activity,e.StartUTC,e.EndUTC, ...
    e.Sunlight_min,e.Eclipse_min,fraction(:,1),fraction(:,2),full, ...
    'VariableNames',{'Case','Retrieval','Transfer','DepartureUTC','ArrivalUTC', ...
    'CoastSunlight_min','CoastEclipse_min','DepartSunFraction','ArriveSunFraction','FullSunThroughout'});
end

function threshold = cfgSunFraction(study)
threshold = study.options.sunRequiredFraction;
end
