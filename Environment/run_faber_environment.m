function out = run_faber_environment(budget,cfg)
% UTC sunlight scheduling and PM4 EPS replay for selected Phase A operations.
% A Phase A result can be supplied directly. Omit it for a reference retrieval.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'ORBMECH'),fullfile(root,'ProximityOps'), ...
    fullfile(root,'Environment'), ...
    fullfile(root,'MonteCarlo','FABER_safety_v3_nine_phase'),'-begin');
if nargin<1 || isempty(budget)
    budget = MA_vehicle_budget_phaseA(struct('N_RPO_MC',30,'captureWait_min',45));
end
if nargin<2, cfg=struct(); end
cfg = inputs(cfg,budget);
mu = 1e9*astroConstants(13);
Re = 1000*astroConstants(23);
solar = FABER_environment_solar(cfg.epochUTC,cfg.maxDays,cfg.useToolbox);

% Preserve the unretimed Phase A reference alongside the gated schedule.
reference = FABER_environment_schedule(budget,cfg,solar,mu,Re,false);
plan = FABER_environment_schedule(budget,cfg,solar,mu,Re,true);
planning = struct('Candidates',0,'Diagnostic',"",'SelectedTrial',NaN,'UsedAlternate',false,'DepartureDelay_min',0,'DepartureTrials',0);
departureBudget = budget;
scheduleCfg = cfg;
if cfg.replanSunlightRPO && any(plan.retrievalSummary.Status~="Delivered")
    [screen,diagnostic] = FABER_environment_rpo_candidates(budget,cfg,mu/1e9,Re/1000);
    planning.Candidates = diagnostic.Candidates;
    planning.Diagnostic = diagnostic.Diagnostic;
    best = rpoProgress(plan);
    if diagnostic.Candidates>0
        for i = diagnostic.CandidateIds(1:min(end,cfg.sunlightRpoSearchLimit))'
            screen.RPO_MC.SelectedCycle = i;
            candidate = FABER_environment_schedule(screen,cfg,solar,mu,Re,true);
            score = rpoProgress(candidate);
            if score>best
                plan = candidate;
                departureBudget = screen;
                best = score;
                planning.SelectedTrial = i;
                planning.UsedAlternate = true;
            end
            if all(candidate.retrievalSummary.Status=="Delivered"), break; end
        end
    end
end
if any(plan.retrievalSummary.Status~="Delivered") && ...
        any(plan.issues.Vehicle=="Taxi" & (contains(plan.issues.Reason,"sunlight", ...
        'IgnoreCase',true) | contains(plan.issues.Reason,"eclipse",'IgnoreCase',true)))
    best = rpoProgress(plan);
    for delay = cfg.taxiDepartureSearch_min(:)'
        trial = cfg;
        trial.taxiInitialWait_s = delay*60;
        candidate = FABER_environment_schedule(departureBudget,trial,solar,mu,Re,true);
        planning.DepartureTrials = planning.DepartureTrials+1;
        score = rpoProgress(candidate);
        if score>best
            plan = candidate;
            best = score;
            planning.DepartureDelay_min = delay;
            scheduleCfg = trial;
        end
        if all(candidate.retrievalSummary.Status=="Delivered"), break; end
    end
end
names = ["Taxi","Atlante","Kiki 1","Kiki 2","Rango 1","Rango 2"];
if cfg.chargeAsNeeded && ~isempty(plan.assemblySorties)
    for attempt = 1:3*numel(cfg.assemblyModuleIndices)
        [powerEdges,powerMid] = sampleTimes(plan,cfg.sampleStep_s);
        powerStates = sampledTracks(plan,powerMid,names);
        [powerSun,powerAspect] = lighting(solar,powerMid,powerStates,cfg,Re);
        [dispatch,powerSummary] = runEPS(cfg,plan,powerEdges,powerSun,powerAspect, ...
            names,budget.cfg.missionYears);
        chargePlan = nextSortieCharge(plan.assemblySorties,dispatch,powerSummary,scheduleCfg);
        if isequal(chargePlan,scheduleCfg.assemblyCharge_min), break; end
        scheduleCfg.assemblyCharge_min = chargePlan;
        plan = FABER_environment_schedule(departureBudget,scheduleCfg,solar,mu,Re,true);
    end
end
cfg.assemblyCharge_min = scheduleCfg.assemblyCharge_min;
[edges,mid] = sampleTimes(plan,cfg.sampleStep_s);
[states,uses,labels] = sampledTracks(plan,mid,names);
[frac,aspect,beta] = lighting(solar,mid,states,cfg,Re);
[toolboxFractions,toolboxIntervals,toolboxDiagnostic] = toolboxAudit(plan,mid,names,cfg);
combined = frac;
hasToolbox = isfinite(toolboxFractions);
combined(hasToolbox) = min(frac(hasToolbox),toolboxFractions(hasToolbox));
scheduledViolations = uses & (~isfinite(combined) | combined<cfg.sunRequiredFraction);
lightTrace = traceTable(cfg.epochUTC,mid,names,labels,states,aspect,beta, ...
    frac,toolboxFractions,combined,uses,scheduledViolations);
[referenceEdges,refMid] = sampleTimes(reference,cfg.sampleStep_s);
[refState,refUses,refLabels] = sampledTracks(reference,refMid,names);
[refFraction,refAspect,refBeta] = lighting(solar,refMid,refState,cfg,Re);
referenceTrace = traceTable(cfg.epochUTC,refMid,names,refLabels,refState, ...
    refAspect,refBeta,refFraction,nan(size(refFraction)),refFraction,refUses, ...
    refUses & (~isfinite(refFraction) | refFraction<cfg.sunRequiredFraction));
lightingSummary = table(names',sum(diff(edges).*(combined>=1-1e-9),1)'/60, ...
    sum(diff(edges).*(combined>0 & combined<1-1e-9),1)'/60, ...
    sum(diff(edges).*(combined==0 & isfinite(frac)),1)'/60, ...
    sum(diff(edges).*~isfinite(frac),1)'/60, ...
    sum(diff(edges).*scheduledViolations,1)'/60, ...
    'VariableNames',{'Vehicle','Sunlight_min','Penumbra_min','Umbra_min', ...
    'UnknownPosition_min','SunRequiredViolation_min'});

[eps,epsSummary,lifetime,cycles,charge,batteryHistogram] = runEPS(cfg,plan,edges, ...
    combined,aspect,names,budget.cfg.missionYears);
out = struct('options',cfg,'lightingSource',solar.source,'toolboxDiagnostic',toolboxDiagnostic, ...
    'events',eventLighting(plan.events,cfg,edges,names,combined,aspect), ...
    'links',withUTC(plan.links,cfg.epochUTC), ...
    'holds',withUTC(plan.holds,cfg.epochUTC),'openItems',plan.issues, ...
    'lightingTrace',lightTrace,'lightingSummary',lightingSummary,'eclipseIntervals', ...
    computedIntervals(edges,combined,frac,names,cfg.epochUTC), ...
        'toolboxEclipseIntervals',toolboxIntervals, ...
    'referenceEvents',eventLighting(reference.events,cfg,referenceEdges,names,refFraction,refAspect), ...
    'referenceLightingTrace',referenceTrace, ...
    'eps',eps,'epsSummary',epsSummary,'dutyCycle',dutyCycle(eps,names,combined), ...
    'chargeSummary',charge,'batteryLifetime',lifetime, ...
    'batteryCycles',cycles,'batteryHistogram',batteryHistogram, ...
    'returnLambert',plan.returnLambert,'holdRecoveryDV_mps',plan.recoveryDV_mps, ...
    'retrievalSummary',retrievalUTC(plan.retrievalSummary,cfg.epochUTC), ...
    'referenceRetrievalSummary',retrievalUTC(reference.retrievalSummary,cfg.epochUTC), ...
    'retrievalArcs',{plan.retrievalArcs}, ...
        'completedAssemblyModules',plan.completedModules,'assemblySorties',plan.assemblySorties, ...
    'rpoPlanning',planning, ...
    'assemblyReferencePaths',{assemblyPaths(budget,cfg.assemblyModuleIndices)});
out.openItems.UTC = cfg.epochUTC+seconds(out.openItems.Epoch_s);
out.scheduleReady = isempty(plan.issues) && ~any(scheduledViolations,'all');
out.epsReady = all(epsSummary.Pass);
out.turnaroundReady = all(epsSummary.EndTargetPass);
out.ready = out.scheduleReady && out.epsReady;
out.figures = gobjects(0);
out.figureTabs = gobjects(0);
if cfg.makePlots
    [out.figures,out.figureTabs] = FABER_environment_plots(out,plan,reference,solar,names,Re);
end
if strlength(cfg.outputFolder)>0
    saveOutputs(out);
end
disp(out.retrievalSummary);
disp(out.lightingSummary);
disp(out.epsSummary);
if ~isempty(out.openItems), disp(out.openItems); end
end

function minutes = nextSortieCharge(sorties,run,summary,cfg)
minutes = cfg.assemblyCharge_min;
floorEnergy = min(run.E,[],1)-summary.MinReserveMargin_kWh';
for i = 1:height(sorties)
    row = sorties(i,:);
    [~,start] = min(abs(run.t-row.LoadStart_s/60));
    [~,finish] = min(abs(run.t-row.BerthEnd_s/60));
    v = [row.Pair+2,row.Pair+4];
    samples = start:finish-1;
    if isempty(samples), continue; end

    % Recover the energy needed for unmet load as well as the reserve.
    % Do not infer the deficit from SoC alone because dispatch clips it at Emin.
    unmet = [zeros(1,2);cumsum(run.unserved(samples,v).*run.dt_h(samples),1)]/0.90;
    required = max(floorEnergy(v)-min(run.E(start:finish,v)-unmet,[],1),0);
    headroom = max(run.EOL(v)-run.E(start,v),0);
    required = min(required,headroom);
    if ~cfg.retainRangoForRelocations, required(2) = 0; end
    if all(required<1e-5), continue; end

    rate = min(max(run.assemblyChargeNet_kW,0),run.chargeLimit(start,v))*0.85;
    chargeable = required>=1e-5 & rate>0;
    % A Rango charge can also replenish Kiki via the Atlante feed.
    if chargeable(2), chargeable(1) = false; end
    if any(chargeable)
        extra = ceil(120*required(chargeable)./rate(chargeable))/2;
        minutes(row.Queue,chargeable) = minutes(row.Queue,chargeable)+extra;
        return;
    end
end
end

function paths = assemblyPaths(budget,indices)
paths = cell(0,1);
for i = indices(:)'
    x = budget.assemblyRelocations.bestRecords{i};
    if isempty(x), continue; end
    for direction = ["outbound","returnTrip"]
        leg = x.(direction);
        paths{end+1,1} = [leg.exitPath(1:3,:), ...
            leg.transferPath(1:3,2:end),leg.finalPath(1:3,2:end)]; %#ok<AGROW>
    end
end
end

function score = rpoProgress(plan)
% Prefer a completed retrieval. Otherwise retain the furthest safely
% scheduled phase, independent of the number of standby intervals.
delivered = nnz(plan.retrievalSummary.Status=="Delivered");
labels = plan.events.Activity(plan.events.Vehicle=="Taxi" & plan.events.Retrieval==delivered+1);
stages = ["Outbound Lambert","Client Homing 1","Client Homing 2", ...
    "Client Closing 1","Client closing and final approach","Capture operations", ...
    "Return Lambert","Atlante Homing 1","Atlante Homing 2","Atlante Closing 1", ...
    "Atlante closing and final approach","Taxi to Atlante handover"];
completed = find(ismember(stages,labels));
score = max([0 completed],[],'all') + 20*delivered;
end

function cfg = inputs(given,budget)
cfg = struct('epochUTC',datetime(2038,1,1,'TimeZone','UTC'), ...
    'useToolbox',true,'sampleStep_s',30,'searchStep_s',30,'sunGuard_s',60, ...
    'maxWait_s',6*3600,'rpoExtraHomingHolds',0,'assemblyDepartureSearchStep_s',900, ...
    'taxiDepartureSearch_min',[15 30 45 60],'taxiInitialWait_s',0,'maxDays',21,'aspectAxisLVLH',[0;1;0], ...
    'sunRequiredFraction',1-1e-8,'rpoHoldKeepOut_m',30,'assemblyHoldKeepOut_m',30, ...
    'replanSunlightRPO',true,'rpoEntryHoldDistance_m',2000,'sunlightRpoTrials',30, ...
    'sunlightRpoSearchLimit',10,'sunlightRpoMaxTransfer_min',52, ...
    'sunlightRpoMaxClosing_min',22,'sunlightRpoMaxFinal_min',15, ...
    'holdErrorLVLH',[0.5;0.5;0.5;0.002;0.002;0.002], ...
    'maxHoldDrift_m',100,'holdContinuationPosition_m',0.01, ...
    'holdContinuationVelocity_mps',0.001,'recoveryTime_s',600,'recoveryIsp_s',293, ...
     'maxRecoveryDV_mps',1,'maxRecoverySpeed_mps',0.3, ...
    'clientHoldGeometry',struct(),'atlanteHoldGeometry',struct(),'assemblyHoldGeometry',struct(), ...
    'assemblyModuleIndices',1,'preloadedModules',0,'modulesPerRetrieval',1, ...
    'retrievals',struct([]),'taxiDeparture_min',12,'taxiDryMass_kg',NaN, ...
     'assemblyPair',[1 2],'captureRequiresSun',true,'captureActive_min',NaN, ...
    'installationRequiresSun',false, ...
    'retainRangoForRelocations',true,'rangoTransfer_min',16.67,'dockedPairSupport',true, ...
    'taxiHandover_min',budget.cfg.handoverWait_min,'kikiLoad_min',budget.cfg.veicLoading_min, ...
    'kikiHandover_min',5,'rangoAssembly_min',32.17,'assemblyBerth_min',5, ...
    'rangoChargeConnect_min',5,'chargeAsNeeded',true,'assemblyCharge_min',[], ...
    'kikiCharge_min',budget.cfg.veicCharge_min, ...
    'rangoCharge_min',98.427554978,'kikiPort_kW',3,'rangoPort_kW',1.5, ...
     'standbyAfter_min',0,'dutyCycleZoom_min',600,'cluster',"ISS",'initialSoC',ones(1,6), ...
    'includeRefueling',false,'taxiRefuel_min',NaN,'atlanteTanker_min',NaN, ...
    'kikiRefuel_min',256.15,'postFluidTopup_min',105,'safeReserve_min',15,'chargeTargetSoC',1, ...
    'pvUsesAspect',false,'journeyFraction',0,'makePlots',true,'outputFolder',"");
% The PM4 reference charge window is used when Phase A has no Kiki charge input.
if isnan(cfg.kikiCharge_min), cfg.kikiCharge_min = 104.417845902; end
if isnan(cfg.kikiLoad_min), cfg.kikiLoad_min = 5; end
for f = string(fieldnames(given))'
    assert(isfield(cfg,f),'Unknown environment setting %s.',f);
    cfg.(f) = given.(f);
end
cfg.epochUTC.TimeZone = 'UTC';
assert(isscalar(cfg.epochUTC) && ~isnat(cfg.epochUTC) && year(cfg.epochUTC)>=2038,'Supply an epoch in 2038 or later.');
cfg.cluster = upper(string(cfg.cluster));
assert(isscalar(cfg.cluster) && any(cfg.cluster==["ISS","SSO"]),'Cluster must be ISS or SSO.');
ix = 1+(cfg.cluster=="SSO");
taxi = [529 539]; tanker = [667 405];
if isnan(cfg.taxiRefuel_min), cfg.taxiRefuel_min=taxi(ix); end
if isnan(cfg.atlanteTanker_min), cfg.atlanteTanker_min=tanker(ix); end
validateattributes(cfg.taxiRefuel_min,{'numeric'},{'scalar','finite','>=',10.84});
validateattributes(cfg.chargeTargetSoC,{'numeric'},{'scalar','finite','>=',0.65,'<=',1});
validateattributes(cfg.aspectAxisLVLH,{'numeric'},{'vector','numel',3,'finite'});
assert(norm(cfg.aspectAxisLVLH)>0,'Aspect axis cannot be zero.');
cfg.aspectAxisLVLH = cfg.aspectAxisLVLH(:)/norm(cfg.aspectAxisLVLH);
validateattributes(cfg.rpoEntryHoldDistance_m,{'numeric'},{'scalar','finite','positive'});
for f = ["sunlightRpoTrials","sunlightRpoSearchLimit"]
    validateattributes(cfg.(f),{'numeric'},{'scalar','integer','finite','positive'});
end
validateattributes(cfg.sunlightRpoMaxTransfer_min,{'numeric'},{'scalar','finite','>=',30});
validateattributes(cfg.sunlightRpoMaxClosing_min,{'numeric'},{'scalar','finite','>=',5});
validateattributes(cfg.sunlightRpoMaxFinal_min,{'numeric'},{'scalar','finite','>=',5});
validateattributes(cfg.captureActive_min,{'numeric'},{'scalar','real'});
if ~isnan(cfg.captureActive_min)
    validateattributes(cfg.captureActive_min,{'numeric'},{'scalar','finite','nonnegative', ...
        '<=',budget.cfg.captureWait_min});
end
validateattributes(cfg.initialSoC,{'numeric'},{'vector','numel',6,'finite','>=',0.65,'<=',1});
cfg.initialSoC = cfg.initialSoC(:)';
validateattributes(cfg.holdErrorLVLH,{'numeric'},{'vector','numel',6,'finite','nonnegative'});
validateattributes(cfg.sunRequiredFraction,{'numeric'},{'scalar','positive','<=',1});
validateattributes(cfg.journeyFraction,{'numeric'},{'scalar','nonnegative','<=',1});
for f = ["sampleStep_s","searchStep_s","maxDays","recoveryTime_s","recoveryIsp_s","dutyCycleZoom_min"]
    validateattributes(cfg.(f),{'numeric'},{'scalar','real','finite','positive'});
end
for f = ["maxWait_s","sunGuard_s","standbyAfter_min","maxHoldDrift_m","rpoHoldKeepOut_m", ...
        "assemblyHoldKeepOut_m","holdContinuationPosition_m","holdContinuationVelocity_mps", ...
        "maxRecoveryDV_mps","maxRecoverySpeed_mps","kikiPort_kW","rangoPort_kW", ...
        "atlanteTanker_min","kikiRefuel_min","safeReserve_min","postFluidTopup_min"]
    validateattributes(cfg.(f),{'numeric'},{'scalar','real','finite','nonnegative'});
end
for f = ["preloadedModules","modulesPerRetrieval"]
    validateattributes(cfg.(f),{'numeric'},{'scalar','integer','finite','nonnegative'});
end
validateattributes(cfg.taxiInitialWait_s,{'numeric'}, ...
    {'scalar','real','finite','nonnegative','<',cfg.maxDays*86400});
validateattributes(cfg.rpoExtraHomingHolds,{'numeric'}, ...
    {'scalar','integer','nonnegative','<=',3});
validateattributes(cfg.assemblyDepartureSearchStep_s,{'numeric'}, ...
    {'scalar','real','finite','positive'});
if ~isempty(cfg.taxiDepartureSearch_min)
    validateattributes(cfg.taxiDepartureSearch_min,{'numeric'}, ...
        {'vector','real','finite','positive','<',cfg.maxDays*1440});
end
validateattributes(cfg.taxiDeparture_min,{'numeric'},{'scalar','finite','positive'});
validateattributes(cfg.taxiDryMass_kg,{'numeric'},{'scalar','real'});
if ~isnan(cfg.taxiDryMass_kg)
    validateattributes(cfg.taxiDryMass_kg,{'numeric'},{'positive','<=',budget.cfg.massTaxi_kg});
end
assert(isstruct(cfg.retrievals),'Supply retrievals as a struct array.');
if ~isempty(cfg.retrievals)
    assert(isvector(cfg.retrievals) && all(isfield(cfg.retrievals, ...
        {'moduleId','phaseAngle_deg','cargoCount'})), ...
        'Each retrieval needs moduleId, phaseAngle_deg and cargoCount.');
    ids = strings(numel(cfg.retrievals),1);
    for trip = 1:numel(cfg.retrievals)
        r = cfg.retrievals(trip);
        moduleId = string(r.moduleId);
        assert(isscalar(moduleId) && ~ismissing(moduleId) && strlength(strtrim(moduleId))>0, ...
            'Supply a nonempty moduleId for each retrieval.');
        ids(trip) = moduleId;
        validateattributes(r.phaseAngle_deg,{'numeric'},{'scalar','real','finite'});
        validateattributes(r.cargoCount,{'numeric'},{'scalar','integer','finite','positive'});
    end
    assert(numel(unique(ids))==numel(ids) && ...
        numel(unique(mod([cfg.retrievals.phaseAngle_deg],360)))==numel(ids), ...
        'Each retrieval needs a distinct module ID and epoch-zero phase angle.');
end
validateattributes(cfg.assemblyBerth_min,{'numeric'},{'scalar','finite','nonnegative'});
validateattributes(cfg.rangoChargeConnect_min,{'numeric'},{'scalar','finite','nonnegative'});
cfg.assemblyCharge_min = zeros(numel(cfg.assemblyModuleIndices),2);
assert(all(ismember(cfg.assemblyPair,[1 2])) && ~isempty(cfg.assemblyPair),'Use assembly pair 1 or 2.');
ids = cfg.assemblyModuleIndices;
assert(all(isfinite(ids) & ids>=1 & ids==floor(ids) & ...
    ids<=numel(budget.assemblyRelocations.bestRecords)) && numel(unique(ids))==numel(ids), ...
    'Assembly indices must refer to unique Phase A module locations.');
end

function [edges,mid] = sampleTimes(plan,step)
edges = unique([0;(0:step:plan.stop_s)';plan.events.Start_s;plan.events.End_s; ...
    plan.links.Start_s;plan.links.End_s;plan.stop_s]);
mid = (edges(1:end-1)+edges(2:end))/2;
end

function trace = traceTable(epoch,time,names,labels,states,aspect,beta, ...
    geometry,toolbox,fraction,required,violations)
trace = table();
for v=1:numel(names)
    x = states(:,:,v)';
    trace = [trace; table(epoch+seconds(time),repmat(names(v),numel(time),1), ...
        labels(:,v),x(:,1),x(:,2),x(:,3),aspect(:,v),beta(:,v), ...
        geometry(:,v),toolbox(:,v),fraction(:,v),required(:,v),violations(:,v), ...
        'VariableNames',{'UTC','Vehicle','Activity','ECI_X_m','ECI_Y_m','ECI_Z_m', ...
        'SunAspect_deg','Beta_deg','GeometrySunFraction','ToolboxSunFraction','SunFraction', ...
        'SunRequired','SunViolation'})]; %#ok<AGROW>
end
end

function [states,sunRequired,labels] = sampledTracks(plan,time,names)
N = numel(time);
states = repmat(plan.atlante(time),1,1,numel(names));
sunRequired = false(N,numel(names));
labels = repmat("Standby",N,numel(names));
states(:,:,1) = plan.taxiPark(time);
if plan.taxiDocked
    mask = time(:)>=plan.taxiEnd_s;
    states(:,mask,1) = plan.taxiBerth(time(mask));
else
    mask = time(:)>plan.taxiEnd_s;
    states(:,mask,1) = NaN;
    labels(mask,1) = "Post-stop state open";
end
for i=1:numel(plan.paths)
    p = plan.paths{i};
    v = find(names==p.vehicle,1);
    mask = time(:)>=p.start & time(:)<p.stop;
    if any(mask)
        states(:,mask,v) = p.state(time(mask));
        sunRequired(mask,v) = p.sunRequired;
        labels(mask,v) = p.label;
    end
end
end

function [fraction,aspect,beta] = lighting(solar,t,states,cfg,Re)
fraction = zeros(numel(t),size(states,3));
aspect = nan(size(fraction));
beta = aspect;
for v=1:size(states,3)
    x=states(:,:,v);
    valid = all(isfinite(x),1);
    fraction(~valid,v)=NaN;
    [f,a,b] = FABER_environment_light(solar,t(valid),x(:,valid),cfg.aspectAxisLVLH,Re);
    fraction(valid,v)=f;
    aspect(valid,v)=a;
    beta(valid,v)=b;
end
end

function [values,intervals,diagnostic] = toolboxAudit(plan,time,names,cfg)
values = nan(numel(time),numel(names));
intervals = table();
diagnostic = "Aerospace Toolbox audit disabled.";
if ~cfg.useToolbox || isempty(time), return; end
diagnostic = "";
for v=1:numel(names)
    try
        horizon = plan.stop_s;
        open = plan.events.Start_s(plan.events.Vehicle==names(v) & plan.events.Activity=="Post-stop state open");
        if ~isempty(open), horizon=min(open); end
        if v==1 && ~plan.taxiDocked, horizon = min(horizon,plan.taxiEnd_s); end
        % Supply position samples including exact maneuver boundaries.
        sample = unique([0:cfg.sampleStep_s:horizon,horizon, ...
            plan.events.Start_s(plan.events.Start_s<horizon)', ...
            plan.events.End_s(plan.events.End_s<horizon)'])';
        if numel(sample)<2, continue; end
        s = sampledTracks(plan,sample,names);
        if (~isempty(open) || (v==1 && ~plan.taxiDocked))
            % Use the left-hand state at the last modeled endpoint.
            p = find(cellfun(@(x) x.vehicle==names(v) && x.stop==horizon,plan.paths),1,'last');
            if ~isempty(p), s(:,end,v) = plan.paths{p}.state(horizon); end
        end
        ta = cfg.epochUTC+seconds(sample);
        sc = satelliteScenario(ta(1),ta(end),cfg.sampleStep_s);
        pos = timetable(ta,s(1:3,:,v)','VariableNames',{'Position'});
        sat = satellite(sc,pos,CoordinateFrame="inertial",Name=names(v));
        ecl = eclipse(sat,IncludeLunarEclipse=false);
        [status,times] = eclipseStatus(ecl);
        values(:,v) = interp1(seconds(times(:)-cfg.epochUTC),status(:),time(:),'linear',NaN);
        intervals = [intervals;eclipseIntervals(ecl)]; %#ok<AGROW>
    catch err
        if strlength(diagnostic)==0, diagnostic=string(err.message); end
        % Missing Toolbox or ephemeris is common to all assets.
        if v==1, break; end
    end
end
if strlength(diagnostic)==0, diagnostic="Aerospace Toolbox eclipse audit completed."; end
end

function result = computedIntervals(edges,fraction,known,names,epoch)
result = table();
for v=1:numel(names)
    shadow = isfinite(known(:,v)) & fraction(:,v)<1-1e-9;
    change = diff([false;shadow;false]);
    a = find(change==1); b = find(change==-1);
    for k=1:numel(a)
        result=[result;table(names(v),epoch+seconds(edges(a(k))),epoch+seconds(edges(b(k))), ...
            (edges(b(k))-edges(a(k)))/60,min(fraction(a(k):b(k)-1,v)), ...
            'VariableNames',{'Vehicle','StartUTC','EndUTC','Duration_min','MinimumSunFraction'})]; %#ok<AGROW>
    end
end
end

function [run,summary,life,cycles,charge,histogram] = runEPS(cfg,plan,edges,fraction,aspect,names,years)
root = fileparts(mfilename('fullpath'));
h = readtable(fullfile(root,'input','hardware.csv'),'TextType','string');
b = readtable(fullfile(root,'input',"power_modes_"+cfg.cluster+".csv"),'TextType','string');
[found,order] = ismember(names,h.Vehicle);
assert(all(found),'EPS hardware table is missing a vehicle.');
h = h(order,:);
c = struct('names',names,'initialSoC',cfg.initialSoC,'solarFlux_W_m2',1361, ...
    'cellEfficiency_BOL',0.32,'packing',0.85,'solarEOL',0.80,'pmadEfficiency',0.90, ...
    'batteryEOL',0.80,'maxDoD',0.35,'batteryChargeEfficiency',0.85, ...
    'batteryDischargeEfficiency',0.90,'interfaceEfficiencyAK',0.90, ...
    'interfaceEfficiencyKR',0.90,'maxActiveKikiPorts',2,'atlanteRelocationPVFraction',0);
if ~isempty(plan.issues), cfg.journeyFraction=0; end
c.life = struct('years',years,'daysPerYear',365.25,'journeyFraction',cfg.journeyFraction, ...
    'highSoC',0.90,'binEdges_pct',[0 5 10 15 20 25 30 35 40 100], ...
    'referenceDoD',[.20 .20 .30 .30 .30 .30]);
events = plan.events;
links = plan.links;
events.Start_min=events.Start_s/60; events.End_min=events.End_s/60;
links.Start_min=links.Start_s/60; links.End_min=links.End_s/60;
fraction(~isfinite(fraction))=0;
if cfg.pvUsesAspect
    projection = max(0,cosd(aspect));
    projection(~isfinite(projection))=0;
    fraction = fraction.*projection;
    h.PointingFactor(:)=1;
end
run = FABER_environment_eps(c,h,b,events,links,edges/60,fraction);
kikiChargeLoad = b.Mean_W(b.Vehicle=="Kiki" & b.Mode=="K7")/1000;
rangoChargeLoad = b.Mean_W(b.Vehicle=="Rango" & b.Mode=="R5")/1000;
rangoTransportLoad = b.Mean_W(b.Vehicle=="Rango" & b.Mode=="R8")/1000;
run.assemblyChargeNet_kW = [cfg.kikiPort_kW-kikiChargeLoad- ...
    cfg.retainRangoForRelocations*rangoTransportLoad/c.interfaceEfficiencyKR, ...
    cfg.rangoPort_kW-rangoChargeLoad];
served = sum(run.served.*run.dt_h,1);
unserved = sum(run.unserved.*run.dt_h,1);
missing = sum(run.missingLink.*run.dt_h,1);
families = ["Taxi","Atlante","Kiki","Kiki","Rango","Rango"];
safeModes = ["T8","A6","K9","K9","R6","R6"];
reserve = zeros(1,numel(names));
for v=1:numel(names)
    safeLoad = b.Mean_W(b.Vehicle==families(v) & b.Mode==safeModes(v))/1000;
    reserve(v) = safeLoad*cfg.safeReserve_min/60/c.batteryDischargeEfficiency;
end
reserveMargin = min(run.E,[],1)-run.Emin-reserve;
summary = table(names',min(run.soc,[],1)',run.soc(end,:)',sum(run.pv.*run.dt_h,1)', ...
    served',unserved',missing',reserveMargin',(unserved<1e-7 & missing<1e-7 & reserveMargin>=-1e-7)', ...
    'VariableNames',{'Vehicle','MinSoC_pct','EndSoC_pct','SolarEnergy_kWh', ...
        'ServedEnergy_kWh','UnservedEnergy_kWh','MissingCharge_kWh','MinReserveMargin_kWh','Pass'});
summary.EndTargetSoC_pct = repmat(100*cfg.chargeTargetSoC,numel(names),1);
summary.EndTargetPass = summary.EndSoC_pct>=summary.EndTargetSoC_pct-1e-7;
mid = (edges(1:end-1)+edges(2:end))/2;
charge = links(links.Kind=="Charge",{'Sender','Receiver','Start_s','End_s','Received_kW'});
charge.Requested_kWh = zeros(height(charge),1); charge.Delivered_kWh = charge.Requested_kWh;
charge.StartSoC_pct = charge.Requested_kWh; charge.EndSoC_pct = charge.Requested_kWh;
charge.TargetShortfall_kWh = charge.Requested_kWh;
charge.BatteryEnergyGained_kWh = charge.Requested_kWh;
charge.AdditionalChargeAtRatedPort_min = nan(height(charge),1);
for i=1:height(charge)
    j = find(names==charge.Receiver(i),1);
    a = find(edges==charge.Start_s(i),1);
    z = find(edges==charge.End_s(i),1);
    keep = mid>=charge.Start_s(i) & mid<charge.End_s(i);
    charge.Requested_kWh(i) = charge.Received_kW(i)*(charge.End_s(i)-charge.Start_s(i))/3600;
    charge.Delivered_kWh(i) = sum(run.rx(keep,j).*run.dt_h(keep));
    charge.StartSoC_pct(i) = run.soc(a,j);
    charge.EndSoC_pct(i) = run.soc(z,j);
    charge.BatteryEnergyGained_kWh(i) = run.E(z,j)-run.E(a,j);
    charge.TargetShortfall_kWh(i) = max(0,cfg.chargeTargetSoC*run.EOL(j)-run.E(z,j));
    net = max(0,charge.Received_kW(i)-mean(run.own(keep,j)));
    rate = c.batteryChargeEfficiency*min(net,h.ChargeC_per_h(j)*run.EOL(j));
    if charge.TargetShortfall_kWh(i)==0
        charge.AdditionalChargeAtRatedPort_min(i) = 0;
    elseif rate>0
        charge.AdditionalChargeAtRatedPort_min(i) = 60*charge.TargetShortfall_kWh(i)/rate;
    end
end
charge.Duration_min = (charge.End_s-charge.Start_s)/60;
charge = withUTC(charge,cfg.epochUTC);
[life,cycles,histogram] = eps_battery_lifetime(c,h,run,"Environment timeline",all(summary.Pass));
end

function events = eventLighting(events,cfg,edges,names,fraction,aspect)
events = withUTC(events,cfg.epochUTC);
dt = diff(edges)/60;
mid = (edges(1:end-1)+edges(2:end))/2;
N = height(events);
events.Sunlight_min = zeros(N,1);
events.Eclipse_min = zeros(N,1);
events.MinSunAspect_deg = nan(N,1);
events.MaxSunAspect_deg = nan(N,1);
for i=1:N
    v = find(names==events.Vehicle(i),1);
    k = mid>=events.Start_s(i) & mid<events.End_s(i);
    f = fraction(k,v);
    d = dt(k);
    events.Sunlight_min(i) = sum(d(f>=1-1e-9));
    events.Eclipse_min(i) = sum(d(f<1-1e-9));
    angles = aspect(k,v);
    if any(isfinite(angles))
        events.MinSunAspect_deg(i) = min(angles,[],'omitnan');
        events.MaxSunAspect_deg(i) = max(angles,[],'omitnan');
    end
end
end

function summary = dutyCycle(run,names,fraction)
summary = table();
for v=1:numel(names)
    for mode = unique(run.mode(:,v),'stable')'
        k = run.mode(:,v)==mode;
        h = run.dt_h(k);
        light = fraction(k,v);
        summary = [summary;table(names(v),mode,sum(h),sum(h(light>=1-1e-9)), ...
            sum(h(light<1-1e-9)),sum(h.*run.own(k,v)),sum(h.*run.pv(k,v)), ...
            sum(h.*run.batOut(k,v)),sum(h.*run.batIn(k,v)), ...
            sum(h.*run.unserved(k,v)),max(run.peak(k,v)), ...
            'VariableNames',{'Vehicle','Mode','Duration_h','Sunlight_h','Eclipse_h', ...
            'OwnLoad_kWh','SolarEnergy_kWh','BatteryOut_kWh','BatteryIn_kWh','Unserved_kWh','PeakModeLoad_kW'})]; %#ok<AGROW>
    end
end
end

function t = withUTC(t,epoch)
t.StartUTC = epoch+seconds(t.Start_s);
t.EndUTC = epoch+seconds(t.End_s);
end

function t = retrievalUTC(t,epoch)
for label = ["Start","Departure","Capture","Handover"]
    t.(label+"UTC") = epoch+seconds(t.(label+"_s"));
end
end

function saveOutputs(out)
folder = out.options.outputFolder;
if ~isfolder(folder), mkdir(folder); end
fields = ["events","links","holds","openItems","lightingTrace","lightingSummary", ...
    "eclipseIntervals","toolboxEclipseIntervals","referenceEvents","referenceLightingTrace", ...
    "epsSummary","dutyCycle", ...
    "chargeSummary","batteryLifetime","batteryCycles","batteryHistogram", ...
    "retrievalSummary","referenceRetrievalSummary","assemblySorties"];
for f=fields
    if width(out.(f))>0
        writetable(out.(f),fullfile(folder,"Environment_"+f+".csv"));
    end
end
if isempty(out.figures)
    return;
end
% Copy each tab into a hidden figure for a separate PNG export.
exportFigure = figure('Visible','off','WindowStyle','normal','Position',[100 100 960 680]);
cleanup = onCleanup(@() delete(exportFigure)); %#ok<NASGU>
theme(exportFigure,'light');
for f = out.figureTabs(:)'
    clf(exportFigure);
    axesAndLegend = [findobj(f,'Type','axes');findobj(f,'Type','legend')];
    copyobj(axesAndLegend,exportFigure);
    exportgraphics(exportFigure,fullfile(folder,string(f.Tag)+".png"),'Resolution',150);
end
end
