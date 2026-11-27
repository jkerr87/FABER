function [study,budget] = main_FABER_kiki_winter(outputFolder,kozFactor,budget,cadenceFolder,tripsPerModule,options)
% Run two winter deliveries and plot the Taxi and assembly duty cycles.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root,fullfile(root,'Environment'));

if nargin<2 || isempty(kozFactor), kozFactor = 1; end

validateattributes(kozFactor,{'numeric'},{'scalar','real','finite','>=',1});

if nargin<1 || isempty(outputFolder)
    outputFolder = fullfile(pwd,sprintf('Kiki_Winter_KOZx%g',kozFactor));
end

if nargin>=4 && ~isempty(cadenceFolder)
    if nargin<5 || isempty(tripsPerModule), tripsPerModule = 4; end
    if nargin<6, options = struct(); end
    validateattributes(tripsPerModule,{'numeric'},{'scalar','integer','positive','finite'});
    [study.roundTrips,study.charging,study.comparison,study.figure,study.figureTabs, ...
        study.events,study.options,study.deliveryCapacity] = ...
        assemblyCadence(cadenceFolder,tripsPerModule,options);
    if ~isfolder(outputFolder), mkdir(outputFolder); end
    if ~isempty(study.figureTabs)
        ax = findobj(study.figureTabs(1),'Type','axes');
        exportgraphics(ax,fullfile(outputFolder,'Environment_EPS_Duty_Cycle.png'),'Resolution',160);
    end
    writetable(study.roundTrips,fullfile(outputFolder,'Assembly_Round_Trips.csv'));
    writetable(study.charging,fullfile(outputFolder,'Assembly_Charging.csv'));
    writetable(study.events,fullfile(outputFolder,'Assembly_Timeline.csv'));
    writetable(study.deliveryCapacity,fullfile(outputFolder,'Assembly_Delivery_Capacity.csv'));
    return;
end

if nargin<3 || isempty(budget)
    budget = MA_vehicle_budget_phaseA(struct('N_RPO_MC',30, ...
        'captureWait_min',45,'handoverWait_min',10,'veicLoading_min',5, ...
        'outboundAtlanteRAAN_deg',180,'runAnnualSK',false, ...
        'includeMomentumUnloading',false));
elseif ~isfield(budget.cfg,'kikiTransferPropulsion')
    defaults = FABER_phaseA_defaults();
    budget.cfg.kikiTransferPropulsion = defaults.kikiTransferPropulsion;
end


cfg = struct('epochUTC',datetime(2038,2,1,'TimeZone','UTC'), ...
    'useToolbox',false,'captureActive_min',5,'assemblyModuleIndices',[1 2], ...
    'assemblyExtraHolds',2,'assemblyHoldKOZFactor',kozFactor, ...
    'assemblyTimelineDetails',true,'makePlots',true,'standbyAfter_min',120, ...
    'outputFolder',outputFolder);

cfg.retrievals = struct('moduleId',{'A','B'},'phaseAngle_deg',{20,10},'cargoCount',1);

study = run_faber_environment(budget,cfg);

mission = findobj(study.figureTabs,'Tag','Environment_EPS_Duty_Cycle');

if isscalar(mission), mission.Parent.SelectedTab = mission; end

end

function [roundTrips,charging,comparison,window,tabs,timeline,options,deliveryCapacity] = assemblyCadence(cadenceFolder,tripsPerModule,options)
% Screen repeated assembly work using the winter schedule and EPS exports.
% This replays durations and battery withdrawals, without propagating new RPO epochs.
readExport = @(folder,name) readtable(fullfile(folder,"Environment_"+name+".csv"), ...
    'TextType','string');
defaults = struct('respectSunlight',true,'leaveRangoOnHost',false, ...
    'rangoAssembly_min',32.17,'rangoChargeConnect_min',5,'rangoWorkBetweenCharges_min',90, ...
    'makePlots',true,'independentKikiPorts',true,'season',"Winter");
for field = string(fieldnames(options))'
    assert(isfield(defaults,field),'Unknown assembly option %s.',field);
    defaults.(field) = options.(field);
end
options = defaults;
options.season = string(validatestring(options.season,{'Winter','Summer'}));
if options.season=="Summer"
    fiveModuleFolder = fullfile(cadenceFolder,'Seasonal_Five','Summer');
    calibrationFolder = fiveModuleFolder;
else
    fiveModuleFolder = fullfile(cadenceFolder,'Winter_On_Demand_Five');
    calibrationFolder = fullfile(cadenceFolder,'Kiki_Winter_Duty_Cycles');
end
baseline = readExport(fiveModuleFolder,"assemblySorties");
events = readExport(fiveModuleFolder,"events");
eps = readExport(calibrationFolder,"epsSummary");
calibrationCharges = readExport(calibrationFolder,"chargeSummary");
validateattributes(options.rangoAssembly_min,{'numeric'},{'scalar','finite','positive'});
validateattributes(options.rangoChargeConnect_min,{'numeric'},{'scalar','finite','nonnegative'});
validateattributes(options.rangoWorkBetweenCharges_min,{'numeric'}, ...
    {'scalar','real','finite','>=',options.rangoAssembly_min});
assert(height(baseline)==5 && isequal(baseline.ModuleIndex,(1:5)'), ...
    'Supply a winter schedule with five completed modules.');
if options.season=="Summer"
    assert(~any(ismember(calibrationCharges.Receiver,["Kiki 2","Rango 2"])), ...
        'Summer calibration requires an uncharged second pair.');
    calibrationTrips = nnz(baseline.Pair==2);
    assert(calibrationTrips>0,'Summer calibration requires completed trips by pair 2.');
else
    assert(isempty(calibrationCharges), ...
        'Use the uncharged two-module winter run to estimate battery withdrawals.');
    calibrationTrips = 1;
end

input = fullfile(fileparts(fileparts(mfilename('fullpath'))),'input');
h = readtable(fullfile(input,'hardware.csv'),'TextType','string');
p = readtable(fullfile(input,'power_modes_ISS.csv'),'TextType','string');
vehicles = ["Kiki 1","Rango 1";"Kiki 2","Rango 2"];
capacity = zeros(2);
draw = zeros(2);
for pair = 1:2
    [~,hardwareRow] = ismember(vehicles(pair,:),h.Vehicle);
    calibrationPair = pair;
    if options.season=="Summer", calibrationPair = 2; end
    [~,epsRow] = ismember(vehicles(calibrationPair,:),eps.Vehicle);
    capacity(pair,:) = 0.80*h.Battery_kWh(hardwareRow)';
    draw(pair,:) = capacity(pair,:).*(1-eps.EndSoC_pct(epsRow)'/100)/calibrationTrips;
end
load = @(mode) p.Mean_W(p.Mode==mode)/1000;
kikiSoloRate = 0.85*(3-load("K7"));
rangoRate = 0.85*(1.5-load("R5"));
hostRate = (load("K6")+1.5/0.90)/0.90;
% Match the environmental EPS floor and 15-minute safe-mode reserve.
reserve = 0.65*capacity+[load("K9") load("R6")]*0.25/0.90;
assert(all(capacity-draw>=reserve,'all'), ...
    'A full battery cannot complete the representative round trip above reserve.');
energy = capacity;
connect_min = options.rangoChargeConnect_min;

replay = events(events.Vehicle=="Taxi",:);
if ~options.respectSunlight, replay = removeSunWaits(replay,"Taxi"); end
handover = replay.End_s(replay.Activity=="Taxi to Atlante handover")/60;
assert(numel(handover)==height(baseline),'Supply a schedule with all Taxi handovers.');
count = height(baseline)*tripsPerModule;
roundTrips = array2table(zeros(count,11),'VariableNames', ...
    {'Module','RoundTrip','Pair','Start_min','BerthEnd_min','ChargeDelay_min', ...
    'InstallationStart_min','InstallationEnd_min','RangoLeftOnHost','WaitForCargo_min','PortionLoadedInParallel_min'});
charging = table('Size',[0 7], ...
    'VariableTypes',{'double','double','double','double','string','string','string'}, ...
    'VariableNames',{'TripIndex','Pair','ChargeStart_min','RechargeEnd_min','Sender','Receiver','Location'});
siteFree = 0;
campaignEnd = 0;
berth = zeros(2,1);
residentSince = nan(2,1);
workSinceCharge = zeros(2,1);
for j = 1:count
    module = ceil(j/tripsPerModule);
    trip = j-(module-1)*tripsPerModule;
    pair = 1+mod(j-1,2);
    row = baseline(module,:);
    source = ["Kiki "+row.Pair,"Rango "+row.Pair];
    target = vehicles(pair,:);
    a = row.LoadStart_s;
    b = row.BerthEnd_s;
    e = events(ismember(events.Vehicle,source) & events.Start_s>=a-1e-6 & ...
        events.End_s<=b+1e-6 & startsWith(events.Activity, ...
        ["Assembly "+module+" ","Attached wait","Eclipse wait"]),:);
    e.Start_s = max(e.Start_s-a,0);
    e.End_s = min(e.End_s-a,b-a);
    if ~options.respectSunlight, e = removeSunWaits(e,source(1)); end

    resident = isfinite(residentSince(pair));
    leave = options.leaveRangoOnHost && trip+2<=tripsPerModule;
    if resident
        deploy = contains(e.Activity,"deploy Rango") & e.Vehicle==source(1);
        deployStart = min(e.Start_s(deploy));
        deployEnd = max(e.End_s(deploy));
        e(e.Vehicle==source(2) & e.End_s<=deployEnd+1e-6,:) = [];
        e = compress(e,deployStart,deployEnd);
        e(contains(e.Activity,"deploy Rango"),:) = [];
    end
    work = e.Mode=="R4";
    assert(nnz(work)==1,'Expected one installation per round-trip template.');
    workStart = e.Start_s(work);
    workEnd = e.End_s(work);
    workDelta = options.rangoAssembly_min-(workEnd-workStart)/60;
    e = delayAfter(e,workEnd,60*workDelta);
    workEnd = workStart+60*options.rangoAssembly_min;
    e.End_s(work | contains(e.Activity,"standby during installation")) = workEnd;

    if leave
        retreat = min(e.Start_s(contains(e.Activity,"unloaded departure")));
        recover = contains(e.Activity,"recover Rango");
        e = delayAfter(e,retreat,workStart-retreat);
        e(recover | contains(e.Activity,"standby during installation") | ...
            (e.Vehicle==source(2) & contains(e.Activity,["unloaded","retained return","retained at Atlante berth"])),:) = [];
        e.End_s(e.Mode=="R4") = workEnd;
    end
    loadAtAtlante = sum(e.End_s(e.Mode=="K3")-e.Start_s(e.Mode=="K3"))/60;
    approach = contains(e.Activity,"loaded terminal arrival") & ...
        ~contains(e.Activity,"unloaded") & e.Vehicle==source(1);
    assert(any(approach),'Missing loaded terminal arrival in assembly template.');
    hostApproachOffset = min(e.Start_s(approach))/60;
    previousFinish = campaignEnd;
    if options.independentKikiPorts
        % Separate port turnaround from the shared host approach reservation.
        available = max(berth(pair),siteFree-hostApproachOffset);
    else
        available = max(berth(pair),siteFree-loadAtAtlante);
    end
    earliest = max(available,handover(module));
    start = earliest;
    demand = max(0,draw(pair,:)+workDelta*[load("K11") load("R4")]/(60*0.90));
    if resident || leave
        % Account for the omitted transport and installation standby using
        % the active modes in this shuttle instead of a carried-pair draw.
        for v = 1:2
            active = e.Vehicle==source(v);
            load_kW = arrayfun(@(mode) load(mode),e.Mode(active));
            demand(v) = sum(load_kW.*max(e.End_s(active)-e.Start_s(active),0))/(3600*0.90);
        end
    end
    % Use measured winter withdrawals as a conservative allowance for shuttles.
    % Also account for Rango's idle load while it remains on the structure.
    rangoBefore = energy(pair,2);
    if resident
        rangoBefore = rangoBefore-(start+workStart/60-residentSince(pair))*load("R7")/(60*0.90);
    end
    targetEnergy = reserve(pair,2)+demand(2);
    if leave
        % Cover the next two shuttle intervals plus a possible Kiki
        % top-up delay while Rango remains on the structure.
        idle_min = 2*max(e.End_s)/60+45;
        targetEnergy = min(capacity(pair,2),targetEnergy+idle_min*load("R7")/(60*0.90));
    end
    if workSinceCharge(pair)+options.rangoAssembly_min>options.rangoWorkBetweenCharges_min
        % Keep installations intact by servicing Rango immediately before
        % the next block if that block would cross the work limit.
        workEnergy = workSinceCharge(pair)*load("R4")/(60*0.90);
        targetEnergy = min(capacity(pair,2),max(targetEnergy,rangoBefore+workEnergy));
    end
    % Kiki must retain enough energy to return to Atlante after supplying
    % Rango. Limit optional host top-up energy to what a full Kiki can supply.
    maxShared = max(0,(60*(capacity(pair,1)-reserve(pair,1)-demand(1)) ...
        -connect_min*load("K6c")/0.90)/hostRate);
    shared = min(maxShared,60*max(0,targetEnergy-rangoBefore)/rangoRate);
    kikiNeed = demand(1)+(shared*hostRate+connect_min*load("K6c")/0.90*(shared>0))/60;
    if energy(pair,1)-kikiNeed<reserve(pair,1)
        topup = ceil(60*(capacity(pair,1)-energy(pair,1))/kikiSoloRate);
        start = max(earliest,berth(pair)+topup);
        if topup>0
            charging(end+1,:) = {j,pair,start-topup,start,"Atlante",target(1),"Atlante berth"};
            addEvent(target(1),"Recharge at Atlante","K7",start-topup,start);
        end
        energy(pair,1) = capacity(pair,1);
        if resident && start>earliest
            extraDraw = (start-earliest)*load("R7")/(60*0.90);
            rangoBefore = rangoBefore-extraDraw;
            if shared>0
                extra = min(maxShared-shared,60*extraDraw/rangoRate);
                shared = shared+extra;
                kikiNeed = kikiNeed+extra*hostRate/60;
            end
        end
    end
    assert(rangoBefore>=reserve(pair,2), ...
        'Rango reaches reserve before Kiki can charge it on module %d trip %d.',module,trip);
    host = start+workStart/60;
    if resident
        receive = min(e.Start_s(e.Mode=="R2"))/60+start;
        addEvent(target(2),"Rango attached standby on host","R7",residentSince(pair),receive);
    end
    if shared>0
        % Work and any departure wait until Kiki finishes the stationary host service.
        span = connect_min+shared;
        e = delayAfter(e,workStart,60*span);
        addEvent(target(1),"Kiki attached at host charging Rango","K6",host,host+span);
        addEvent(target(2),"Rango charging from Kiki at host","R5",host+connect_min,host+span);
        charging(end+1,:) = {j,pair,host+connect_min,host+span,target(1),target(2),"Host structure"};
        rangoBefore = min(capacity(pair,2),rangoBefore+shared*rangoRate/60);
        workSinceCharge(pair) = 0;
    end
    workSinceCharge(pair) = workSinceCharge(pair)+options.rangoAssembly_min;
    assert(workSinceCharge(pair)<=options.rangoWorkBetweenCharges_min+1e-9, ...
        'Kiki cannot recharge Rango before the assembly limit on module %d trip %d.',module,trip);
    energy(pair,:) = [energy(pair,1)-kikiNeed,rangoBefore-demand(2)];
    assert(all(energy(pair,:)>=reserve(pair,:)-1e-9), ...
        ['Battery reserve failed on module %d trip %d. Kiki %.2f kWh (floor %.2f), ' ...
        'Rango %.2f kWh (floor %.2f).'],module,trip,energy(pair,1),reserve(pair,1), ...
        energy(pair,2),reserve(pair,2));
    work = e.Mode=="R4";
    installStart = start+e.Start_s(work)/60;
    installEnd = start+e.End_s(work)/60;
    finish = start+max(e.End_s(e.Vehicle==source(1)))/60;
    berth(pair) = finish;
    residentSince(pair) = NaN;
    if leave, residentSince(pair) = installEnd; end
    roundTrips(j,:) = {module,trip,pair,start,finish,start-earliest,installStart,installEnd,leave, ...
        max(0,handover(module)-available),max(0,min(loadAtAtlante,previousFinish-start))};
    if options.independentKikiPorts
        exit = endsWith(e.Activity,"unloaded departure") & e.Vehicle==source(1);
        assert(any(exit),'Missing host departure in assembly template.');
        siteFree = max(start+max(e.End_s(exit))/60,installEnd);
    else
        siteFree = max(finish,installEnd);
    end
    campaignEnd = max([campaignEnd,finish,installEnd]);
    for v = 1:2
        e.Vehicle(e.Vehicle==source(v)) = target(v);
    end
    e.Start_s = e.Start_s+60*start;
    e.End_s = e.End_s+60*start;
    replay = [replay;e(e.End_s>e.Start_s,:)]; %#ok<AGROW>
end
comparison = table(count,campaignEnd/1440,sum(roundTrips.ChargeDelay_min),height(charging), ...
    sum(charging.Location=="Host structure"),sum(roundTrips.InstallationEnd_min-roundTrips.InstallationStart_min), ...
    'VariableNames',{'ScreenRoundTrips','ScreenCompletion_days','AtlanteChargeDelay_min', ...
    'RechargeSessions','RangoHostCharges','RangoInstallation_min'});
disp(comparison);
deliveryCapacity = table();
for module = 1:numel(handover)-1
    rows = roundTrips(roundTrips.Module==module,:);
    complete = max(rows.BerthEnd_min,rows.InstallationEnd_min);
    cutoff = handover(module+1);
    line = table(module,height(rows),nnz(complete<=cutoff),cutoff-max(complete), ...
        sum(rows.ChargeDelay_min), 'VariableNames',{'Module','ScheduledTrips', ...
        'TripsBeforeNextHandover','SlackToNextHandover_min','ChargeDelay_min'});
    deliveryCapacity = [deliveryCapacity;line]; %#ok<AGROW>
end
out = struct('events',replay,'roundTrips',roundTrips,'charging',charging, ...
    'handover_min',handover,'season',options.season);
window = gobjects(0);
tabs = gobjects(0);
if options.makePlots, [window,tabs] = FABER_environment_plots(out); end
timeline = sortrows(replay(:,{'Vehicle','Activity','Mode','Start_s','End_s'}),'Start_s');
timeline.Duration_min = (timeline.End_s-timeline.Start_s)/60;

    function addEvent(vehicle,activity,mode,start,finish)
        if finish<=start, return; end
        event = events(1,:);
        event.Vehicle = vehicle;
        event.Activity = activity;
        event.Mode = mode;
        event.Start_s = start*60;
        event.End_s = finish*60;
        event.SunRequired = false;
        event.DeltaV_mps = 0;
        replay = [replay;event];
    end
end

function events = removeSunWaits(events,vehicle)
waits = startsWith(events.Activity,["Attached wait","Eclipse wait"]);
rows = find(waits & events.Vehicle==vehicle);
for i = flip(rows)'
    events = compress(events,events.Start_s(i),events.End_s(i));
end
events(waits,:) = [];
events.SunRequired(:) = false;
end

function events = compress(events,first,last)
events.Start_s = events.Start_s-min(max(events.Start_s-first,0),last-first);
events.End_s = events.End_s-min(max(events.End_s-first,0),last-first);
end

function events = delayAfter(events,cut,delta)
events.Start_s(events.Start_s>=cut-1e-7) = events.Start_s(events.Start_s>=cut-1e-7)+delta;
events.End_s(events.End_s>cut+1e-7) = events.End_s(events.End_s>cut+1e-7)+delta;
end
