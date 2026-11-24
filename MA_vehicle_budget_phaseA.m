function out = MA_vehicle_budget_phaseA(cfg)
% Build the per-event maneuver ledger and report stationkeeping over the configured exposure.
root = fileparts(mfilename('fullpath'));
addpath(genpath(root));
addpath(fullfile(root,'MonteCarlo','FABER_safety_v3_nine_phase'),'-begin');
if nargin < 1 || isempty(cfg)
    cfg = FABER_phaseA_defaults();
else
    defaults = FABER_phaseA_defaults();
    if isfield(cfg,'orbitRadialAmplitude_m')
        defaults.taxiParkingEllipse.ellipseRadialAmplitude_m = cfg.orbitRadialAmplitude_m;
        defaults.atlanteParkingEllipse.ellipseRadialAmplitude_m = cfg.orbitRadialAmplitude_m;
    end

    for name = {'taxiParkingEllipse','atlanteParkingEllipse'}
        field = name{1};
        if isfield(cfg,field)
            values = fieldnames(cfg.(field));
            for k = 1:numel(values)
                defaults.(field).(values{k}) = cfg.(field).(values{k});
            end

            cfg = rmfield(cfg,field);
        end
    end

    fields = fieldnames(cfg);
    for k = 1:numel(fields)
        defaults.(fields{k}) = cfg.(fields{k});
    end

    if ~isfield(cfg,'taxiParkingYears')
        defaults.taxiParkingYears = defaults.missionYears;
    end

    cfg = defaults;
end

if exist('astroConstants','file') == 2
    mu = astroConstants(13);
    Re = astroConstants(23);
else
    mu = 398600.4418;
    Re = 6378.1363;
end

r0 = Re + cfg.h_target_km;
assert(exist('truss_platform','file') ~= 0, ...
    'Missing truss_platform.m');
run(which('truss_platform'));
Nreloc = size(modulePoints,1);
timesMin = [5 10 15 20 25 30 35 40 45 50];
hostPoint = [0;-250;0];

% Assemble individual maneuvers before applying grouped event margins.
rows = cell(0,12);
[dva,dvb,tinj] = hohmann(mu, ...
    r0 + cfg.launchAltError_km,r0);
vInjection = sqrt(mu/(r0+cfg.launchAltError_km))*1000;
vTarget = sqrt(mu/r0)*1000;
DV_alt = abs(dva)+abs(dvb);
DV_inc = 2*vInjection* ...
    sind(cfg.launchIncError_deg/2);
planeAngle = 2*asind(abs( ...
    sind(cfg.i_target_deg)* ...
    sind(cfg.launchRAANError_deg/2)));
DV_RAAN = 2*vTarget*sind(planeAngle/2);
rows = add(rows,'Taxi','Deployment', ...
    'Injection altitude','Baseline',1, ...
    DV_alt,tinj/60,'Estimated', ...
    'MA_deltav',true);
rows = add(rows,'Taxi','Deployment', ...
    'Injection inclination','Baseline',1, ...
    DV_inc,0,'Estimated', ...
    'MA_deltav',true);
rows = add(rows,'Taxi','Deployment', ...
    'Injection RAAN','Baseline',1, ...
    DV_RAAN,0,'Estimated', ...
    'MA_deltav',true);
rows = add(rows,'ATLANTE','Deployment', ...
    'Stack injection altitude','Baseline',1, ...
    DV_alt,tinj/60,'Estimated', ...
    'ATLANTE + Kiki + Rango',true);
rows = add(rows,'ATLANTE','Deployment', ...
    'Stack injection inclination','Baseline',1, ...
    DV_inc,0,'Estimated', ...
    'ATLANTE + Kiki + Rango',true);
rows = add(rows,'ATLANTE','Deployment', ...
    'Stack injection RAAN','Baseline',1, ...
    DV_RAAN,0,'Estimated', ...
    'ATLANTE + Kiki + Rango',true);


[ca,cb,tc] = hohmann(mu, ...
    r0+cfg.hClientSpread_km, ...
    r0-cfg.hClientSpread_km);
iHigh = cfg.i_target_deg + ...
    cfg.iClientSpread_deg + cfg.launchIncError_deg;
iLow = cfg.i_target_deg - ...
    cfg.iClientSpread_deg - cfg.launchIncError_deg;
cosGamma = cosd(iHigh)*cosd(iLow) + ...
    sind(iHigh)*sind(iLow)* ...
    cosd(2*cfg.clientRAANSpread_deg);
gamma = acosd(max(-1,min(1,cosGamma)));
vClient = sqrt(mu/(r0-cfg.hClientSpread_km))*1000;
DV_plane = 2*vClient*sind(gamma/2);
DV_altClient = abs(ca)+abs(cb);
for vehicle = {'Taxi','ATLANTE'}
    veh = vehicle{1};
    rows = add(rows,veh,'Inter-client', ...
        'Extreme altitude transfer','Baseline', ...
        1,DV_altClient,tc/60, ...
        'Estimated','Worst-case client envelope',true);
    rows = add(rows,veh,'Inter-client', ...
        'Extreme plane transfer','Baseline', ...
        1,DV_plane,0, ...
        'Estimated','Combined plane-angle change',true);
    rows = add(rows,veh,'Inter-client', ...
        'RAAN reserve','Baseline',1,0,0, ...
        'Allocated','PM3 RAAN reserve per inter-client transfer',true);
end

% Use the same module orbit for both RPO and Lambert targeting.
hModule = cfg.h_target_km + cfg.moduleAltitudeOffset_km;
iModule = cfg.i_target_deg + cfg.moduleInclinationOffset_deg;
validateattributes(cfg.captureWait_min,{'numeric'},{'scalar','real'});
assert(isnan(cfg.captureWait_min) || (isfinite(cfg.captureWait_min) && ...
    cfg.captureWait_min>=0),'Capture dwell must be nonnegative or NaN.');
validateattributes(cfg.farFieldIsp_s,{'numeric'},{'scalar','real','finite','positive'});
outbound = FABER_outbound_lambert(mu,Re,cfg);
rows = add(rows,'Taxi','Outbound far', ...
    'Lambert departure burn','Baseline',1, ...
    norm(outbound.departureBurn_mps),outbound.transferTime_s/60, ...
    'Estimated','Two-body Lambert transfer from parking handoff',true);
rows = add(rows,'Taxi','Outbound far', ...
    'RPO entry velocity match','Baseline',1, ...
    norm(outbound.arrivalBurn_mps),0, ...
    'Estimated','Inertial velocity match at the RPO gate',true);

% Initialize paired RPO from the propagated outbound gate state.
mcCfg = FABER_phaseA_rpo_config(cfg,mu,Re,outbound);
RPO_MC = FABER_RPO_9phase_MC(mcCfg);
phases = RPO_MC.phaseNames;
valid = find(RPO_MC.valid(:,1));
idxCycle = NaN;
phaseDV = nan(2,numel(phases));
phaseTOF = phaseDV;
status = 'Open';
clientStatus = status;
returnTransfer = [];
if ~isempty(valid)
    cycleDV = RPO_MC.CycleDV_mps(valid,1);
    [~,nearP95] = min(abs(cycleDV-prctile(cycleDV,95)));
    idxCycle = valid(nearP95);
    phaseDV = squeeze(RPO_MC.DV_phase_mps(idxCycle,1,:,:));
    phaseTOF = squeeze(RPO_MC.TOF_phase_min(idxCycle,1,:,:));
    status = 'Paired MC P95-DV representative';
    clientStatus = status;
    returnTransfer = RPO_MC.records{idxCycle,1}.ReturnTransfer;
else
    % Preserve an outbound estimate when no complete retrieval is feasible.
    clients = RPO_MC.clientCandidates(:,1);
    clients = clients(~cellfun(@isempty,clients));
    if ~isempty(clients)
        dv = cellfun(@(client) client.dv,clients);
        [~,idxClient] = min(abs(dv-prctile(dv,95)));
        phaseDV(1,:) = clients{idxClient}.dvPhase;
        phaseTOF(1,:) = clients{idxClient}.tofPhase;
        clientStatus = 'Outbound-only MC P95-DV representative';
    end
end

DV_client_phase = phaseDV(1,:);
DV_atlante_phase = phaseDV(2,:);
TOF_client_phase = phaseTOF(1,:);
TOF_atlante_phase = phaseTOF(2,:);
for k = 1:numel(phases)
    rows = add(rows,'Taxi','Outbound RPO',phases{k},'Baseline', ...
        1,DV_client_phase(k),TOF_client_phase(k), ...
        clientStatus,'Nominal client approach feasibility checks',true);
    rows = add(rows,'Taxi','Inbound RPO',phases{k},'Baseline', ...
        1,DV_atlante_phase(k),TOF_atlante_phase(k), ...
        status,'Same paired cycle as outbound',true);
end

% Keep return costs open until a complete pair has been selected.
returnTOF_min = NaN;
returnDV = [NaN NaN];
if ~isempty(returnTransfer)
    returnTOF_min = returnTransfer.transferTime_s/60;
    returnDV = [norm(returnTransfer.departureBurn_mps),norm(returnTransfer.arrivalBurn_mps)];
end

rows = add(rows,'Taxi','Inbound far', ...
    'Lambert departure burn','Baseline',1, ...
    returnDV(1),returnTOF_min,status,'Captured stack after dwell',true);
rows = add(rows,'Taxi','Inbound far', ...
    'RPO entry velocity match','Baseline',1, ...
    returnDV(2),0,status,'Inertial velocity match relative to Atlante',true);

RPO_MC_Table = table( ...
    string(phases(:)), ...
    DV_client_phase(:), ...
    TOF_client_phase(:), ...
    DV_atlante_phase(:), ...
    TOF_atlante_phase(:), ...
    'VariableNames',{ ...
        'Phase', ...
        'Client_DV_mps', ...
        'Client_TOF_min', ...
        'ATLANTE_DV_mps', ...
        'ATLANTE_TOF_min'});
fprintf('\n========================================\n');
fprintf('TAXI MONTE CARLO RPO PHASE BUDGET\n');
fprintf('========================================\n');
disp(RPO_MC_Table);
RPO_MC.SelectedCycle = idxCycle;
RPO_MC.PhaseTable = RPO_MC_Table;
rows = add(rows,'Taxi','Payload', ...
    'Capture module','Baseline',1, ...
    NaN,cfg.captureWait_min, ...
    'Open','Capture dynamics not modelled',true);
rows = add(rows,'Taxi','Payload', ...
    'ATLANTE handover','Baseline',1, ...
    NaN,cfg.handoverWait_min, ...
    'Open','Handover dynamics not modelled',true);
opts = struct();
opts.methods = {'GeneralHCW','Vbar','Rbar'};
opts.gateDistance_m = 40;
timesMin = cfg.assemblyTransferTimes_min;
opts.exitTimesMin = 5:5:30;
opts.finalTimesMin = 5:5:30;
opts.maxExitSpeed_mps = 0.1;
opts.maxTransferSpeed_mps = 0.3;
opts.maxFinalSpeed_mps = 0.1;
opts.maxExitAccel_mps2 = inf;
opts.maxFinalAccel_mps2 = inf;
opts.printSummary = false;
Npts = 201;
% Assembly relocation optimization
opts.objective = 'weighted';

% Equal normalized importance of DV and TOF
opts.weightDV = 0.5;
opts.weightTOF = 0.5;

% Reference scales for objective normalization
opts.dvReference_mps = 0.1;
opts.tofReference_min = 60;

A = assembly_relocations_v4( ...
    mu,Re,cfg.h_target_km,201, ...
    timesMin,hostPoint,modulePoints, opts);
vehicles = {'Kiki','Rango'};

%% ================================================================
% ASSEMBLY RELOCATIONS - SELECTED OPTIMIZATION RESULTS
% ================================================================

% A contains the selected trajectories from assembly_relocations_v4.
% Selection uses opts.objective ('dv', 'tof', or 'weighted').

Nreloc = A.Nmodules;

TOF_outbound_min = A.totalTimePerDirection_min(:,1);
TOF_return_min = A.totalTimePerDirection_min(:,2);
DV_outbound_mps = A.DV_outbound;
DV_return_mps = A.DV_return;
tol = 1e-6;
assert(isequaln(isnan(TOF_outbound_min + TOF_return_min), ...
    isnan(A.TOF_perModule_min)));
assert(all(abs(TOF_outbound_min + TOF_return_min - ...
    A.TOF_perModule_min) < tol | isnan(A.TOF_perModule_min)));

%% VEHICLE-LEVEL RESULTS

vehicles = {'Kiki','Rango'};

AssemblyVehicleTable = table();

for j = 1:2

    indices = 1:Nreloc;
    if j == 2
        indices = [];
    end

    count = numel(indices);
    unresolved = nnz(~isfinite(A.DV_perModule(indices)));

    % Delta-V from selected weighted trajectories
    dvOut  = sum(DV_outbound_mps(indices));
    dvBack = sum(DV_return_mps(indices));

    % Complete TOF from selected weighted trajectories
    tofOut  = sum(TOF_outbound_min(indices));
    tofBack = sum(TOF_return_min(indices));

    % Total vehicle relocation requirements
    dvTotal  = dvOut + dvBack;
    tofTotal = tofOut + tofBack;

    % Store a separate results table
    newRow = table( ...
        string(vehicles{j}), ...
        count, ...
        unresolved, ...
        dvOut, ...
        dvBack, ...
        dvTotal, ...
        tofOut, ...
        tofBack, ...
        tofTotal, ...
        'VariableNames',{ ...
        'Vehicle', ...
        'Nmodules', ...
        'UnresolvedModules', ...
        'DV_Outbound_mps', ...
        'DV_Return_mps', ...
        'DV_Total_mps', ...
        'TOF_Outbound_min', ...
        'TOF_Return_min', ...
        'TOF_Total_min'});

    AssemblyVehicleTable = [AssemblyVehicleTable;newRow];

    %% ADD TO PHASE A BUDGET

    if count == 0
        rows = add(rows,vehicles{j},'Assembly', ...
            'Surface mobility and installation','Baseline',1, ...
            NaN,NaN,'Open','Rango rides Kiki for orbital transfers',true);
        continue;
    end

    if unresolved > 0
        relocationStatus = 'Open';
    else
        relocationStatus = 'Calculated';
    end

    rows = add(rows,vehicles{j},'Assembly', ...
        'Loaded relocation','Baseline',1, ...
        dvOut/count,tofOut/count,relocationStatus, ...
        'assembly_relocations_v4, Kiki fleet',true);
    rows = add(rows,vehicles{j},'Assembly', ...
        'Unloaded return','Baseline',1, ...
        dvBack/count,tofBack/count,relocationStatus, ...
        'assembly_relocations_v4, Kiki fleet',true);
    rows = add(rows,vehicles{j},'Operations', ...
        'Docked loading','Baseline',1, ...
        0,cfg.veicLoading_min,'Assumed','Docked on Atlante',true);
    rows = add(rows,vehicles{j},'Operations', ...
        'Recharging','Baseline',1, ...
        0,cfg.veicCharge_min,'Assumed','Docked on Atlante',true);
    rows = add(rows,vehicles{j},'Operations', ...
        'Free-flight Vbar holding','Baseline',1, ...
        NaN,cfg.veicHold_min,'Open','Perturbed hold unmodelled',true);
    rows = add(rows,vehicles{j},'Assembly', ...
        'Host access transport with Rango','Baseline',1, ...
        NaN,NaN,'Open','Requires host access and transported mass model',true);

end

%% PRINT SELECTED ASSEMBLY RESULTS

fprintf('\n========================================\n');
fprintf('ASSEMBLY RELOCATIONS - SELECTED OBJECTIVE\n');
fprintf('========================================\n');

fprintf('Objective: %s\n',string(opts.objective));
fprintf('Total DV: %.5f m/s\n',A.DV);
fprintf('Total TOF: %.2f min\n',A.TOF_min);

disp(AssemblyVehicleTable);

if A.feasible
    assert(abs(sum(AssemblyVehicleTable.DV_Total_mps)-A.DV) < tol);
    assert(abs(sum(AssemblyVehicleTable.TOF_Total_min)-A.TOF_min) < tol);
end

missionMin = cfg.missionYears * 365.25 * 24 * 60;
validateattributes(cfg.taxiParkingYears,{'numeric'}, ...
    {'scalar','real','finite','nonnegative','<=',cfg.missionYears});
taxiIdleMin = cfg.taxiParkingYears*365.25*24*60;
fprintf('\n========================================\n');
fprintf('MISSION PARKING TIME ESTIMATES\n');
fprintf('========================================\n');
fprintf('Assumed Taxi parking exposure:   %.2f days\n',taxiIdleMin/1440);
fprintf('ATLANTE SK:  %.2f days\n',missionMin/1440);
hSK = cfg.h_target_km - ...
      cfg.launchAltError_km - cfg.hClientSpread_km;
rhoSK = 9.158e-12 * exp(-(hSK-350)/53.298);
fprintf('SK altitude: %.2f km\n',hSK);
fprintf('Atmospheric density: %.5e kg/m^3\n',rhoSK);
vehicleNames = {'Taxi','ATLANTE'};
SK = table();
relativeSK = struct();
atlanteStack = struct('mass_kg',cfg.massAtlante_kg+cfg.massKiki_kg+cfg.massRango_kg, ...
    'area_m2',cfg.areaAtlante_m2,'Cd',cfg.CdAtlante);
host = struct('mass_kg',cfg.massHost_kg,'area_m2',cfg.areaHost_m2,'Cd',cfg.CdHost);
relativeDensity = cfg.relativeSKDensity_kgm3;
if isnan(relativeDensity)
    relativeDensity = 9.158e-12*exp(-(cfg.h_target_km-350)/53.298);
end

relativeOptions = struct('altitude_m',1000*cfg.h_target_km, ...
    'rho_kgm3',relativeDensity,'correctionInterval_hr',cfg.relativeSKInterval_hr, ...
    'AZmargin_m',cfg.relativeSKAZmargin_m);
for j = 1:numel(vehicleNames)
    veh = vehicleNames{j};
    if strcmp(veh,'Taxi')
        mass = cfg.massTaxi_kg;
        area = cfg.areaTaxi_m2;
        Cd = cfg.CdTaxi;
        distance = cfg.taxiDistance_m;
        occupancyMin = taxiIdleMin;
        xEntry = [0;-100;0;0;0;0];
        ellipse = cfg.taxiParkingEllipse;
        departureState = cfg.taxiParkingDepartureState;
        departureWait_s = cfg.taxiParkingWaitBeforeDeparture_s;
        reference = atlanteStack;
        relativeOptions.AZ_m = cfg.taxiParkingAZ_m;
        referenceName = 'Atlante stack';
    else
        mass = cfg.massAtlante_kg + ...
               cfg.massKiki_kg + cfg.massRango_kg;
        area = cfg.areaAtlante_m2;
        Cd = cfg.CdAtlante;
        distance = cfg.atlanteDistance_m;
        occupancyMin = missionMin;
        xEntry = [0;-100;0;0;0;0];
        ellipse = cfg.atlanteParkingEllipse;
        departureState = cfg.atlanteParkingDepartureState;
        departureWait_s = cfg.atlanteParkingWaitBeforeDeparture_s;
        reference = host;
        relativeOptions.AZ_m = cfg.atlanteParkingAZ_m;
        referenceName = 'host structure';
    end

    if isnan(ellipse.ellipseCenter_m)
        ellipse.ellipseCenter_m = distance;
    end

    parking = ellipse;
    parking.mass_kg = mass;
    parking.area_m2 = area;
    parking.Cd = Cd;
    parking.staticDistance_m = distance;
    relAnnual = NaN;
    relRange = nan(2,1);
    azViolations = nan(2,1);
    relativeStatus = 'Open';
    if cfg.runAnnualSK
        relativeSK.(veh) = run_FABER_Taxi_idle_SK(parking,reference,relativeOptions);
        relAnnual = relativeSK.(veh).annualDV_mps;
        relRange = relativeSK.(veh).summary.MinRange_m;
        azViolations = relativeSK.(veh).summary.AZViolationSteps;
        relativeStatus = 'Estimated';
    end

    occupancyYears = occupancyMin/(365.25*24*60);
    DVrelative = relAnnual*occupancyYears;
    rows = add(rows,veh,'Relative stationkeeping', ...
        'Differential-drag compensation','Baseline',1,DVrelative,NaN, ...
        relativeStatus,['Periodic two-impulse HCW relative to ' referenceName],true);
    dragAnnual = station_keeping( ...
        mass,rhoSK,area,Cd,hSK,1);
    dragLifetime = station_keeping( ...
        mass,rhoSK,area,Cd,hSK,occupancyYears);
    for a = 1:2
        if a == 1
            architecture = 'Vbar';

        else
            architecture = 'orbit';

        end

        label = [veh ' ' architecture];

        % Apply absolute drag compensation once for both parking options.
        if a == 1
            rows = add(rows,veh,'Stationkeeping', ...
                'Orbital stationkeeping','Baseline',1, ...
                dragLifetime.DV_SK,0, ...
                'Calculated','station_keeping',true);
        end

        P = parkingAcquisitionHCW( ...
            mu,Re,cfg.h_target_km, ...
            xEntry,distance,architecture, ...
            ellipse,cfg.parkingEntryTime_s);
        rows = add(rows,veh,'Parking', ...
            'Parking acquisition',label,1, ...
            P.DV_mps,P.TOF_min, ...
            'Estimated','Two-impulse HCW transfer',false);

        rows = add(rows,veh,'Parking', ...
            'Hohmann departure burn (already allocated)', ...
            label,1,0,0, ...
            'Included elsewhere', ...
            'Far-field/inter-client orbital transfer',false);
        n = sqrt(mu/(Re+cfg.h_target_km)^3);
        parkingStart = parkingState(n,distance,architecture,ellipse, ...
            cfg.parkingEntryTime_s+departureWait_s);
        D = parkingTransferHCW(n,parkingStart,departureState,cfg.parkingDepartureTime_s);
        rows = add(rows,veh,'Parking', ...
            'Departure state alignment correction', ...
            label,1,D.DV_mps,D.TOF_min, ...
            'Estimated','Two-impulse HCW transfer to pre-burn handoff',false);
        new = table( ...
            string(veh),string(architecture), ...
            occupancyYears, ...
            dragAnnual.DV_SK, ...
            dragLifetime.DV_SK, ...
            relAnnual, ...
            DVrelative, ...
            relRange(a), ...
            azViolations(a), ...
            P.DV_mps, ...
            P.TOF_min, ...
            D.DV_mps,D.TOF_min, ...
            'VariableNames',{ ...
            'Vehicle', ...
            'Architecture', ...
            'ParkingYears', ...
            'DragAnnual_mps', ...
            'DragLifetime_mps', ...
            'RelativeAnnual_mps', ...
            'RelativeLifetime_mps', ...
            'RelativeMinRange_m', ...
            'RelativeAZViolationSteps', ...
            'EntryDV_mps', ...
            'EntryTOF_min', ...
            'DepartureDV_mps','DepartureTOF_min'});
        SK = [SK;new];
    end
end

fprintf('\n========================================\n');
fprintf('PARKING ARCHITECTURE RESULTS\n');
fprintf('========================================\n');
disp(SK);

%% STATIONKEEPING BUDGET OUTPUT

Tsk = cell2table(rows, ...
    'VariableNames',{ ...
    'Vehicle','Segment','Phase','Architecture', ...
    'Count','DV_perEvent_mps','TOF_perEvent_min', ...
    'Status','Basis','IncludeInBaseline', ...
    'DV_lifetime_mps','ActivityTime_min'});

%% ================================================================
% ORBITAL STATIONKEEPING - ANNUAL AND LIFETIME RESULTS
% ================================================================

fprintf('\n========================================\n');
fprintf('ORBITAL STATIONKEEPING\n');
fprintf('========================================\n');

% Extract one row per vehicle
SK_orbital = SK(SK.Architecture == "Vbar", ...
    {'Vehicle','ParkingYears', ...
     'DragAnnual_mps','DragLifetime_mps'});

disp(SK_orbital);

fprintf('\n========================================\n');
fprintf('RELATIVE STATIONKEEPING\n');
fprintf('========================================\n');

SK_relative = SK(:, ...
    {'Vehicle','Architecture','ParkingYears', ...
     'RelativeAnnual_mps','RelativeLifetime_mps'});

disp(SK_relative);

MU = table();
if cfg.includeMomentumUnloading
    rhoMU = rhoSK;
    hMU = hSK;
    CdMU = 2.2;
    CrMU = 1.5;
    kikiActive_min = mean(TOF_outbound_min+TOF_return_min);
    vehicleMU = { ...
        'Taxi'; ...
        'ATLANTE'; ...
        'Kiki'; ...
        'Rango'};
    massMU = [ ...
        cfg.massTaxi_kg; ...
        cfg.massAtlante_kg + ...
            cfg.massKiki_kg + cfg.massRango_kg; ...
        cfg.massKiki_kg; ...
        cfg.massRango_kg];
    durationYearsMU = [ ...
        cfg.missionYears; ...
        cfg.missionYears; ...
        kikiActive_min/(365.25*24*60); ...
        NaN];
    configurationMU = { ...
        'Taxi independent ADCS'; ...
        'ATLANTE with docked payloads'; ...
        'Kiki free flight per relocation round trip'; ...
        'Rango attached operations unresolved'};
    for j = 1:numel(vehicleMU)
        m = massMU(j);
        years = durationYearsMU(j);
        L = (m/100)^(1/3);
        area = L^2;
        I = (1/6)*m*L^2;
        offset = 0.1*L;
        deltaI = 0.1*I;
        leverArm = L/2;
        U = momentum_unloading( ...
            m,rhoMU,hMU,area,CdMU, ...
            offset,CrMU,deltaI, ...
            years,leverArm);
        rows = add(rows,vehicleMU{j}, ...
            'ADCS', ...
            'Momentum unloading', ...
            'Baseline', ...
            1, ...
            U.DV, ...
            NaN, ...
            'Estimated', ...
            configurationMU{j}, ...
            true);
        newRow = table( ...
            string(vehicleMU{j}), ...
            string(configurationMU{j}), ...
            m, ...
            L, ...
            area, ...
            years, ...
            years*365.25*24*60, ...
            U.DV, ...
            'VariableNames',{ ...
                'Vehicle', ...
                'Configuration', ...
                'Mass_kg', ...
                'CharacteristicLength_m', ...
                'ReferenceArea_m2', ...
                'OperatingYears', ...
                'OperatingTime_min', ...
                'BudgetDV_mps'});
        MU = [MU;newRow];
    end
end

fprintf('\n========================================\n');
fprintf('MOMENTUM UNLOADING - PER VEHICLE\n');
fprintf('========================================\n');
disp(MU);
[dvDeorbit,~,tofDeorbit] = hohmann( ...
    mu,r0,Re+cfg.disposalPerigee_km);
for veh = {'Taxi','ATLANTE'}
    rows = add(rows,veh{1},'End of life', ...
        'Deorbit','Baseline',1, ...
        abs(dvDeorbit),tofDeorbit/60, ...
        'Estimated','MA_deltav disposal model',true);
end

rows = add(rows,'Taxi','Contingency', ...
    'Aborted approach','Baseline', ...
    1, ...
    cfg.abortDV_mps,NaN,'Assumed', ...
    'MA_deltav allowance',true);
rows = add(rows,'Taxi','Contingency', ...
    'Collision avoidance','Baseline', ...
    1, ...
    cfg.CAM_DV_mps,NaN,'Assumed', ...
    'MA_deltav allowance',true);
for veh = {'Taxi','ATLANTE'}
    rows = add(rows,veh{1},'Contingency', ...
        'Safety contingency and recovery','Baseline',1,NaN,NaN, ...
        'Open','Fault and recovery allocation unresolved',true);
end

rows = add(rows,'ATLANTE','Contingency','Aborted approach','Baseline', ...
    NaN,NaN,NaN,'Open','Atlante-specific contingency unresolved',true);
rows = add(rows,'ATLANTE','Contingency','Collision avoidance','Baseline', ...
    NaN,NaN,NaN,'Open','Atlante-specific contingency unresolved',true);
T = cell2table(rows,'VariableNames', ...
    {'Vehicle','Segment','Phase','Architecture', ...
    'Count','DV_perEvent_mps','TOF_perEvent_min', ...
    'Status','Basis','IncludeInBaseline', ...
    'DV_lifetime_mps','ActivityTime_min'});
T.Vehicle = string(T.Vehicle);
T.Architecture = string(T.Architecture);
phase = string(T.Phase);
segment = string(T.Segment);
stochastic = ismember(segment,["Outbound RPO","Inbound RPO","Contingency"]);
margin = @(dv,isStochastic) FABER_apply_DV_margin(dv,isStochastic, ...
    cfg.marginFraction,cfg.marginMinimum_mps);
T.MarginedDV_perEvent_mps = margin(T.DV_perEvent_mps,stochastic);
% Apply the floor to the whole event, then allocate proportionally to
% the detailed component rows. Unknown components leave the event open.
for vehicle = ["Taxi","ATLANTE"]
    for group = ["Deployment","Inter-client","Outbound far","Inbound far"]
        idx = T.Vehicle==vehicle & segment==group;
        eventDV = sum(T.DV_perEvent_mps(idx));
        if isnan(eventDV)
            T.MarginedDV_perEvent_mps(idx) = NaN;
        elseif eventDV>0
            assert(all(T.Count(idx)==T.Count(find(idx,1))));
            T.MarginedDV_perEvent_mps(idx) = ...
                T.DV_perEvent_mps(idx)*(margin(eventDV,false)/eventDV);
        end
    end
end

T.MarginedDV_perEvent_mps(T.Vehicle=="Kiki" & phase=="Loaded relocation") = ...
    mean(margin(DV_outbound_mps,false));
T.MarginedDV_perEvent_mps(T.Vehicle=="Kiki" & phase=="Unloaded return") = ...
    mean(margin(DV_return_mps,false));
validateattributes(cfg.interClientRAANReserve_mps,{'numeric'}, ...
    {'scalar','real','finite','nonnegative'});
% Add the reserve after nominal maneuver margins to avoid margin on margin.
T.MarginedDV_perEvent_mps(segment=="Inter-client" & phase=="RAAN reserve") = ...
    cfg.interClientRAANReserve_mps;
T.MarginedDV_mps = T.Count.*T.MarginedDV_perEvent_mps;
T.BudgetUnit = repmat("Per event",height(T),1);
T.BudgetUnit(ismember(segment,["Stationkeeping","Relative stationkeeping"]) | ...
    (segment=="ADCS" & ismember(T.Vehicle,["Taxi","ATLANTE"]))) = "Operating period";
T.BudgetUnit(segment=="ADCS" & T.Vehicle=="Kiki") = "Per relocation round trip";
vehicles = ["Taxi";"Kiki";"Rango";"ATLANTE"];
V = table();
S = table();
for veh = vehicles'
    baselineMask = T.Vehicle==veh & T.IncludeInBaseline;
    V = [V; summarizeBudget(T,baselineMask,veh,"Baseline")];
    if ismember(veh,["Taxi","ATLANTE"])
        for arch = ["Vbar","orbit"]
            mask = baselineMask | (T.Vehicle==veh & T.Architecture==veh+" "+arch);
            S = [S; summarizeBudget(T,mask,veh,arch)];
        end
    end
end

out = struct();
out.cfg = cfg;
out.phaseTable = T;
out.vehicleTotals = V;
out.vehicleScenarioTotals = S;
out.annualSK = SK;
out.relativeSK = relativeSK;
out.momentumUnloading = MU;
out.openItems = T(isnan(T.DV_lifetime_mps) | ...
    isnan(T.ActivityTime_min),:);
out.massAssumptions = table( ...
    ["Taxi";"Kiki";"Rango";"ATLANTE";"Host"], ...
    [cfg.massTaxi_kg;cfg.massKiki_kg; ...
    cfg.massRango_kg;cfg.massAtlante_kg; ...
    cfg.massHost_kg], ...
    'VariableNames',{'Vehicle','Mass_kg'});
out.outboundLambert = outbound;
out.returnLambert = returnTransfer;
out.RPO_MC = RPO_MC;
out.RPO_MC_PhaseTable = RPO_MC_Table;

% Report baseline rows separately from alternative parking architectures.
baseline = T(T.IncludeInBaseline,:);
BudgetTable = table( ...
    baseline.Vehicle,string(baseline.Segment),string(baseline.Phase), ...
    baseline.BudgetUnit,baseline.DV_perEvent_mps,baseline.MarginedDV_perEvent_mps, ...
    baseline.TOF_perEvent_min, ...
    'VariableNames',{'Vehicle','Segment','Maneuver','BudgetUnit', ...
    'RawDV_mps','MarginedDV_mps','ModeledTime_min'});
BudgetTable.Status = string(baseline.Status);
BudgetTable.Basis = string(baseline.Basis);
BudgetTable.Status(isnan(BudgetTable.MarginedDV_mps)) = "Open";
out.budgetTable = BudgetTable;
out.vehicleBudgetTables = struct();
if ~isempty(cfg.vehicleBudgetFolder) && ~isfolder(cfg.vehicleBudgetFolder)
    mkdir(cfg.vehicleBudgetFolder);
end

for j = 1:numel(vehicles)
    name = char(vehicles(j));
    detail = BudgetTable(BudgetTable.Vehicle==name,:);
    out.vehicleBudgetTables.(name) = detail;
    if ~isempty(cfg.vehicleBudgetFolder)
        writetable(detail,fullfile(cfg.vehicleBudgetFolder,[name '.csv']));
    end
end

if ~isempty(cfg.vehicleBudgetFolder)
    writetable(V,fullfile(cfg.vehicleBudgetFolder,'Vehicle_Budget_Summary.csv'));
    writetable(SK,fullfile(cfg.vehicleBudgetFolder,'Parking_Stationkeeping.csv'));
end

out.retrievalAssumptions = struct('moduleAltitude_km',hModule, ...
    'moduleInclination_deg',iModule,'initialRAANOffset_deg',cfg.moduleRAANOffset_deg, ...
    'modeledFarFieldCoastDays',(outbound.transferTime_s/60+returnTOF_min)/1440);

% SAVE ASSEMBLY OPTIMIZATION RESULTS

out.assemblyRelocations = A;

out.assemblyVehicleTable = AssemblyVehicleTable;

out.assemblyTOF_outbound_min = TOF_outbound_min;
out.assemblyTOF_return_min = TOF_return_min;

out.assemblyDV_outbound_mps = DV_outbound_mps;
out.assemblyDV_return_mps = DV_return_mps;

out.assemblyObjective = string(opts.objective);

fprintf('\n====================================================\n');
fprintf('FABER PHASE A - INDIVIDUAL DELTA V BUDGET\n');
fprintf('====================================================\n');
disp(BudgetTable);
fprintf('\n====================================================\n');
fprintf('FABER PHASE A - BUDGET BY UNIT\n');
fprintf('====================================================\n');
disp(V);
fprintf('\n====================================================\n');
fprintf('FABER PHASE A - PARKING SCENARIOS\n');
fprintf('====================================================\n');
disp(S);
end

function summary = summarizeBudget(T,mask,vehicle,architecture)
summary = table();
for unit = unique(T.BudgetUnit(mask),'stable')'
    rows = T(mask & T.BudgetUnit==unit,:);
    summary = [summary; table(vehicle,architecture,unit, ...
        sum(rows.DV_perEvent_mps,'omitnan'),sum(rows.MarginedDV_perEvent_mps,'omitnan'), ...
        nnz(isnan(rows.DV_perEvent_mps) | isnan(rows.MarginedDV_perEvent_mps)), ...
        'VariableNames',{'Vehicle','Architecture','BudgetUnit','KnownRawDV_mps', ...
        'KnownMarginedDV_mps','UnresolvedDV'})];
end
end

function rows = add(rows,veh,seg,phase,arch, ...
    count,dv,tof,status,basis,baseline)
if isfinite(dv)
    lifeDV = count*dv;
else
    lifeDV = NaN;
end

if isfinite(tof)
    lifeTime = count*tof;
else
    lifeTime = NaN;
end

rows(end+1,:) = {veh,seg,phase,arch, ...
    count,dv,tof,status,basis,logical(baseline), ...
    lifeDV,lifeTime};
end

function dv = circularPlaneChange(mu,r,iStart,iEnd,raanOffset)
% Exact angle between orbital normals, expressed as a circular speed change.
sinHalfSquared = sind((iEnd-iStart)/2)^2 + ...
    sind(iStart)*sind(iEnd)*sind(raanOffset/2)^2;
dv = 2000*sqrt(mu/r)*sqrt(max(0,min(1,sinHalfSquared)));
end

function [dv1,dv2,tof] = hohmann(mu,r1,r2)
a = (r1+r2)/2;
dv1 = (sqrt(mu*(2/r1-1/a))-sqrt(mu/r1))*1000;
dv2 = (sqrt(mu/r2)-sqrt(mu*(2/r2-1/a)))*1000;
tof = pi*sqrt(a^3/mu);
end

function P = parkingAcquisitionHCW( ...
    mu,Re,h_km,x0,distance_m,architecture, ...
    ellipse,TOF_s)
n = sqrt(mu/(Re+h_km)^3);
xf = parkingState(n,distance_m,architecture,ellipse,TOF_s);
P = parkingTransferHCW(n,x0,xf,TOF_s);
end

function xf = parkingState(n,distance_m,architecture,ellipse,elapsed_s)
validateattributes(elapsed_s,{'numeric'},{'scalar','real','finite','nonnegative'});
validateattributes(distance_m,{'numeric'},{'scalar','real','finite','positive'});
if strcmpi(architecture,'Vbar')
    xf = [0;-distance_m;0;0;0;0];
elseif strcmpi(architecture,'orbit')
    amp = ellipse.ellipseRadialAmplitude_m;
    theta = n*elapsed_s + ellipse.ellipsePhase_rad;
    xf = [ ...
        amp*cos(theta);
        -ellipse.ellipseCenter_m-2*amp*sin(theta);
        0;
        -amp*n*sin(theta);
        -2*amp*n*cos(theta);
        0];
else
    error('Unknown parking architecture.');
end
end

function P = parkingTransferHCW(n,x0,xf,TOF_s)
validateattributes(x0,{'numeric'},{'vector','numel',6,'real','finite'});
validateattributes(xf,{'numeric'},{'vector','numel',6,'real','finite'});
validateattributes(TOF_s,{'numeric'},{'scalar','real','finite','positive'});
x0 = x0(:);
xf = xf(:);
A = [ ...
    0 0 0 1 0 0;
    0 0 0 0 1 0;
    0 0 0 0 0 1;
    3*n^2 0 0 0 2*n 0;
    0 0 0 -2*n 0 0;
    0 0 -n^2 0 0 0];
Phi = expm(A*TOF_s);
Phi_rr = Phi(1:3,1:3);
Phi_rv = Phi(1:3,4:6);
Phi_vr = Phi(4:6,1:3);
Phi_vv = Phi(4:6,4:6);
r0 = x0(1:3);
v0 = x0(4:6);
rf = xf(1:3);
vf = xf(4:6);
if rcond(Phi_rv) < 1e-10
    error('Ill-conditioned HCW transfer duration.');
end

vTransfer = Phi_rv\(rf-Phi_rr*r0);
dv1 = vTransfer-v0;
vArrival = Phi_vr*r0 + Phi_vv*vTransfer;
dv2 = vf-vArrival;
P.DV_mps = norm(dv1)+norm(dv2);
P.TOF_min = TOF_s/60;
P.dv1 = dv1;
P.dv2 = dv2;
P.targetState = xf;
end
