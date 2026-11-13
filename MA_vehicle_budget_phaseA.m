function out = MA_vehicle_budget_phaseA(cfg)
% FABER Phase A - Vehicle-level Delta-V and TOF budget.
%
% Main outputs:
%   out.phaseTable
%   out.vehicleTotals
%   out.vehicleScenarioTotals
%   out.annualSK
%   out.massAssumptions
%   out.openItems
%
% Requires the existing FABER functions:
%   truss_platform
%   phasing_manoeuvre
%   proximity_operations_v2
%   assembly_relocations_v2
%   momentum_unloading
%   ASV_FormationKeeping_Annual

root = fileparts(mfilename('fullpath'));
addpath(genpath(root));

if nargin < 1 || isempty(cfg)
    cfg = defaultInputs();
else
    defaults = defaultInputs();
    fields = fieldnames(cfg);
    for k = 1:numel(fields)
        defaults.(fields{k}) = cfg.(fields{k});
    end
    cfg = defaults;
end

%% CONSTANTS

if exist('astroConstants','file') == 2
    mu = astroConstants(13);
    Re = astroConstants(23);
else
    mu = 398600.4418;
    Re = 6378.1363;
end

r0 = Re + cfg.h_target_km;

%% STRUCTURE CONFIGURATION

assert(exist('truss_platform','file') ~= 0, ...
    'Missing truss_platform.m');

run(which('truss_platform'));

Nreloc = size(modulePoints,1);

timesMin = [5 10 15 20];
hostPoint = [0;-250;0];

% Preliminary retrieval assumption:
% Five deliveries per year for seven years.
Ndeliveries = cfg.Ndeliveries;

%% RESULTS LEDGER

rows = cell(0,12);

% Vehicle, Segment, Phase, Architecture, Count,
% DV/event, TOF/event, Status, Basis, Included,
% Lifetime DV, Lifetime activity time.

%% INJECTION CORRECTIONS

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

% Taxi injected independently.
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

% ATLANTE carries both Veicolettos.
% Shared maneuver counted once for the stack.
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

%% WORST-CASE INTER-CLIENT ORBIT TRANSFER

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
        cfg.Nclients,DV_altClient,tc/60, ...
        'Estimated','Worst-case client envelope',true);

    rows = add(rows,veh,'Inter-client', ...
        'Extreme plane transfer','Baseline', ...
        cfg.Nclients,DV_plane,0, ...
        'Estimated','Combined plane-angle change',true);
end

%% TAXI FAR-FIELD PHASING

PH = phasing_manoeuvre( ...
    cfg.h_target_km,cfg.h_target_km, ...
    cfg.launchAltError_km+cfg.hClientSpread_km, ...
    cfg.phaseAngle_deg,cfg.phasingRevs,Ndeliveries);

% Outbound
rows = add(rows,'Taxi','Outbound far', ...
    'Altitude correction','Baseline',Ndeliveries, ...
    DV_alt,tinj/60,'Estimated','MA_deltav',true);

rows = add(rows,'Taxi','Outbound far', ...
    'Inclination correction','Baseline',Ndeliveries, ...
    DV_inc,0,'Estimated','MA_deltav',true);

rows = add(rows,'Taxi','Outbound far', ...
    'RAAN correction','Baseline',Ndeliveries, ...
    DV_RAAN,0,'Estimated','MA_deltav',true);

rows = add(rows,'Taxi','Outbound far', ...
    'Phasing','Baseline',Ndeliveries, ...
    PH.DV_phase1_module,PH.time_phase1_min, ...
    'Calculated','phasing_manoeuvre',true);

% Return
[ra,rb,tr] = hohmann(mu, ...
    Re+PH.worst_h_launcher,r0);

DV_returnAlt = abs(ra)+abs(rb);
DV_returnInc = DV_inc;
DV_returnRAAN = DV_RAAN;

rows = add(rows,'Taxi','Inbound far', ...
    'Altitude to ATLANTE','Baseline',Ndeliveries, ...
    DV_returnAlt,tr/60,'Estimated','MA_deltav',true);

rows = add(rows,'Taxi','Inbound far', ...
    'Inclination','Baseline',Ndeliveries, ...
    DV_returnInc,0,'Estimated','MA_deltav',true);

rows = add(rows,'Taxi','Inbound far', ...
    'RAAN','Baseline',Ndeliveries, ...
    DV_returnRAAN,0,'Estimated','MA_deltav',true);

rows = add(rows,'Taxi','Inbound far', ...
    'Phasing to ATLANTE','Baseline',Ndeliveries, ...
    PH.DV_phase2_base,PH.time_phase2_min, ...
    'Calculated','phasing_manoeuvre',true);

% %% TAXI CLOSE-FIELD RPO
% 
% rpo = configureTaxiRPO();
% 
% R1 = proximity_operations_v2( ...
%     mu,Re,cfg.h_target_km,cfg.x0_client,rpo);
% 
% R2 = proximity_operations_v2( ...
%     mu,Re,cfg.h_target_km,cfg.x0_atlante,rpo);
% 
% phases = {'Homing','Hold point 1','Closing', ...
%     'Hold point 2','Final approach'};
% 
% fieldsDV = {'DV_homing','DV_holding1','DV_closing', ...
%     'DV_holding2','DV_final'};
% 
% fieldsTOF = {'TOF_homing','TOF_holding1','TOF_closing', ...
%     'TOF_holding2','TOF_final'};
% 
% for k = 1:numel(phases)
% 
%     rows = add(rows,'Taxi','Outbound RPO', ...
%         phases{k},'Baseline',Ndeliveries, ...
%         R1.(fieldsDV{k}),R1.(fieldsTOF{k}), ...
%         'Calculated','proximity_operations_v2',true);
% 
%     rows = add(rows,'Taxi','Inbound RPO', ...
%         phases{k},'Baseline',Ndeliveries, ...
%         R2.(fieldsDV{k}),R2.(fieldsTOF{k}), ...
%         'Estimated','ATLANTE arrival assumed',true);
% end


%% ================================================================
% TAXI CLOSE-FIELD RPO - TWO MONTE CARLO RUNS
% ================================================================
%
% Run 1: Taxi -> Client
% Run 2: Taxi -> ATLANTE
%
% For each approach:
% - Monte Carlo samples feasible trajectories
% - Select the trajectory at the upper empirical P95 DV
% - Extract DV and TOF for its five phases
%
% The five phase values belong to ONE sampled trajectory.
% Therefore their sums equal the selected leg's DV and TOF.
%
% P95 here represents the sampled design-space distribution,
% NOT mission reliability.

mcCfg = ASV_RPO_MC_config();

mcCfg.mu_km3s2 = mu;
mcCfg.Re_km = Re;
mcCfg.h_km = cfg.h_target_km;

mcCfg.x0_client = cfg.x0_client(:);
mcCfg.x0_host = cfg.x0_atlante(:);

mcCfg.N = cfg.N_RPO_MC;
mcCfg.makePlots = false;
mcCfg.saveCSV = false;

%% RUN 1 - OUTBOUND TO CLIENT

fprintf('\n========================================\n');
fprintf('MONTE CARLO 1 - TAXI TO CLIENT\n');
fprintf('========================================\n');

mcCfg.seed = cfg.RPO_seed_client;

MC_client = ASV_RPO_MonteCarlo(mcCfg);

% MC run 1: outbound approach to client
DV_client_trials = MC_client.DV_phase_mps(:,1:7);
TOF_client_trials = MC_client.TOF_phase_min(:,1:7);

% Total DV and TOF of every accepted client approach
DV_client_total = sum(DV_client_trials,2);
TOF_client_total = sum(TOF_client_trials,2);

% Select the actual sampled trajectory at empirical P95 DV.
[DV_client_sorted,idxClientSorted] = ...
    sort(DV_client_total,'ascend');

idx95_client = max(1, ...
    ceil(0.95*numel(DV_client_sorted)));

idxClient = idxClientSorted(idx95_client);

% Phase results for the selected trajectory
DV_client_phase = DV_client_trials(idxClient,:);
TOF_client_phase = TOF_client_trials(idxClient,:);

% Independent P95 total TOF (diagnostic only)
P95_TOF_client = prctile(TOF_client_total,95);

%% RUN 2 - INBOUND TO ATLANTE

fprintf('\n========================================\n');
fprintf('MONTE CARLO 2 - TAXI TO ATLANTE\n');
fprintf('========================================\n');

mcCfg.seed = cfg.RPO_seed_atlante;

MC_atlante = ASV_RPO_MonteCarlo(mcCfg);

%  MC run 2: inbound approach to ATLANTE
DV_atlante_trials = MC_atlante.DV_phase_mps(:,8:14);
TOF_atlante_trials = MC_atlante.TOF_phase_min(:,8:14);

DV_atlante_total = sum(DV_atlante_trials,2);
TOF_atlante_total = sum(TOF_atlante_trials,2);

% Select the actual sampled trajectory at empirical P95 DV.
[DV_atlante_sorted,idxAtlanteSorted] = ...
    sort(DV_atlante_total,'ascend');

idx95_atlante = max(1, ...
    ceil(0.95*numel(DV_atlante_sorted)));

idxAtlante = idxAtlanteSorted(idx95_atlante);

DV_atlante_phase = DV_atlante_trials(idxAtlante,:);
TOF_atlante_phase = TOF_atlante_trials(idxAtlante,:);

P95_TOF_atlante = prctile(TOF_atlante_total,95);

%% FIVE RPO PHASES

phases = { ...
    'Homing 1', ...
    'Hold point 1', ...
    'Homing 2', ...
    'Hold point 2', ...
    'Closing', ...
    'Hold point 3', ...
    'Final approach'};

%% ADD OUTBOUND AND INBOUND TO MISSION BUDGET

for k = 1:numel(phases)

    rows = add(rows,'Taxi','Outbound RPO', ...
        phases{k},'Baseline',Ndeliveries, ...
        DV_client_phase(k), ...
        TOF_client_phase(k), ...
        'MC P95-DV design', ...
        'MC run 1: selected client approach',true);

    rows = add(rows,'Taxi','Inbound RPO', ...
        phases{k},'Baseline',Ndeliveries, ...
        DV_atlante_phase(k), ...
        TOF_atlante_phase(k), ...
        'MC P95-DV design', ...
        'MC run 2: selected ATLANTE approach',true);
end

%% SUMMARY

DV_RPO_client = sum(DV_client_phase);
DV_RPO_atlante = sum(DV_atlante_phase);

TOF_RPO_client = sum(TOF_client_phase);
TOF_RPO_atlante = sum(TOF_atlante_phase);

DV_RPO_cycle = DV_RPO_client + DV_RPO_atlante;
TOF_RPO_cycle = TOF_RPO_client + TOF_RPO_atlante;

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

fprintf('\nOutbound to client:\n');
fprintf('  Selected DV:  %.4f m/s\n',DV_RPO_client);
fprintf('  Selected TOF: %.2f min\n',TOF_RPO_client);
fprintf('  Independent P95 TOF: %.2f min\n', ...
    P95_TOF_client);

fprintf('\nInbound to ATLANTE:\n');
fprintf('  Selected DV:  %.4f m/s\n',DV_RPO_atlante);
fprintf('  Selected TOF: %.2f min\n',TOF_RPO_atlante);
fprintf('  Independent P95 TOF: %.2f min\n', ...
    P95_TOF_atlante);

fprintf('\nComplete retrieval RPO:\n');
fprintf('  Selected DV:  %.4f m/s\n',DV_RPO_cycle);
fprintf('  Selected TOF: %.2f min\n',TOF_RPO_cycle);

fprintf('\nLifetime RPO (%d retrievals):\n',Ndeliveries);
fprintf('  DV: %.4f m/s\n', ...
    Ndeliveries*DV_RPO_cycle);

%% STORE MC RESULTS

RPO_MC = struct();

RPO_MC.Client = MC_client;
RPO_MC.ATLANTE = MC_atlante;

RPO_MC.SelectedClientTrial = idxClient;
RPO_MC.SelectedATLANTETrial = idxAtlante;

RPO_MC.PhaseTable = RPO_MC_Table;

RPO_MC.CycleDV_mps = DV_RPO_cycle;
RPO_MC.CycleTOF_min = TOF_RPO_cycle;

RPO_MC.P95ClientTOF_min = P95_TOF_client;
RPO_MC.P95ATLANTETOF_min = P95_TOF_atlante;


%% CAPTURE AND HANDOVER

rows = add(rows,'Taxi','Payload', ...
    'Capture module','Baseline',Ndeliveries, ...
    NaN,cfg.captureWait_min, ...
    'Open','Capture dynamics not modelled',true);

rows = add(rows,'Taxi','Payload', ...
    'ATLANTE handover','Baseline',Ndeliveries, ...
    NaN,cfg.handoverWait_min, ...
    'Open','Handover dynamics not modelled',true);

%% VEICOLETTO ASSEMBLY RELOCATIONS

A = assembly_relocations_v2( ...
    mu,Re,cfg.h_target_km,101, ...
    timesMin,hostPoint,modulePoints);

vehicles = {'Kiki','Rango'};

for j = 1:2

    indices = j:2:Nreloc;
    count = numel(indices);

    dvOut = sum(A.DV_outbound(indices));
    dvBack = sum(A.DV_return(indices));

    tofOut = sum(A.timePerLeg_min(indices));
    tofBack = tofOut;

    rows = add(rows,vehicles{j},'Assembly', ...
        'Loaded relocation','Baseline',1, ...
        dvOut,tofOut,'Calculated', ...
        'assembly_relocations_v2',true);

    rows = add(rows,vehicles{j},'Assembly', ...
        'Unloaded return','Baseline',1, ...
        dvBack,tofBack,'Calculated', ...
        'assembly_relocations_v2',true);

    rows = add(rows,vehicles{j},'Operations', ...
        'Docked loading','Baseline',count, ...
        0,cfg.veicLoading_min,'Assumed', ...
        'No independent translation while docked',true);

    rows = add(rows,vehicles{j},'Operations', ...
        'Recharging','Baseline',count, ...
        0,cfg.veicCharge_min,'Assumed', ...
        'Docked charging',true);

    rows = add(rows,vehicles{j},'Operations', ...
        'Free-flight Vbar holding','Baseline',count, ...
        NaN,cfg.veicHold_min,'Open', ...
        'Perturbed control not modelled',true);
end


%% ================================================================
% TAXI AND ATLANTE STATION KEEPING - PHASE A
% ================================================================
%
% Components:
% 1. Absolute orbital drag compensation
% 2. HCW parking acquisition
% 3. Relative maintenance allocation (literature-informed)
% 4. Parking departure: Hohmann burn already budgeted elsewhere
%
% Parking architectures:
% - Vbar: stationary negative V-bar holding point
% - orbit: bounded HCW relative ellipse
%
% Taxi idle: seven-year duration minus modelled retrieval time
% ATLANTE: seven-year operating lifetime
%
% The two architectures are ALTERNATIVES, not summed.
% No annual numerical feedback simulation is required.

missionMin = cfg.missionYears * 365.25 * 24 * 60;

%% TAXI IDLE DURATION

taxiBusyMin = 0;

for k = 1:size(rows,1)

    isTaxi = strcmp(rows{k,1},'Taxi');
    isBaseline = strcmp(rows{k,4},'Baseline');

    % Only count known operational activity.
    % Do not subtract injection, inter-client or deorbit
    % unless those are separately accounted for in schedule.
    isActivity = ismember(rows{k,2}, ...
        {'Outbound far','Inbound far', ...
         'Outbound RPO','Inbound RPO','Payload'});

    if isTaxi && isBaseline && isActivity

        t = rows{k,12};

        if isnumeric(t) && isfinite(t)
            taxiBusyMin = taxiBusyMin + t;
        end
    end
end

taxiIdleMin = max(0,missionMin-taxiBusyMin);

fprintf('\n========================================\n');
fprintf('MISSION PARKING TIME ESTIMATES\n');
fprintf('========================================\n');
fprintf('Taxi idle:   %.2f days\n',taxiIdleMin/1440);
fprintf('ATLANTE SK:  %.2f days\n',missionMin/1440);

%% ATMOSPHERIC ASSUMPTIONS FROM MA_deltav

% Conservative altitude selection inherited from MA_deltav.
% Must be treated as a scenario, not a guaranteed worst case.

hSK = cfg.h_target_km - ...
      cfg.launchAltError_km - cfg.hClientSpread_km;

rhoSK = 9.158e-12 * exp(-(hSK-350)/53.298);

fprintf('SK altitude: %.2f km\n',hSK);
fprintf('Atmospheric density: %.5e kg/m^3\n',rhoSK);

%% PARKING CONFIGURATIONS

vehicleNames = {'Taxi','ATLANTE'};

SK = table();

for j = 1:numel(vehicleNames)

    veh = vehicleNames{j};

    if strcmp(veh,'Taxi')

        mass = cfg.massTaxi_kg;
        area = cfg.areaTaxi_m2;
        Cd = cfg.CdTaxi;

        distance = cfg.taxiDistance_m;

        occupancyMin = taxiIdleMin;

        % Nominal free-flight state after return/handover.
        % Replace with actual post-handover state later.
        xEntry = [0;-100;0;0;0;0];

        % Preliminary: one parking entry after each retrieval.
        Nentries = Ndeliveries;

    else

        % ATLANTE with Kiki and Rango docked.
        % Representative configuration, not constant real mass.
        mass = cfg.massAtlante_kg + ...
               cfg.massKiki_kg + cfg.massRango_kg;

        area = cfg.areaAtlante_m2;
        Cd = cfg.CdAtlante;

        distance = cfg.atlanteDistance_m;

        occupancyMin = missionMin;

        % Representative state before acquiring host-relative parking.
        xEntry = [0;-100;0;0;0;0];

        % Initial acquisition plus subsequent client transitions.
        % Replace with exact mission sequence when known.
        Nentries = cfg.Nclients + 1;
    end

    %% ABSOLUTE ORBITAL DRAG COMPENSATION

    occupancyYears = occupancyMin/(365.25*24*60);

    dragAnnual = station_keeping( ...
        mass,rhoSK,area,Cd,hSK,1);

    dragLifetime = station_keeping( ...
        mass,rhoSK,area,Cd,hSK,occupancyYears);

    %% ALTERNATIVE ARCHITECTURES

    for a = 1:2

        if a == 1
            architecture = 'Vbar';
            relAnnual = cfg.relativeSK_Vbar_mpsYr;
        else
            architecture = 'orbit';
            relAnnual = cfg.relativeSK_orbit_mpsYr;
        end

        label = [veh ' ' architecture];

        %% 1. ORBITAL DRAG COMPENSATION

        rows = add(rows,veh,'Parking', ...
            'Orbital drag compensation',label,1, ...
            dragLifetime.DV_SK,occupancyMin, ...
            'Estimated','MA_deltav drag model',false);

        %% 2. PARKING ACQUISITION (HCW)

        P = parkingAcquisitionHCW( ...
            mu,Re,cfg.h_target_km, ...
            xEntry,distance,architecture, ...
            cfg.orbitRadialAmplitude_m, ...
            cfg.parkingEntryTime_s);

        rows = add(rows,veh,'Parking', ...
            'Parking acquisition',label,Nentries, ...
            P.DV_mps,P.TOF_min, ...
            'Estimated','Two-impulse HCW transfer',false);

        %% 3. RELATIVE FORMATION MAINTENANCE

        % Preliminary allowance independent of absolute
        % orbital drag compensation.
        %
        % This represents differential disturbances,
        % navigation/control corrections and geometry keeping.
        % Validate later using perturbed relative dynamics.

        DVrelative = relAnnual * occupancyYears;

        rows = add(rows,veh,'Parking', ...
            'Relative maintenance allocation',label,1, ...
            DVrelative,0, ...
            'Assumed','Formation-flying trade allowance',false);

        %% 4. PARKING DEPARTURE

        % The nominal departure Hohmann burn is already
        % counted under the transfer budget.
        %
        % Do NOT add the first Hohmann burn again.
        %
        % An additional relative-state alignment correction
        % may still be required before initiating the burn.

        rows = add(rows,veh,'Parking', ...
            'Hohmann departure burn (already allocated)', ...
            label,Nentries,0,0, ...
            'Included elsewhere', ...
            'Far-field/inter-client orbital transfer',false);

        rows = add(rows,veh,'Parking', ...
            'Departure state alignment correction', ...
            label,Nentries,NaN,NaN, ...
            'Open', ...
            'Requires departure state and burn geometry',false);

        %% RESULTS TABLE

        new = table( ...
            string(veh),string(architecture), ...
            occupancyYears, ...
            dragAnnual.DV_SK, ...
            dragLifetime.DV_SK, ...
            relAnnual, ...
            DVrelative, ...
            P.DV_mps, ...
            P.TOF_min, ...
            Nentries, ...
            Nentries*P.DV_mps, ...
            'VariableNames',{ ...
            'Vehicle', ...
            'Architecture', ...
            'ParkingYears', ...
            'DragAnnual_mps', ...
            'DragLifetime_mps', ...
            'RelativeAnnual_mps', ...
            'RelativeLifetime_mps', ...
            'EntryDV_mps', ...
            'EntryTOF_min', ...
            'Nentries', ...
            'EntryLifetime_mps'});

        SK = [SK;new];

    end
end

fprintf('\n========================================\n');
fprintf('PARKING ARCHITECTURE RESULTS\n');
fprintf('========================================\n');

disp(SK);


%% ================================================================
% MOMENTUM UNLOADING - PER VEHICLE
% ================================================================

% Phase A assumptions:
% - Taxi: independent ADCS throughout mission
% - ATLANTE: independent ADCS throughout mission
% - Kiki/Rango: independent ADCS during free flight only
% - Docked Veicolettos are managed by ATLANTE ADCS
% - Cubic spacecraft approximation for geometry
% - Disturbance environment held constant at hSK
%
% All Delta-V values are preliminary estimates.

MU = table();

if cfg.includeMomentumUnloading

    % Effective atmospheric density and altitude
    rhoMU = rhoSK;
    hMU = hSK;

    % Geometry and aerodynamic assumptions
    CdMU = 2.2;
    CrMU = 1.5;

    % Kiki and Rango active durations [min]
    % Each relocation includes outbound and return.

    kikiIndices = 1:2:Nreloc;
    rangoIndices = 2:2:Nreloc;

    kikiActive_min = ...
        2*sum(A.timePerLeg_min(kikiIndices));

    rangoActive_min = ...
        2*sum(A.timePerLeg_min(rangoIndices));

    % Vehicle configuration
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
        rangoActive_min/(365.25*24*60)];

    configurationMU = { ...
        'Taxi independent ADCS'; ...
        'ATLANTE with docked payloads'; ...
        'Kiki free flight'; ...
        'Rango free flight'};

    % Individual calculations

    for j = 1:numel(vehicleMU)

        m = massMU(j);
        years = durationYearsMU(j);

        % Cubic spacecraft approximation
        L = (m/100)^(1/3);
        area = L^2;

        % Approximate inertia
        I = (1/6)*m*L^2;

        % Disturbance and unloading assumptions
        offset = 0.1*L;
        deltaI = 0.1*I;
        leverArm = L/2;

        U = momentum_unloading( ...
            m,rhoMU,hMU,area,CdMU, ...
            offset,CrMU,deltaI, ...
            years,leverArm);

        % Add to mission Delta-V ledger
        % Count = 1 because U.DV already covers
        % the entire allocated operating duration.

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

        % Detailed results for Excel

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
                'LifetimeDV_mps'});

        MU = [MU;newRow];

    end

end

fprintf('\n========================================\n');
fprintf('MOMENTUM UNLOADING - PER VEHICLE\n');
fprintf('========================================\n');

disp(MU);

%% DEORBIT

[dvDeorbit,~,tofDeorbit] = hohmann( ...
    mu,r0,Re+cfg.disposalPerigee_km);

for veh = {'Taxi','ATLANTE'}

    rows = add(rows,veh{1},'End of life', ...
        'Deorbit','Baseline',1, ...
        abs(dvDeorbit),tofDeorbit/60, ...
        'Estimated','MA_deltav disposal model',true);
end

%% CONTINGENCY

rows = add(rows,'Taxi','Contingency', ...
    'Aborted approach','Baseline', ...
    cfg.abortFraction*Ndeliveries, ...
    cfg.abortDV_mps,NaN,'Assumed', ...
    'MA_deltav allowance',true);

rows = add(rows,'Taxi','Contingency', ...
    'Collision avoidance','Baseline', ...
    cfg.CAMperYear*cfg.missionYears, ...
    cfg.CAM_DV_mps,NaN,'Assumed', ...
    'MA_deltav allowance',true);

%% CREATE RESULTS TABLES

T = cell2table(rows,'VariableNames', ...
    {'Vehicle','Segment','Phase','Architecture', ...
    'Count','DV_perEvent_mps','TOF_perEvent_min', ...
    'Status','Basis','IncludeInBaseline', ...
    'DV_lifetime_mps','ActivityTime_min'});

T.Vehicle = string(T.Vehicle);
T.Architecture = string(T.Architecture);

T.MarginedDV_mps = T.DV_lifetime_mps * ...
    (1+cfg.marginFraction);

vehicles = ["Taxi";"Kiki";"Rango";"ATLANTE"];

V = table();

for j = 1:numel(vehicles)

    mask = T.Vehicle==vehicles(j) & T.IncludeInBaseline;

    sub = T(mask,:);

    new = table(vehicles(j), ...
        sum(sub.DV_lifetime_mps,'omitnan'), ...
        sum(sub.MarginedDV_mps,'omitnan'), ...
        sum(sub.ActivityTime_min,'omitnan'), ...
        sum(isnan(sub.DV_lifetime_mps)), ...
        'VariableNames', ...
        {'Vehicle','KnownDV_mps', ...
        'KnownMarginedDV_mps','ActivityTime_min', ...
        'UnresolvedDV'});

    V = [V;new];
end

S = table();

for j = [1 4]

    veh = vehicles(j);

    for arch = ["Vbar","orbit"]

        mask = T.Vehicle==veh & ...
            (T.IncludeInBaseline | ...
            T.Architecture==veh+" "+arch);

        sub = T(mask,:);

        new = table(veh,arch, ...
            sum(sub.DV_lifetime_mps,'omitnan'), ...
            sum(sub.MarginedDV_mps,'omitnan'), ...
            sum(isnan(sub.DV_lifetime_mps)), ...
            'VariableNames', ...
            {'Vehicle','Architecture','KnownDV_mps', ...
            'KnownMarginedDV_mps','UnresolvedDV'});

        S = [S;new];
    end
end

%% OUTPUT

out = struct();

out.cfg = cfg;
out.phaseTable = T;
out.vehicleTotals = V;
out.vehicleScenarioTotals = S;
out.annualSK = SK;
out.momentumUnloading = MU;
out.openItems = T(isnan(T.DV_lifetime_mps) | ...
    isnan(T.ActivityTime_min),:);

out.massAssumptions = table( ...
    ["Taxi";"Kiki";"Rango";"ATLANTE";"Host"], ...
    [cfg.massTaxi_kg;cfg.massKiki_kg; ...
    cfg.massRango_kg;cfg.massAtlante_kg; ...
    cfg.massHost_kg], ...
    'VariableNames',{'Vehicle','Mass_kg'});

out.RPO_MC = RPO_MC;
out.RPO_MC_PhaseTable = RPO_MC_Table;

disp(V);
disp(S);

if cfg.saveCSV

    if ~exist(cfg.outputFolder,'dir')
        mkdir(cfg.outputFolder);
    end

    writetable(T,fullfile(cfg.outputFolder, ...
        'FABER_phase_budget.csv'));

    writetable(V,fullfile(cfg.outputFolder, ...
        'FABER_vehicle_totals.csv'));

    writetable(S,fullfile(cfg.outputFolder, ...
        'FABER_parking_scenarios.csv'));

    writetable(out.massAssumptions, ...
        fullfile(cfg.outputFolder, ...
        'FABER_masses.csv'));
end

end

%% LOCAL FUNCTIONS

function cfg = defaultInputs()

cfg = struct();

cfg.h_target_km = 457.25;
cfg.i_target_deg = 51.02;

cfg.missionYears = 7;
cfg.Nclients = 3;
cfg.Ndeliveries = 35;

cfg.launchAltError_km = 15;
cfg.launchIncError_deg = 0.1;
cfg.launchRAANError_deg = 0.1;

cfg.hClientSpread_km = 67.47;
cfg.iClientSpread_deg = 2.81;
cfg.clientRAANSpread_deg = 0;

cfg.phaseAngle_deg = 20;
cfg.phasingRevs = 35;
cfg.disposalPerigee_km = 250;

cfg.x0_client = [50;-2000;20;0.01;0.04;-0.01];
cfg.x0_atlante = cfg.x0_client;

cfg.abortFraction = 0.5;
cfg.abortDV_mps = 1;
cfg.CAMperYear = 4;
cfg.CAM_DV_mps = 1;

cfg.marginFraction = 0.05;

% Vehicle mass placeholders [kg]
cfg.massTaxi_kg = 3503.59;
cfg.massKiki_kg = 1813.62;
cfg.massRango_kg = 1813.62;
cfg.massAtlante_kg = 4805.47;  
cfg.massHost_kg = 1000;    

% Area assumptions [m^2]
cfg.areaTaxi_m2 = (cfg.massTaxi_kg/100)^(2/3);
% Phase A: consistent equivalent-cube geometry
% ATLANTE + two docked Veicolettos
massAtlanteStack = cfg.massAtlante_kg + ...
                   cfg.massKiki_kg + cfg.massRango_kg;

cfg.areaAtlante_m2 = (massAtlanteStack/100)^(2/3);
cfg.areaHost_m2 = 8;

cfg.CdTaxi = 2.2;
cfg.CdAtlante = 2.2;
cfg.CdHost = 2.2;

% Parking
cfg.taxiDistance_m = 1000;
cfg.taxiDeadband_m = 25;

cfg.atlanteDistance_m = 250;
cfg.atlanteDeadband_m = 10;

cfg.orbitRadialAmplitude_m = 20;

cfg.skStep_s = 60;
cfg.rhoRef_kgm3 = 1e-12;
cfg.scaleHeight_m = 60000;

% Operations
cfg.captureWait_min = NaN;
cfg.handoverWait_min = NaN;
cfg.veicLoading_min = NaN;
cfg.veicCharge_min = NaN;
cfg.veicHold_min = NaN;

cfg.runAnnualSK = true;
cfg.includeMomentumUnloading = true;


% PARKING SETTINGS

% Nominal parking distances [m]
cfg.taxiDistance_m = 1000;
cfg.atlanteDistance_m = 250;

% HCW relative orbit radial amplitude [m]
% Along-track amplitude is twice this value.
cfg.orbitRadialAmplitude_m = 20;

% Time to acquire parking point/orbit [s]
cfg.parkingEntryTime_s = 3600;

% RELATIVE STATION-KEEPING ALLOWANCES
% [m/s/year]
%
% Literature-informed Phase A sensitivity assumptions.
% These are NOT validated FABER requirements.
%
% Suggested sensitivity:
% Low:          1 m/s/year
% Reference:   10 m/s/year
% Conservative:40 m/s/year
%
% V-bar and orbit values may be independently varied.

cfg.relativeSK_Vbar_mpsYr = 0;
cfg.relativeSK_orbit_mpsYr = 0;


%% TAXI RPO MONTE CARLO SETTINGS

% Number of accepted samples per Monte Carlo run
cfg.N_RPO_MC = 1000;

% Independent random seeds
cfg.RPO_seed_client = 42;
cfg.RPO_seed_atlante = 2026;



cfg.saveCSV = true;
cfg.outputFolder = 'FABER_vehicle_budget_results';

end

function rpo = configureTaxiRPO()

rpo = struct();
rpo.N = 600;

rpo.homing.targets = [0;-1000;0;0;0;0];
rpo.homing.times_s = 3600;

rpo.holding1.mode = 'point';
rpo.holding1.param.yHold = -1000;
rpo.holding1.duration_s = 30*60;

rpo.closing.targets = [0;-100;0;0;0;0];
rpo.closing.methods = {'Vbar'};
rpo.closing.times_s = 1800;

rpo.holding2.mode = 'point';
rpo.holding2.param.yHold = -100;
rpo.holding2.duration_s = 10*60;

rpo.final.target = [0;-5;0;0;0;0];
rpo.final.method = 'Vbar';
rpo.final.duration_s = 600;

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

function [dv1,dv2,tof] = hohmann(mu,r1,r2)

a = (r1+r2)/2;

dv1 = (sqrt(mu*(2/r1-1/a))-sqrt(mu/r1))*1000;
dv2 = (sqrt(mu/r2)-sqrt(mu*(2/r2-1/a)))*1000;

tof = pi*sqrt(a^3/mu);

end


function P = parkingAcquisitionHCW( ...
    mu,Re,h_km,x0,distance_m,architecture, ...
    radialAmp_m,TOF_s)
% PARKINGACQUISITIONHCW
%
% Calculates a two-impulse transfer from an initial
% host-relative state to:
%
% 1. Stationary negative V-bar point
% 2. Bounded in-plane HCW relative orbit
%
% Inputs:
% mu            [km^3/s^2]
% Re            [km]
% h_km          [km]
% x0            [x;y;z;vx;vy;vz] SI units
% distance_m    [m]
% radialAmp_m   [m]
% TOF_s         [s]
%
% Outputs:
% P.DV_mps
% P.TOF_min
% P.dv1
% P.dv2
% P.targetState
%
% Important:
% - Ideal HCW model only.
% - No finite-thrust constraints.
% - No KOZ/passive safety enforcement.
% - No parking departure maneuver included.

n = sqrt(mu/(Re+h_km)^3);

A = [ ...
    0 0 0 1 0 0;
    0 0 0 0 1 0;
    0 0 0 0 0 1;
    3*n^2 0 0 0 2*n 0;
    0 0 0 -2*n 0 0;
    0 0 -n^2 0 0 0];

x0 = x0(:);

assert(numel(x0)==6, ...
    'Initial state must contain six components.');

assert(TOF_s > 0);
assert(distance_m > 0);

if strcmpi(architecture,'Vbar')

    xf = [0;-distance_m;0;0;0;0];

elseif strcmpi(architecture,'orbit')

    % Target point on a natural in-plane HCW ellipse.
    %
    % The ellipse is centred at -Vbar distance.
    % It has radial amplitude A and along-track
    % amplitude 2A.
    %
    % The phase at the arrival time is selectable.
    % Here zero phase is assumed at the start.

    amp = radialAmp_m;

    theta = n*TOF_s;

    xf = [ ...
        amp*cos(theta);
        -distance_m-2*amp*sin(theta);
        0;
        -amp*n*sin(theta);
        -2*amp*n*cos(theta);
        0];

else
    error('Unknown parking architecture.');
end

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

% First impulse
vTransfer = Phi_rv\(rf-Phi_rr*r0);

dv1 = vTransfer-v0;

% Second impulse to match the parking velocity
vArrival = Phi_vr*r0 + Phi_vv*vTransfer;

dv2 = vf-vArrival;

P.DV_mps = norm(dv1)+norm(dv2);

P.TOF_min = TOF_s/60;

P.dv1 = dv1;
P.dv2 = dv2;

P.targetState = xf;

end

