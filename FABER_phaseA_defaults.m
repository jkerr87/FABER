function cfg = FABER_phaseA_defaults()
% Return shared assumptions for the Phase A budget and round-trip study.
cfg = struct();

% cfg.h_target_km = 680.12;
% cfg.i_target_deg = 98.12;
% cfg.hClientSpread_km = 72.84;
% cfg.iClientSpread_deg = 0.31;

cfg.h_target_km = 457.25;
cfg.i_target_deg = 51.02;
cfg.hClientSpread_km = 67.47;
cfg.iClientSpread_deg = 2.81;

cfg.missionYears = 7;
cfg.taxiParkingYears = cfg.missionYears; % Continuous-parking exposure until a duty cycle is supplied.

% Deployment dispersions and the inter-client planning reserve.
cfg.launchAltError_km = 15;
cfg.launchIncError_deg = 0.1;
cfg.launchRAANError_deg = 0.1;
cfg.clientRAANSpread_deg = 0;
cfg.interClientRAANReserve_mps = 100;

% Nominal retrieval geometry and Lambert search settings.
cfg.phaseAngle_deg = 20;
cfg.phasingRevs = 10;
cfg.outboundLambertTOF_s = []; % Empty searches around the nominal revolution count.
cfg.outboundLambertRevolutions = [];
cfg.returnLambertTOF_s = [];
cfg.returnLambertRevolutions = [];

cfg.farFieldIsp_s = 315.5;
% Installed hardware from the vehicle propulsion table. The firing-time
% estimates assume one main engine or four co-aligned RCS thrusters.
cfg.atlanteTransferPropulsion = struct('type',"MAIN - Nammo LEROS 4-ET", ...
    'nInstalled',5,'nActiveThrusters',1,'thrust_N',1310,'Isp_s',321.0,'thrustEfficiency',1);
cfg.atlanteParkingPropulsion = struct('type',"RCS - IHI HBT-1", ...
    'nInstalled',20,'nActiveThrusters',4,'thrust_N',124.4,'Isp_s',274.9,'thrustEfficiency',1);
cfg.kikiTransferPropulsion = struct('type',"RCS/translation - Moog DST-12", ...
    'nInstalled',16,'nActiveThrusters',4,'thrust_N',22,'Isp_s',297.0,'thrustEfficiency',1);
cfg.outboundMinimumPerigee_km = 160;
cfg.outboundDepartureEpoch_s = 0; % Seconds on the scenario clock.
cfg.outboundDepartureArgumentOfLatitude_deg = 0;
cfg.outboundAtlanteRAAN_deg = 0;

cfg.moduleAltitudeOffset_km = -15;
cfg.moduleInclinationOffset_deg = 0.1;
cfg.moduleRAANOffset_deg = 0.1;
cfg.payloadMass_kg = 2000;

cfg.assemblyTransferTimes_min = 5:5:130;
cfg.vehicleBudgetFolder = '';
cfg.vehicleReportTemplate = "";
cfg.vehicleReportWorkbook = "";
% PassiveDeorbit reference assumptions for end-of-life configurations.
cfg.disposalPerigee_km = 160; % Drag integration cutoff.
cfg.disposalDensityProfile = "Low";
cfg.disposalTaxiMass_kg = 1059.88;
cfg.disposalTaxiArea_m2 = 15;
cfg.disposalAtlanteMass_kg = 4805.47+2*606.58+2*230; % Attached Kiki/Rango pairs.
cfg.disposalAtlanteArea_m2 = 100;

% Target-relative gate states in LVLH meters and meters per second.
cfg.x0_client = [50;-2000;20;0.01;0.04;-0.01];
cfg.x0_atlante = cfg.x0_client;

cfg.abortDV_mps = 1;
cfg.CAM_DV_mps = 1;
cfg.marginFraction = 0.05;
cfg.marginMinimum_mps = 10;

% Vehicle properties used in drag and propulsion estimates.
cfg.massTaxi_kg = 3503.59;
cfg.massKiki_kg = 1813.62;
cfg.massRango_kg = 1813.62;
cfg.massAtlante_kg = 4805.47;
cfg.massHost_kg = 1000;
cfg.areaTaxi_m2 = (cfg.massTaxi_kg/100)^(2/3);
massAtlanteStack = cfg.massAtlante_kg + ...
                   cfg.massKiki_kg + cfg.massRango_kg;
cfg.areaAtlante_m2 = (massAtlanteStack/100)^(2/3);
cfg.areaHost_m2 = 8;
cfg.CdTaxi = 2.2;
cfg.CdAtlante = 2.2;
cfg.CdHost = 2.2;

cfg.taxiDistance_m = 1000;
cfg.taxiDeadband_m = 25;
cfg.atlanteDistance_m = 250;
cfg.atlanteDeadband_m = 10;
cfg.skStep_s = 60;

cfg.rhoRef_kgm3 = 1e-12;
cfg.scaleHeight_m = 60000;

% Default capture dwell for the coupled return. Other service times stay open.
cfg.captureWait_min = 5;
cfg.handoverWait_min = NaN;
cfg.veicLoading_min = NaN;
cfg.veicCharge_min = NaN;
cfg.veicHold_min = NaN;

cfg.runAnnualSK = true;
cfg.includeMomentumUnloading = true;
cfg.taxiDistance_m = 1000;
cfg.atlanteDistance_m = 250;
cfg.parkingEntryTime_s = 3600;
cfg.parkingDepartureTime_s = 3600;
% LVLH pre-burn handoff [m; m/s], with V-bar behind the reference.
cfg.taxiParkingDepartureState = [0;-100;0;0;0;0];
cfg.atlanteParkingDepartureState = [0;-100;0;0;0;0];
cfg.taxiParkingWaitBeforeDeparture_s = 0;
cfg.atlanteParkingWaitBeforeDeparture_s = 0;
defaultEllipse = struct('ellipseCenter_m',NaN, ...
    'ellipseRadialAmplitude_m',20,'ellipsePhase_rad',0);
cfg.taxiParkingEllipse = defaultEllipse;
cfg.atlanteParkingEllipse = defaultEllipse;
cfg.relativeSKDensity_kgm3 = NaN; % Evaluate the density model at parking altitude.
cfg.relativeSKInterval_hr = 6;

cfg.taxiParkingAZ_m = 300;
cfg.atlanteParkingAZ_m = NaN; % Host approach-zone radius requires a mission input.
cfg.relativeSKAZmargin_m = 0;
cfg.N_RPO_MC = 1000;
cfg.RPO_seed_client = 42;
% Finite-thrust replay of the budget representative, or all nominal pairs.
cfg.rpoFiniteValidation = "representative";
cfg.rpoClientGeometry = struct();
cfg.rpoAtlanteGeometry = struct();
cfg.rpoTerminalCorridorHalfWidth_m = Inf;
cfg.rpoMaxPropellantCycle_kg = Inf;
cfg.rpoMaxDVPerCycle_mps = Inf;
cfg.rpoMaxBurnDuration_s = Inf;
cfg.rpoMaxPeakAccel_mps2 = Inf;
end
