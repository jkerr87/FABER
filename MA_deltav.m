%% ASV MISSION ANALYSIS - DELTA V BUDGET AND TRANSFER VISUALIZATION
%
% All orbital distances: km
% Velocities and reported delta-V: m/s
% Earth gravitational parameter: km^3/s^2
%
% Transfer visualizations are ideal two-body impulsive trajectories.

clear; clc; close all;

%% ================================================================
%  1. CONSTANTS AND MISSION INPUTS
% ================================================================
% Constants
addpath('ORBMECH');
addpath('MonteCarlo');
addpath('config');
addpath('ProximityOps');

mu = astroConstants(13);
R_E = astroConstants(23);

% Target orbit --- CHANGE PER RUN ========================================
h_target = 457.25;           % [km]
%h_target = 680.12;
i_target = 51.02;            % [deg]
%i_target = 98.12;
RAAN_target = 0;             % [deg]
h_disposal = 250;

% Choose one: 'modular_station' 'large_antenna' 'solar_array' 'orbital_computing' 'truss_platform'
structureType = 'solar_array';

r_target = R_E + h_target;

% Target ranges -- CHANGE PER RUN
h_SD = 67.47;
%h_SD = 72.84;
i_SD = 2.81;
%i_SD = 0.31;
RAAN_SD = 0;

% Mission configuration
mission_years = 7;
modules_per_year = 5;

% Planning assumption only
N_modules_planned = mission_years*modules_per_year;

N_clients = 3;

%% LOAD SELECTED STRUCTURE CONFIGURATION

configFile = which([structureType '.m']);

assert(~isempty(configFile), ...
    'Structure configuration file was not found.');

run(configFile);

N_deliveries = N_containers;
N_modules = N_deliveries;


% Launcher injection errors -- FALCON 9
alt_error_launch = 15;    % [km]
inc_error_launch = 0.1;   % [deg]
RAAN_error_launch = 0.1;  % [deg]

% Ranges of target orbit for clients + launcher error
alt_err = alt_error_launch + h_SD;
inc_err = inc_error_launch + i_SD;
RAAN_err = RAAN_error_launch + RAAN_SD;

% Values for the far/close field rendezvous
N_far_legs = 1;
abort_fraction = 0.5;
CAM_per_year = 4;
N_proximity_legs = 1;
DV_abort = 1; % LITERATURE VAL = 1m/s
DV_CAM = 1; % 1 m/s LIT

% MARGINS ===============================================================
margin_fraction = 0.05;
margin_minimum = 10;
abort_margin_factor = margin_fraction;
CAM_margin_factor = margin_fraction;
momentum_margin_factor = margin_fraction;

% PROXIMITY OPERATIONS ==================================================
% VEICOLETTO PROXIMITY OPS
envelope = [40; 40; 20];
clearance = 2;
asvRadius = 3;
timesMin = [5 10 15 20];

hostPoint = [0;-50;0];       % Host location in LVLH [m] --> veicoletto starts at ASV/ container

% Values for Station Keeping
Cd_Taxi = 2.2;                     % Aerodynamic drag coefficient
mass_Taxi_dry = 8256.93;           % [kg] Estimated dry mass of the empty Taxi
mass_Veicoletto = 875.17;          % [kg] Estimated Veicoletto mass
mass_idle = mass_Taxi_dry + mass_Veicoletto;% Mass in Idle mode, with all vehicles stowed inside Taxi
density_sat = 100;                 % estimated total S/C density [kg/m^3]
vol_sat = mass_idle / density_sat;
length_sat = vol_sat^(1/3);
area_drag = length_sat^2;
rho_orb = 9.158e-12;
%rho_orb = 1.454e-13;
H_orb = 53.298;
%H_orb = 71.835;
ref_altitude_density = 9.158e-12;
%ref_altitude_density = 1.187e-14;
h_ref = 350;
%h_ref = 600;
rho = ref_altitude_density * exp(-((h_target - alt_err) - h_ref) / H_orb);
Isp_thruster_sk = 297;

% Values for momentum unloading
I_approx = (1/6) * mass_idle * length_sat^2;             % Cube inertia
delta_I = 0.1 * I_approx;                               % 10% of asymmetry as worst case
offset_cg_cp = 0.1 * length_sat;
arm_thrusters = length_sat / 2;
Cr = 1.5;                  % reflection coefficient


%% ================================================================
%  1. LAUNCH INJECTION CORRECTIONS
% ================================================================
% Workbook assumption:
% A launcher inserts the spacecraft at target altitude + error.
%
% Correction consists of:
% 1. Hohmann altitude correction
% 2. Inclination correction
% 3. RAAN correction
%
% Corrections are added separately, conservatively.
% This is not an optimized combined plane change.

r_injection = R_E + h_target + alt_error_launch;

% Altitude correction
[dva,dvb,TOF_injection] = hohmannTransfer(mu,r_injection,r_target);

DV_alt = abs(dva)+abs(dvb);

% Inclination correction
v_injection = sqrt(mu/r_injection)*1000;

DV_inc = 2*v_injection*sind(abs(inc_error_launch/2));

% RAAN correction at nominal inclination
v_target = sqrt(mu/r_target)*1000;

plane_angle = 2*asind(abs(sind(i_target)*sind(RAAN_error_launch/2)));

DV_RAAN = 2*v_target*sind(plane_angle/2);

% Separate corrections added linearly
DV_injection = DV_alt+DV_inc+DV_RAAN;

fprintf('\n====================================================\n');
fprintf('LAUNCHER INJECTION CORRECTION\n');
fprintf('====================================================\n');
fprintf('Delta V Altitude: %.2f\n',DV_alt);
fprintf('Delta V Inclination: %.2f\n',DV_inc);
fprintf('Delta V RAAN: %.2f\n',DV_RAAN);
fprintf('Delta V TOT Injection: %.2f\n',DV_injection);


%% ================================================================
%  3. DEORBIT TRANSFER
% ================================================================
% Single retrograde burn from the circular target orbit
% to an elliptical orbit with a lower perigee.
%
% No circularization at the disposal perigee.

r_initial = R_E+h_target;
r_perigee = R_E+h_disposal;

[dv_deorbit,~,tof_deorbit] = ...
    hohmannTransfer(mu,r_initial,r_perigee);

DV_deorbit = abs(dv_deorbit);

fprintf('\n====================================================\n');
fprintf('DEORBIT CALCULATION\n');
fprintf('====================================================\n');
fprintf('Initial altitude:      %.1f km\n',h_target);
fprintf('Disposal perigee:      %.1f km\n',h_disposal);
fprintf('Retrograde delta-V:   %.4f m/s\n',DV_deorbit);
fprintf('Time to perigee:      %.2f minutes\n',tof_deorbit/60);

% %% ================================================================
% %  4. VISUALIZE DEORBIT TRANSFER
% % ================================================================
% theta = linspace(0,2*pi,1000);
% 
% figure('Color','w','Name','Deorbit Transfer');
% 
% hold on; grid on; axis equal;
% 
% % Earth
% fill(R_E*cos(theta),R_E*sin(theta), ...
%      [0.25 0.50 0.80], ...
%      'EdgeColor','none','DisplayName','Earth');
% 
% % Initial orbit
% plot(r_initial*cos(theta),r_initial*sin(theta), ...
%     'b--','LineWidth',1.8, ...
%     'DisplayName','Initial circular orbit');
% 
% % Deorbit ellipse
% a_deorbit = (r_initial+r_perigee)/2;
% e_deorbit = (r_initial-r_perigee)/(r_initial+r_perigee);
% p_deorbit = a_deorbit*(1-e_deorbit^2);
% 
% nu_deorbit = linspace(pi,2*pi,500);
% 
% r_d = p_deorbit./(1+e_deorbit*cos(nu_deorbit));
% 
% x_d = r_d.*cos(nu_deorbit);
% y_d = r_d.*sin(nu_deorbit);
% 
% plot(x_d,y_d,'r','LineWidth',2.5, ...
%     'DisplayName','Deorbit transfer');
% 
% % Final perigee radius for reference
% plot(r_perigee*cos(theta),r_perigee*sin(theta), ...
%     'k:','LineWidth',1.5, ...
%     'DisplayName','Disposal perigee radius');
% 
% plot(x_d(1),y_d(1),'ko','MarkerFaceColor','r', ...
%     'MarkerSize',8,'DisplayName','Deorbit burn');
% 
% xlabel('X [km]');
% ylabel('Y [km]');
% title('ASV End-of-Life Deorbit Transfer');
% 
% legend('Location','bestoutside');

%% ================================================================
%  CLIENT ORBIT TRANSFER - ONCE PER CLIENT
% ================================================================

% Extreme client orbit altitudes
% h_SD is treated as the +/- variation around h_target

h_client_high = h_target + h_SD;
h_client_low  = h_target - h_SD;

r_client_high = R_E + h_client_high;
r_client_low  = R_E + h_client_low;

% Hohmann correction between the extreme circular orbits
[dv1_client,dv2_client,TOF_client] = ...
    hohmannTransfer(mu,r_client_high,r_client_low);

DV_alt_client = abs(dv1_client) + abs(dv2_client);

% Inclination difference between the extremes
i_client_high = i_target + i_SD + inc_error_launch;
i_client_low  = i_target - i_SD - inc_error_launch;

delta_i_client = abs(i_client_high-i_client_low);

v_client = sqrt(mu/r_client_low)*1000;

DV_inc_client = 2*v_client*sind(delta_i_client/2);

% RAAN difference between extreme client orbits
delta_RAAN_client = 2*RAAN_SD;

cos_gamma = cosd(i_client_high)*cosd(i_client_low) + ...
    sind(i_client_high)*sind(i_client_low)* ...
    cosd(delta_RAAN_client);

cos_gamma = max(-1,min(1,cos_gamma));

gamma_client = acosd(cos_gamma);

% Total plane change between orbital planes
DV_plane_client = 2*v_client*sind(gamma_client/2);

% Avoid double-counting the inclination contribution:
% DV_plane_client ALREADY includes inclination + RAAN differences.
DV_client_transfer = DV_alt_client + DV_plane_client;

fprintf('\n====================================================\n');
fprintf('WORST CASE DELTA V - NEW CLIENT\n');
fprintf('====================================================\n');
fprintf('Hohmann DV:         %.3f m/s\n',DV_alt_client);
fprintf('Plane change DV:    %.3f m/s\n',DV_plane_client);
fprintf('Total client DV:    %.3f m/s\n',DV_client_transfer);
fprintf('Hohmann TOF:        %.2f min\n',TOF_client/60);

%% ================================================================
%  FAR-FIELD RENDEZVOUS - TO CLIENT
% ================================================================
v_circular = sqrt(mu/(R_E+h_target))*1000;
RV_RAAN_err = RAAN_error_launch;

DV_far_RAAN = 2*v_circular*abs(sind(i_target))*sind(abs(RV_RAAN_err)/2);

DV_far_one_leg = DV_alt+DV_inc+DV_far_RAAN;

%DV_far_total = DV_far_one_leg*N_far_legs*2;

TOF_far = TOF_injection;  % [s]

% MATTEO'S CODE
% True Anomaly Phasing Manoeuvres

% Parameters for the corners logic
h_base = h_target;                     % Host Structure altitude
h_launcher_nominal = h_target;         % Nominal module release altitude (modify if different)
dh_err = alt_error_launch + h_SD;      % Total altitude uncertainty
phase_deg = 180;                       % Maximum true anomaly phase angle to recover 
N_rev = 100;                           % Number of orbits allowed for the maneuver

% Phasing manoeuvre function call
PH = phasing_manoeuvre(h_base, h_launcher_nominal, dh_err, phase_deg, N_rev, N_modules);

% Total Delta_V
DV_phasing = PH.DV_total_lifetime/N_modules;
TOF_phasing = PH.time_cycle_min;


%% ================================================================
% INDIVIDUAL OUTBOUND FAR-FIELD RESULTS
% ================================================================

% Preserve the existing outbound methodology
DV_alt_client_RV  = DV_alt;
DV_inc_client_RV  = DV_inc;
DV_RAAN_client_RV = DV_far_RAAN;

% Phasing to client/module
DV_phasing_to_client = PH.DV_phase1_module;

% Transfer time from current far-field estimate [min]
% Assigned to altitude-transfer row to avoid counting it
% again in the inclination and RAAN rows.
TOF_alt_client_RV = TOF_far/60;

% Impulsive plane corrections: ideal burn duration not modelled
TOF_inc_client_RV = 0;
TOF_RAAN_client_RV = 0;

TOF_phasing_to_client = PH.time_phase1_min;

fprintf('\n========================================\n');
fprintf('PHASING TO CLIENT / MODULE\n');
fprintf('========================================\n');

fprintf('Delta V: %.5f m/s\n',DV_phasing_to_client);
fprintf('TOF: %.2f min\n',TOF_phasing_to_client);
fprintf('TOF: %.2f days\n',TOF_phasing_to_client/1440);

% PRINT OUTBOUND RESULTS

fprintf('\n========================================\n');
fprintf('FAR FIELD - TO CLIENT\n');
fprintf('========================================\n');

fprintf('Altitude DV:     %.5f m/s\n',DV_alt_client_RV);
fprintf('Inclination DV:  %.5f m/s\n',DV_inc_client_RV);
fprintf('RAAN DV:         %.5f m/s\n',DV_RAAN_client_RV);
fprintf('Phasing DV:      %.5f m/s\n',DV_phasing_to_client);

fprintf('Far-field TOF estimate: %.2f min\n', ...
    TOF_alt_client_RV);

%% ================================================================
%  CLOSE-FIELD RENDEZVOUS -- TO CLIENT
% ================================================================
% TAXI PROXIMITY OPS --> Capture

% Initial state from far-field rendezvous
% [x;y;z;vx;vy;vz] in [m;m/s]
X0_RPO = [50;-2000;20;0.01;0.04;-0.01];

rpo = struct();
rpo.N = 600;

% HOMING
% End of far-field rendezvous
X0_RPO = [50;-2000;20;0.01;0.04;-0.01];

% One maneuver to stationary holding point at -1000 m
rpo.homing.targets = [0;-1000;0;0;0;0];

% Duration of the single homing maneuver
rpo.homing.times_s = 3600;  % [s]

% % Each COLUMN is the desired terminal state of one maneuver
% rpo.homing.targets = [ ...
%     0       0;
%     -1500   -1000;
%     0       0;
%     0       0;
%     0       0;
%     0       0];
% 
% % Duration of each homing maneuver [s]
% rpo.homing.times_s = [1800 2400];

% HOLDING POINT 1

rpo.holding1.mode = 'point';
rpo.holding1.param.yHold = -1000;
rpo.holding1.duration_s = 30*60;

% CLOSING
rpo.closing.targets = [0;-100;0;0;0;0];

rpo.closing.methods = {'Vbar'};

rpo.closing.times_s = 1800;

% HOLDING POINT 2

rpo.holding2.mode = 'point';
rpo.holding2.param.yHold = -100;
rpo.holding2.duration_s = 10*60;

% FINAL APPROACH

rpo.final.target = [0;-5;0;0;0;0];
rpo.final.method = 'Vbar';
rpo.final.duration_s = 600;

% RUN RPO ANALYSIS

result_rendezvous = proximity_operations_v2( ...
    mu,R_E,h_target,X0_RPO,rpo);



%% VALIDATION - SINGLE-TRANSFER HCW

cfgValidation = ASV_RPO_MC_config();

cfgValidation.mu_km3s2 = mu;
cfgValidation.Re_km = R_E;
cfgValidation.h_km = h_target;

cfgValidation.x0_client = X0_RPO;

V = FABER_RPO_validate_baseline( ...
    cfgValidation,result_rendezvous);

%% INDIVIDUAL HOMING MANEUVERS

for i = 1:result_rendezvous.N_homing

    H = result_rendezvous.homing(i);

    fprintf('Homing %d: DV = %.5f m/s, TOF = %.2f min\n', ...
        i,H.DV_mps,H.TOF_min);

end

%% INDIVIDUAL CLOSING MANEUVERS

for i = 1:result_rendezvous.N_closing

    C = result_rendezvous.closing(i);

    fprintf('Closing %d (%s): DV = %.5f m/s, TOF = %.2f min\n', ...
        i,C.method,C.DV_mps,C.TOF_min);

end

% TOTALS

DV_homing   = result_rendezvous.DV_homing;
DV_holding1 = result_rendezvous.DV_holding1;
DV_closing  = result_rendezvous.DV_closing;
DV_holding2 = result_rendezvous.DV_holding2;
DV_final    = result_rendezvous.DV_final;

TOF_homing   = result_rendezvous.TOF_homing;
TOF_holding1 = result_rendezvous.TOF_holding1;
TOF_closing  = result_rendezvous.TOF_closing;
TOF_holding2 = result_rendezvous.TOF_holding2;
TOF_final    = result_rendezvous.TOF_final;



%% RPO PHASE RESULTS

RPO_Phase = [
    "Homing"
    "Holding point 1"
    "Closing"
    "Holding point 2"
    "Final approach"
    ];

RPO_DV = [
    DV_homing
    DV_holding1
    DV_closing
    DV_holding2
    DV_final
    ];

RPO_TOF = [
    TOF_homing
    TOF_holding1
    TOF_closing
    TOF_holding2
    TOF_final
    ];

RPO_Table = table(RPO_Phase,RPO_DV,RPO_TOF, ...
    'VariableNames', ...
    {'Phase','DV_mps','TOF_min'});

fprintf('\n========================================\n');
fprintf('TAXI CLOSE-FIELD RENDEZVOUS\n');
fprintf('========================================\n');

disp(RPO_Table);

DV_RPO_total = sum(RPO_DV);
TOF_RPO_total = sum(RPO_TOF);

fprintf('Total RPO DV: %.5f m/s\n',DV_RPO_total);
fprintf('Total RPO TOF: %.2f min\n',TOF_RPO_total);





%% ================================================================
% FAR-FIELD RENDEZVOUS - RETURN TO HOST
% ================================================================
% Taxi carrying captured module -> Host structure
%
% 1. Estimate J2 differential RAAN precession
% 2. Calculate altitude correction
% 3. Calculate inclination correction
% 4. Calculate RAAN correction using propagated RAAN
%
% Assumptions:
% - Circular reference orbits
% - Constant orbital elements except J2 RAAN precession
% - Both initial RAAN values refer to the same epoch
% - Impulsive, separate plane corrections

% 1. DEFINE ORBITS
J2 = 1.08262668e-3;

% HOST ORBIT
h_host = h_target;
i_host = i_target;
RAAN_host_0 = RAAN_target;

% CAPTURED MODULE ORBIT
% Worst-case altitude selected by Matteo's model
h_module = PH.worst_h_launcher;

% Example orbital-plane offsets for captured module.
% Replace with actual module orbital elements.
i_module = i_target + inc_error_launch;
RAAN_module_0 = RAAN_target + RAAN_error_launch;

% Circular-orbit assumption
e_host = 0;
e_module = 0;

a_host = R_E + h_host;       % [km]
a_module = R_E + h_module;   % [km]

% 2. ESTIMATE DIFFERENTIAL J2 PRECESSION
n_host = sqrt(mu/a_host^3);
n_module = sqrt(mu/a_module^3);

% RAAN precession rates [rad/s]
OmegaDot_host = ...
    -1.5*J2*n_host*(R_E/a_host)^2 * ...
    cosd(i_host)/(1-e_host^2)^2;

OmegaDot_module = ...
    -1.5*J2*n_module*(R_E/a_module)^2 * ...
    cosd(i_module)/(1-e_module^2)^2;

% Convert to degrees per second
OmegaDot_host_deg = rad2deg(OmegaDot_host);
OmegaDot_module_deg = rad2deg(OmegaDot_module);

% ESTIMATED ELAPSED TIME TO RETURN DEPARTURE

T_outbound = ...
    TOF_alt_client_RV*60 + ...
    TOF_phasing_to_client*60 + ...
    result_rendezvous.TOF_total*60;

% Capture / handling time is not modelled yet
T_capture = 0;

T_returnDeparture = T_outbound + T_capture;

% Preliminary J2 sensitivity assumption (time taken to capture module)
% Both orbits are treated as fixed for this interval
T_separation = T_returnDeparture; % [s]

fprintf('Elapsed time to return departure: %.2f days\n', ...
    T_returnDeparture/86400);

% Propagate RAAN to return departure epoch
RAAN_host_return = mod( ...
    RAAN_host_0 + OmegaDot_host_deg*T_separation,360);

RAAN_module_return = mod( ...
    RAAN_module_0 + OmegaDot_module_deg*T_separation,360);

% Shortest signed RAAN difference [-180,180]
delta_RAAN_return = mod( ...
    RAAN_host_return - RAAN_module_return + 180,360) - 180;

% 3. ALTITUDE CORRECTION
r_module = R_E + h_module;
r_host = R_E + h_host;

[dv1_host,dv2_host,TOF_alt_host_s] = ...
    hohmannTransfer(mu,r_module,r_host);

DV_alt_host = abs(dv1_host) + abs(dv2_host);

TOF_alt_host = TOF_alt_host_s/60; % [min]

% 4. INCLINATION CORRECTION
delta_i_return = i_host - i_module;

% Approximate sequential plane correction at host radius
v_host = sqrt(mu/r_host)*1000; % [m/s]

DV_inc_host = ...
    2*v_host*sind(abs(delta_i_return)/2);

TOF_inc_host = 0; % Ideal impulsive burn duration

% 5. RAAN CORRECTION USING J2 DRIFT
% RAAN-only plane rotation after inclination correction
% to the host inclination.

gamma_RAAN = 2*asind( ...
    min(1,abs(sind(i_host)* ...
    sind(delta_RAAN_return/2))));

DV_RAAN_host = ...
    2*v_host*sind(gamma_RAAN/2);

TOF_RAAN_host = 0; % Ideal impulsive burn duration

% 6. PRINT INDIVIDUAL RESULTS
fprintf('\n========================================\n');
fprintf('RETURN TO HOST - ORBITAL DIFFERENCES\n');
fprintf('========================================\n');

fprintf('Module altitude: %.3f km\n',h_module);
fprintf('Host altitude:   %.3f km\n',h_host);
fprintf('Altitude difference: %.3f km\n',h_host-h_module);

fprintf('\nModule inclination: %.4f deg\n',i_module);
fprintf('Host inclination:   %.4f deg\n',i_host);
fprintf('Inclination difference: %.4f deg\n',delta_i_return);

fprintf('\nHost RAAN rate:   %.5f deg/day\n', ...
    OmegaDot_host_deg*86400);

fprintf('Module RAAN rate: %.5f deg/day\n', ...
    OmegaDot_module_deg*86400);

fprintf('Return epoch: %.3f days after reference\n', ...
    T_separation/86400);

fprintf('Host RAAN at return:   %.5f deg\n', ...
    RAAN_host_return);

fprintf('Module RAAN at return: %.5f deg\n', ...
    RAAN_module_return);

fprintf('Relative RAAN at return: %.5f deg\n', ...
    delta_RAAN_return);

fprintf('\n========================================\n');
fprintf('RETURN TO HOST - DELTA V\n');
fprintf('========================================\n');

fprintf('Altitude correction:    %.5f m/s\n',DV_alt_host);
fprintf('Inclination correction: %.5f m/s\n',DV_inc_host);
fprintf('RAAN correction:        %.5f m/s\n',DV_RAAN_host);



% RETURN PHASING - TO HOST STRUCTURE
% Matteo's phasing function
DV_phasing_to_host = PH.DV_phase2_base;
TOF_phasing_to_host = PH.time_phase2_min;

%% ================================================================
% FAR-FIELD RETURN TOTALS
% ================================================================

% Altitude, inclination and RAAN corrections
DV_far_to_host = ...
    DV_alt_host + DV_inc_host + DV_RAAN_host;

% Total return including phasing
DV_return_to_host = ...
    DV_far_to_host + DV_phasing_to_host;

% Altitude transfer TOF was already computed
% using hohmannTransfer:
% TOF_alt_host = TOF_alt_host_s/60;

% Impulsive plane corrections
TOF_inc_host = 0;
TOF_RAAN_host = 0;

% ================================================================
% OUTBOUND VALUES - FOR COMPARISON
% ================================================================

DV_far_to_client = DV_far_one_leg;
TOF_far_to_client = TOF_alt_client_RV;

DV_phasing_to_client = PH.DV_phase1_module;
TOF_phasing_to_client = PH.time_phase1_min;

DV_outbound_to_client = ...
    DV_far_to_client + DV_phasing_to_client;

% ================================================================
% COMPLETE RETRIEVAL CYCLE - SUMMARY ONLY
% ================================================================

DV_far_cycle = ...
    DV_far_to_client + DV_far_to_host;

DV_phasing_cycle = ...
    DV_phasing_to_client + DV_phasing_to_host;

DV_transfer_cycle = ...
    DV_far_cycle + DV_phasing_cycle;

TOF_far_cycle = ...
    TOF_alt_client_RV + TOF_alt_host;

TOF_phasing = PH.time_cycle_min;

% ================================================================
% PRINT RETURN RESULTS
% ================================================================

fprintf('\n========================================\n');
fprintf('FAR FIELD - RETURN TO HOST\n');
fprintf('========================================\n');

fprintf('Altitude DV:      %.5f m/s\n',DV_alt_host);
fprintf('Inclination DV:   %.5f m/s\n',DV_inc_host);
fprintf('RAAN DV:          %.5f m/s\n',DV_RAAN_host);
fprintf('Phasing DV:       %.5f m/s\n',DV_phasing_to_host);

fprintf('\nAltitude TOF:     %.2f min\n',TOF_alt_host);
fprintf('Phasing TOF:      %.2f min\n',TOF_phasing_to_host);
fprintf('Phasing TOF:      %.2f days\n', ...
    TOF_phasing_to_host/1440);

fprintf('\nTotal return DV:  %.5f m/s\n',DV_return_to_host);

% ================================================================
% PRINT PHASING CYCLE
% ================================================================

fprintf('\n========================================\n');
fprintf('PHASING - COMPLETE RETRIEVAL CYCLE\n');
fprintf('========================================\n');

fprintf('To client DV: %.5f m/s\n',DV_phasing_to_client);
fprintf('To client TOF: %.2f days\n', ...
    TOF_phasing_to_client/1440);

fprintf('To host DV: %.5f m/s\n',DV_phasing_to_host);
fprintf('To host TOF: %.2f days\n', ...
    TOF_phasing_to_host/1440);

fprintf('Complete phasing TOF: %.2f days\n', ...
    TOF_phasing/1440);


%% ================================================================
%  CLOSE FIELD RENDEZVOUS - TO HOST STRUCTURE
% ================================================================
% TAXI PROXIMITY OPS --> Back to host structure + holding pt/orbit

% TAXI PROXIMITY OPERATIONS --> RETURN TO HOST
%
% Preliminary baseline:
% - Same homing/closing profile as container approach
% - All coordinates relative to HOST
% - Different initial state can be specified
% - Host final approach terminates at a stationary hold

%% HOST-RELATIVE INITIAL STATE

% Placeholder arrival state after far-field return.
% Replace with actual far-field terminal state once available.
X0_host = [50;-2000;20;0.01;0.04;-0.01];

%% HOST APPROACH CONFIGURATION

rpo_host = rpo;

% Homing to first holding point
rpo_host.homing.targets = [0;-1000;0;0;0;0];
rpo_host.homing.times_s = 3600;

% Holding point 1
rpo_host.holding1.mode = 'point';
rpo_host.holding1.param.yHold = -1000;
rpo_host.holding1.duration_s = 30*60;

% Closing to second holding point
rpo_host.closing.targets = [0;-100;0;0;0;0];
rpo_host.closing.methods = {'Vbar'};
rpo_host.closing.times_s = 1800; % Ignored for Vbar

% Holding point 2
rpo_host.holding2.mode = 'point';
rpo_host.holding2.param.yHold = -100;
rpo_host.holding2.duration_s = 10*60;

% Final approach to host capture/docking-ready point
rpo_host.final.target = [0;-5;0;0;0;0];
rpo_host.final.method = 'Vbar';
rpo_host.final.duration_s = 600;

%% RUN HOST RENDEZVOUS

result_host = proximity_operations_v2( ...
    mu,R_E,h_target,X0_host,rpo_host);

%% INDIVIDUAL HOST HOMING MANEUVERS

for i = 1:result_host.N_homing

    H = result_host.homing(i);

    fprintf('Host homing %d: DV = %.5f m/s, TOF = %.2f min\n', ...
        i,H.DV_mps,H.TOF_min);

end

%% INDIVIDUAL HOST CLOSING MANEUVERS

for i = 1:result_host.N_closing

    C = result_host.closing(i);

    fprintf('Host closing %d (%s): DV = %.5f m/s, TOF = %.2f min\n', ...
        i,C.method,C.DV_mps,C.TOF_min);

end

%% HOST PHASE TOTALS

DV_host_homing   = result_host.DV_homing;
DV_host_holding1 = result_host.DV_holding1;
DV_host_closing  = result_host.DV_closing;
DV_host_holding2 = result_host.DV_holding2;
DV_host_final    = result_host.DV_final;

TOF_host_homing   = result_host.TOF_homing;
TOF_host_holding1 = result_host.TOF_holding1;
TOF_host_closing  = result_host.TOF_closing;
TOF_host_holding2 = result_host.TOF_holding2;
TOF_host_final    = result_host.TOF_final;



%% ================================================================
% TOTAL RPO - CLIENT AND HOST
% ================================================================

%% CLIENT PHASE TOTALS

DV_client_homing   = result_rendezvous.DV_homing;
DV_client_holding1 = result_rendezvous.DV_holding1;
DV_client_closing  = result_rendezvous.DV_closing;
DV_client_holding2 = result_rendezvous.DV_holding2;
DV_client_final    = result_rendezvous.DV_final;

TOF_client_homing   = result_rendezvous.TOF_homing;
TOF_client_holding1 = result_rendezvous.TOF_holding1;
TOF_client_closing  = result_rendezvous.TOF_closing;
TOF_client_holding2 = result_rendezvous.TOF_holding2;
TOF_client_final    = result_rendezvous.TOF_final;

%% ROUND-TRIP PHASE TOTALS

DV_homing = DV_client_homing + DV_host_homing;
DV_holding1 = DV_client_holding1 + DV_host_holding1;
DV_closing = DV_client_closing + DV_host_closing;
DV_holding2 = DV_client_holding2 + DV_host_holding2;
DV_final = DV_client_final + DV_host_final;

TOF_homing = TOF_client_homing + TOF_host_homing;
TOF_holding1 = TOF_client_holding1 + TOF_host_holding1;
TOF_closing = TOF_client_closing + TOF_host_closing;
TOF_holding2 = TOF_client_holding2 + TOF_host_holding2;
TOF_final = TOF_client_final + TOF_host_final;

%% RPO TOTALS

DV_RPO_client = result_rendezvous.DV_total;
DV_RPO_host = result_host.DV_total;

TOF_RPO_client = result_rendezvous.TOF_total;
TOF_RPO_host = result_host.TOF_total;

DV_RPO_cycle = DV_RPO_client + DV_RPO_host;
TOF_RPO_cycle = TOF_RPO_client + TOF_RPO_host;

%% PRINT EACH PHASE

Phase = [
    "Homing"
    "Holding point 1"
    "Closing"
    "Holding point 2"
    "Final approach"
];

DV_Client = [
    DV_client_homing
    DV_client_holding1
    DV_client_closing
    DV_client_holding2
    DV_client_final
];

DV_Host = [
    DV_host_homing
    DV_host_holding1
    DV_host_closing
    DV_host_holding2
    DV_host_final
];

TOF_Client = [
    TOF_client_homing
    TOF_client_holding1
    TOF_client_closing
    TOF_client_holding2
    TOF_client_final
];

TOF_Host = [
    TOF_host_homing
    TOF_host_holding1
    TOF_host_closing
    TOF_host_holding2
    TOF_host_final
];

RPO_CycleTable = table( ...
    Phase,DV_Client,DV_Host,DV_Client+DV_Host, ...
    TOF_Client,TOF_Host,TOF_Client+TOF_Host, ...
    'VariableNames', { ...
    'Phase','DV_Client_mps','DV_Host_mps', ...
    'DV_Cycle_mps','TOF_Client_min','TOF_Host_min', ...
    'TOF_Cycle_min'});

fprintf('\n========================================\n');
fprintf('TAXI RPO - COMPLETE RETRIEVAL CYCLE\n');
fprintf('========================================\n');

disp(RPO_CycleTable);

fprintf('Client RPO DV: %.5f m/s\n',DV_RPO_client);
fprintf('Host RPO DV:   %.5f m/s\n',DV_RPO_host);
fprintf('Total RPO DV:  %.5f m/s\n',DV_RPO_cycle);

fprintf('Client RPO TOF: %.2f min\n',TOF_RPO_client);
fprintf('Host RPO TOF:   %.2f min\n',TOF_RPO_host);
fprintf('Total RPO TOF:  %.2f min\n',TOF_RPO_cycle);

%% ================================================================
%  13. TAXI HOLD WHILE ASSEMBLING (assumed elliptical)
% ================================================================

%% TAXI HOLDING DURING ASSEMBLY

n = sqrt(mu/(R_E+h_target)^3);

param = struct();
param.yCenter = -1000;
param.radAmp = 75;
param.phase = 0;

T_assemblyHolding = 7*24*3600; % Example: seven days

holding = rpo_holdingpt( ...
    n,T_assemblyHolding,'ellipse',param,600);

DV_Taxi_holding = holding.DV_holding_mps;
TOF_Taxi_holding = holding.TOF_holding_min;


%% ================================================================
%  STRUCTURE CONFIGURATION - CHANGE PER RUN
% ================================================================

% One installation relocation per structural element
N_moduleVisits = size(modulePoints,1);

fprintf('\n========================================\n');
fprintf('SELECTED ASSEMBLY ARCHITECTURE\n');
fprintf('========================================\n');

fprintf('Structure: %s\n',structureName);
fprintf('Envelope volume: %.1f m^3\n',structureVolume);

fprintf('Large modules: %d\n',N_largeModules);
fprintf('Small elements: %d\n',N_smallElements);

fprintf('Small elements/container: %d (capacity)\n', ...
    N_smallPerFairing);

fprintf('Launches/containers: %d\n',N_containers);
fprintf('Assembly relocation trips: %d\n',N_moduleVisits);

% Container contents

for k = 1:N_containers

    fprintf('Container %d: %d elements\n', ...
        k,elementsPerContainer(k));

end
%% ================================================================
%  ASSEMBLY RELOCATIONS
% ================================================================
% Veicoletto PROXIMITY OPS --> just consistent "final approaches"

% VEICOLETTO PROXIMITY OPS
assemblyResult = assembly_relocations_v2( ...
    mu, R_E, h_target, 101, ...
    timesMin, hostPoint, modulePoints);


%% ASSEMBLY RELOCATION RESULTS

% Aggregate delta-V for all elements
DV_assembly_out = sum(assemblyResult.DV_outbound);
DV_assembly_return = sum(assemblyResult.DV_return);

DV_assembly_relocation = assemblyResult.DV;

% Aggregate transfer times
% Both legs currently use the same duration
TOF_assembly_out = sum(assemblyResult.timePerLeg_min);
TOF_assembly_return = sum(assemblyResult.timePerLeg_min);

TOF_assembly_relocation = assemblyResult.TOF_min;

% Assembly operations
N_assemblyOperations = assemblyResult.Nmodules;

% PRINT ASSEMBLY SUMMARY

fprintf('\n========================================\n');
fprintf('VEICOLETTO ASSEMBLY RELOCATIONS\n');
fprintf('========================================\n');

fprintf('Number of elements: %d\n',N_assemblyOperations);

fprintf('\nOutbound transfers:\n');
fprintf('DV: %.5f m/s\n',DV_assembly_out);
fprintf('TOF: %.2f min\n',TOF_assembly_out);

fprintf('\nReturn transfers:\n');
fprintf('DV: %.5f m/s\n',DV_assembly_return);
fprintf('TOF: %.2f min\n',TOF_assembly_return);

fprintf('\nTOTAL ASSEMBLY TRANSFERS\n');
fprintf('DV: %.5f m/s\n',DV_assembly_relocation);
fprintf('TOF: %.2f min\n',TOF_assembly_relocation);


%% ================================================================
%  12. STATION KEEPING (Taxi in Parking Orbit - Worst Case)
% ================================================================

% Station keeping orbit altitude
SK_altitude = h_target - alt_err;

% Function call
sk_results = station_keeping(mass_idle, rho, area_drag, Cd_Taxi, SK_altitude, mission_years);

% Assignment for the final budget
DV_SK = sk_results.DV_SK;
DV_stationkeeping = DV_SK;
% Results output
fprintf('\n====================================================\n');
fprintf('STATION KEEPING CALCULATION\n');
fprintf('====================================================\n');
fprintf('Total system mass in idle mode: %.2f kg\n', mass_idle);
fprintf('Frontal area: %.4f m^2\n', area_drag);
fprintf('Air density: %.4f kg/m^3\n', rho);
fprintf('Station Keeping DeltaV (Worst Case, %d years): %.2f m/s\n', mission_years, DV_SK);


%% ================================================================
%  12. MOMENTUM UNLOADING
% ================================================================

% Function call
DV_mom_unloading = momentum_unloading(mass_idle, rho, SK_altitude, area_drag, Cd_Taxi, offset_cg_cp, Cr, delta_I, mission_years, arm_thrusters);

% Results
fprintf('Max Disturbance Torque:     %.3e Nm\n', DV_mom_unloading.T_total_disturbance);
fprintf('Delta_V (%d years):         %.3f m/s\n', mission_years, DV_mom_unloading.DV);

DV_momentum_unloading = DV_mom_unloading.DV;



%% ================================================================
% 13. SEVEN-YEAR DELTA V AND TOF BUDGET
% ================================================================
% All maneuvers are listed INDIVIDUALLY.
%
% Outbound = Taxi to client/container
% Inbound  = Taxi return to host
%
% DV [m/s]
% TOF [min]
%
% RPO entries represent one approach in the named direction.
% Assembly entries already aggregate all element transfers.

%% EVENT COUNTS

N_aborts = abort_fraction*N_deliveries;
N_CAM = CAM_per_year*mission_years;

N_assemblyCampaigns = 1;

TOF_abort = NaN;
TOF_CAM = NaN;

%% ================================================================
% MANEUVER NAMES
% ================================================================

maneuver_names = {
    'Initial ASV injection correction'    % 1
    'Client orbit transfer'              % 2

    'To client - Altitude correction'     % 3
    'To client - Inclination correction'  % 4
    'To client - RAAN correction'         % 5
    'To client - Phasing'                 % 6

    'To client - Homing'                  % 7
    'To client - Holding point 1'         % 8
    'To client - Closing'                 % 9
    'To client - Holding point 2'         % 10
    'To client - Final approach'          % 11

    'To host - Altitude correction'       % 12
    'To host - Inclination correction'    % 13
    'To host - RAAN correction'           % 14
    'To host - Phasing'                   % 15

    'To host - Homing'                    % 16
    'To host - Holding point 1'           % 17
    'To host - Closing'                   % 18
    'To host - Holding point 2'           % 19
    'To host - Final approach'            % 20

    'Aborts'                              % 21
    'Collision avoidance'                 % 22

    'Assembly outbound'                   % 23
    'Assembly return'                     % 24

    'Stationkeeping'                      % 25
    'Momentum unloading'                  % 26
    'Disposal'                            % 27
};

%% ================================================================
% RAW DELTA V [m/s]
% ================================================================

DV_raw = [
    DV_injection
    DV_client_transfer

    % TO CLIENT - FAR FIELD
    DV_alt_client_RV
    DV_inc_client_RV
    DV_RAAN_client_RV
    DV_phasing_to_client

    % TO CLIENT - CLOSE FIELD
    DV_client_homing
    DV_client_holding1
    DV_client_closing
    DV_client_holding2
    DV_client_final

    % TO HOST - FAR FIELD
    DV_alt_host
    DV_inc_host
    DV_RAAN_host
    DV_phasing_to_host

    % TO HOST - CLOSE FIELD
    DV_host_homing
    DV_host_holding1
    DV_host_closing
    DV_host_holding2
    DV_host_final

    % CONTINGENCIES
    DV_abort
    DV_CAM

    % ASSEMBLY
    DV_assembly_out
    DV_assembly_return

    % LIFETIME OPERATIONS AND DISPOSAL
    DV_stationkeeping
    DV_momentum_unloading
    DV_deorbit
];

%% ================================================================
% TIME OF FLIGHT [min]
% ================================================================

TOF_perEvent_min = [
    TOF_injection/60
    TOF_client/60

    % TO CLIENT - FAR FIELD
    TOF_alt_client_RV
    TOF_inc_client_RV
    TOF_RAAN_client_RV
    TOF_phasing_to_client

    % TO CLIENT - CLOSE FIELD
    TOF_client_homing
    TOF_client_holding1
    TOF_client_closing
    TOF_client_holding2
    TOF_client_final

    % TO HOST - FAR FIELD
    TOF_alt_host
    TOF_inc_host
    TOF_RAAN_host
    TOF_phasing_to_host

    % TO HOST - CLOSE FIELD
    TOF_host_homing
    TOF_host_holding1
    TOF_host_closing
    TOF_host_holding2
    TOF_host_final

    % CONTINGENCIES
    TOF_abort
    TOF_CAM

    % ASSEMBLY
    TOF_assembly_out
    TOF_assembly_return

    % LIFETIME OPERATIONS AND DISPOSAL
    NaN
    NaN
    tof_deorbit/60
];

%% ================================================================
% NUMBER OF EVENTS
% ================================================================

N_events = [
    1                    % Injection
    N_clients            % Client orbit transfer

    % TO CLIENT - FAR FIELD
    N_deliveries
    N_deliveries
    N_deliveries
    N_deliveries

    % TO CLIENT - CLOSE FIELD
    N_deliveries
    N_deliveries
    N_deliveries
    N_deliveries
    N_deliveries

    % TO HOST - FAR FIELD
    N_deliveries
    N_deliveries
    N_deliveries
    N_deliveries

    % TO HOST - CLOSE FIELD
    N_deliveries
    N_deliveries
    N_deliveries
    N_deliveries
    N_deliveries

    % CONTINGENCIES
    N_aborts
    N_CAM

    % ASSEMBLY (already aggregated)
    N_assemblyCampaigns
    N_assemblyCampaigns

    % LIFETIME OPERATIONS AND DISPOSAL
    1
    1
    1
];

%% ================================================================
% VALIDATION
% ================================================================

assert(numel(maneuver_names)==numel(DV_raw), ...
    'Maneuver names and Delta V lengths do not match.');

assert(numel(DV_raw)==numel(TOF_perEvent_min), ...
    'Delta V and TOF lengths do not match.');

assert(numel(DV_raw)==numel(N_events), ...
    'Delta V and event-count lengths do not match.');

assert(all(isfinite(DV_raw)), ...
    'Undefined Delta V detected.');

assert(all(DV_raw>=0), ...
    'Delta V must be non-negative.');

%% ================================================================
% DELTA V MARGINS
% ================================================================
% Default fractional margin for each maneuver.
% Minimum allocation applied only to selected major events.
%
% The chosen policy is preliminary and not Monte Carlo-derived.

DV_margin = DV_raw.*(1+margin_fraction);

majorNames = [
    "Initial ASV injection correction"
    "Client orbit transfer"
    "Disposal"
];

majorIdx = ismember(string(maneuver_names),majorNames);

DV_margin(majorIdx) = DV_raw(majorIdx) + ...
    max(margin_fraction.*DV_raw(majorIdx), ...
        margin_minimum);

% Special contingency multipliers
abortIdx = strcmp(maneuver_names,'Aborts');
momentumIdx = strcmp(maneuver_names,'Momentum unloading');

DV_margin(abortIdx) = 2*DV_raw(abortIdx);
DV_margin(momentumIdx) = 2*DV_raw(momentumIdx);

% Ideal zero-DV holds remain zero.
% Real holding control must be allocated separately.

%% ================================================================
% LIFETIME TOTALS
% ================================================================

DV_lifetime_raw = DV_raw.*N_events;

DV_lifetime_margin = DV_margin.*N_events;

TOF_lifetime_min = TOF_perEvent_min.*N_events;

%% ================================================================
% BUDGET TABLE
% ================================================================

BudgetTable = table( ...
    string(maneuver_names), ...
    N_events, ...
    DV_raw, ...
    DV_margin, ...
    DV_lifetime_raw, ...
    DV_lifetime_margin, ...
    TOF_perEvent_min, ...
    TOF_lifetime_min, ...
    'VariableNames', { ...
    'Maneuver', ...
    'Count', ...
    'RawPerEvent_mps', ...
    'MarginedPerEvent_mps', ...
    'RawLifetime_mps', ...
    'MarginedLifetime_mps', ...
    'TOF_perEvent_min', ...
    'TOF_lifetime_min'});

%% ================================================================
% PRINT EACH EVENT AND MISSION TOTAL
% ================================================================

fprintf('\n====================================================\n');
fprintf('ASV MISSION - INDIVIDUAL DELTA V BUDGET\n');
fprintf('====================================================\n');

disp(BudgetTable);

DV_total_raw = sum(DV_lifetime_raw);
DV_total_margin = sum(DV_lifetime_margin);

% TOF excludes phases with unknown duration
knownTOF = TOF_lifetime_min(isfinite(TOF_lifetime_min));

TOF_modelled_min = sum(knownTOF);
TOF_modelled_hr = TOF_modelled_min/60;
TOF_modelled_days = TOF_modelled_min/(60*24);

fprintf('\n========================================\n');
fprintf('MISSION TOTALS\n');
fprintf('========================================\n');

fprintf('Raw lifetime DV:      %.5f m/s\n',DV_total_raw);
fprintf('Margined lifetime DV: %.5f m/s\n',DV_total_margin);

fprintf('Modelled TOF: %.2f hours\n',TOF_modelled_hr);
fprintf('Modelled TOF: %.2f days\n',TOF_modelled_days);

if any(isnan(TOF_lifetime_min))
    fprintf(['WARNING: Some maneuver durations are ', ...
        'not included in the modelled TOF.\n']);
end

%% ================================================================
% OPTIONAL: OUTBOUND VS INBOUND SUMMARY
% ================================================================
% This is for comparison only.
% These values are ALREADY included in the main BudgetTable.

outboundIdx = 3:11;
inboundIdx = 12:20;

DV_outbound_cycle = sum(DV_raw(outboundIdx));
DV_inbound_cycle = sum(DV_raw(inboundIdx));

fprintf('\n========================================\n');
fprintf('TAXI TRANSFER COMPARISON - PER RETRIEVAL\n');
fprintf('========================================\n');

fprintf('Outbound to client: %.5f m/s\n', ...
    DV_outbound_cycle);

fprintf('Inbound to host:    %.5f m/s\n', ...
    DV_inbound_cycle);

%% ================================================================
%  14. PLOT SEVEN-YEAR DELTA V BUDGET
% ================================================================

figure('Color','w','Name','Mission Delta-V Budget');

barh(DV_lifetime_margin);

set(gca,'YTick',1:numel(maneuver_names), ...
    'YTickLabel',maneuver_names);

set(gca,'YDir','reverse');

xlabel('Seven-year allocated delta-V [m/s]');

title('ASV Mission Delta-V Budget - Workbook Assumptions');

grid on;

%% ================================================================
%  15. COMPARISON: RAW VS MARGINED BUDGET
% ================================================================

figure('Color','w','Name','Raw vs Margined Budget');

barh([DV_lifetime_raw DV_lifetime_margin]);

set(gca,'YTick',1:numel(maneuver_names), ...
    'YTickLabel',maneuver_names, ...
    'YDir','reverse');

xlabel('Lifetime allocation [m/s]');

title('ASV Delta-V Budget: Raw vs Margined');

legend('Raw','Margined','Location','best');

grid on;

%% ================================================================
%  16. SAVE SUMMARY TABLES
% ================================================================

writetable(BudgetTable, ...
    'ASV_SevenYear_DeltaV_Budget.csv');




%% ================================================================
% RPO MONTE CARLO - OPTIMIZATION STUDY
% ================================================================

cfg = ASV_RPO_MC_config();

% Same orbital parameters as mission analysis
cfg.mu_km3s2 = mu;
cfg.Re_km = R_E;
cfg.h_km = h_target;

% Same initial conditions
cfg.x0_client = X0_RPO;
cfg.x0_host = X0_host;

% Run
MC_RPO = ASV_RPO_MonteCarlo(cfg);

% Analyze selected transfer counts and Pareto frontier
A = FABER_RPO_analyze_results(MC_RPO);

% P95 results
DV_P95 = MC_RPO.P95_DV_cycle_mps;
TOF_P95 = MC_RPO.P95_TOF_cycle_min;


%% ANALYZE TRANSFER COUNT VS DELTA V

T = MC_RPO.trialTable;

fprintf('\nHOMING TRANSFER COUNT ANALYSIS\n');

for Nh = 1:4

    idx = T.ClientN_Homing == Nh;

    if any(idx)

        fprintf(['Nh = %d | Samples = %d | ', ...
            'Mean Cycle DV = %.3f m/s | ', ...
            'Mean Cycle TOF = %.1f min\n'], ...
            Nh,nnz(idx), ...
            mean(T.CycleDV_mps(idx)), ...
            mean(T.CycleTOF_min(idx)));

    end


end

% Plot transfer count vs complete-cycle Delta V
figure('Color','w');

scatter(T.ClientN_Homing,T.CycleDV_mps,35,'filled');

xlabel('Number of client homing transfers');
ylabel('Complete retrieval cycle Delta V [m/s]');
title('Influence of Homing Transfer Count on Delta V');
grid on;

T = MC_RPO.trialTable;
DV = MC_RPO.DV_phase_mps;

fprintf('\nCLIENT HOMING COUNT VS HOMING DV\n');

for Nh = 1:4
    idx = T.ClientN_Homing == Nh;

    if any(idx)
        fprintf('Nh=%d | N=%d | Mean DV=%.4f m/s\n', ...
            Nh, nnz(idx), mean(DV(idx,1)));
    end
end


%% ================================================================
%  LOCAL FUNCTIONS
% ================================================================

function [dv1,dv2,tof] = hohmannTransfer(mu,r1,r2)
% HOHMANNTRANSFER
% Inputs:
%   mu       gravitational parameter [km^3/s^2]
%   r1, r2   initial and final circular-orbit radii [km]
%
% Outputs:
%   dv1      signed first burn [m/s]
%   dv2      signed second burn [m/s]
%   tof      transfer time [s]

    if r1 <= 0 || r2 <= 0
        error('Orbit radii must be positive.');
    end

    a = (r1+r2)/2;

    v_c1 = sqrt(mu/r1);
    v_c2 = sqrt(mu/r2);

    v_t1 = sqrt(mu*(2/r1-1/a));
    v_t2 = sqrt(mu*(2/r2-1/a));

    dv1 = (v_t1-v_c1)*1000;
    dv2 = (v_c2-v_t2)*1000;

    tof = pi*sqrt(a^3/mu);

end
