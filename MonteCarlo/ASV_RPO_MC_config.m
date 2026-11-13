
function cfg = ASV_RPO_MC_config()
% ASV_RPO_MC_CONFIG
% FABER Taxi - Phase A RPO Monte Carlo
%
% Far-field rendezvous ends at approximately 2 km.
% The following seven phases are modelled:
%
% 1. Homing 1:      ~2 km -> HP1
% 2. Holding 1:     HP1
% 3. Homing 2:      HP1 -> HP2
% 4. Holding 2:     HP2
% 5. Closing:       HP2 -> HP3
% 6. Holding 3:     HP3
% 7. Final approach: HP3 -> capture/docking
%
% All distance and speed bounds are preliminary.
% Units: m, m/s, min unless otherwise indicated.

%% 1. MONTE CARLO

cfg.N = 100;
cfg.seed = 42;
cfg.maxAttemptsFactor = 30;
cfg.nTrajectorySamples = 101;

cfg.outputFolder = 'ASV_RPO_MC_results';
cfg.makePlots = false;
cfg.saveCSV = false;

%% 2. OPTIMIZATION SETTINGS (FUTURE USE)

% Current implementation samples feasible trajectories.
% It does not perform objective-based optimization.

cfg.objective = 'dv';
cfg.methods = {'Vbar'};
cfg.candidatesPerLeg = 3;

cfg.weightDV = 0.5;
cfg.weightTOF = 0.5;
cfg.refDV_mps = 2.0;
cfg.refTOF_min = 400;

%% 3. ORBIT

cfg.mu_km3s2 = 398600.4418;
cfg.Re_km = 6378.1363;
cfg.h_km = 457.25;

%% 4. INITIAL RELATIVE STATES

% Far-field to RPO handover at approximately 2 km.
% LVLH: [x;y;z;vx;vy;vz]
% Position in m, velocity in m/s.

cfg.x0_client = [50;-2000;20;0.01;0.04;-0.01];
cfg.x0_host   = [50;-2000;20;0.01;0.04;-0.01];

% Initial state dispersions
cfg.dx0_client = [100;300;100;0.05;0.05;0.05];
cfg.dx0_host   = [100;300;100;0.05;0.05;0.05];

%% 5. HOLDING POINT GEOMETRY

% Keep your original waypoint distances.

cfg.hp1_m = [800 1200];  % HP1, approximately 1 km
cfg.hp2_m = [350 450];
cfg.hp3_m = [20 20];     % HP3, terminal gate
cfg.final_m = [5 5];    % Capture/docking standoff

% Safety geometry
cfg.AZ_m = 300;
cfg.KOZ_m = 30;
cfg.corridorHalfWidth_m = 5;
cfg.maxClosingExcursionFactor = 1.50;

%% 6. MANEUVER TIME BOUNDS [min]

% HOMING 1: approximately 2 km -> HP1
cfg.homing_min = [30 100];

% HOLDING POINT 1
cfg.hold1_min = [10 60];

% HOMING 2: HP1 -> HP2
%
% Existing Monte Carlo internal name: closing_min
cfg.closing_min = [20 130];

% HOLDING POINT 2
cfg.hold2_min = [5 45];

% CLOSING: HP2 -> HP3
%
% Existing Monte Carlo internal name: controlled_min
cfg.controlled_min = [10 70];

% HOLDING POINT 3
cfg.hold3_min = [5 15];

% FINAL APPROACH: HP3 -> 5 m
cfg.final_min = [5 30];

% Minimum transfer duration per HCW arc [s]
cfg.minTransferTime_s = 120;

%% 7. NUMBER OF HCW ARCS

% Exactly TWO homing arcs in total:
%
% Arc 1: initial state -> HP1
% Arc 2: HP1 -> HP2
%
% One HCW transfer in each phase.
% No random intermediate HCW waypoints for now.

cfg.N_homing = [1 1];

% Legacy internal name for the SECOND homing arc.
cfg.N_closing = [1 1];

%% 8. RELATIVE VELOCITY LIMITS

% Maximum total LVLH relative-speed magnitude.
% Preliminary literature-informed screening limits.

cfg.maxSpeed_mps = 1.0;

% HOMING 1: 2 km -> HP1
cfg.maxSpeed_homing_mps = 1.0;

% HOMING 2: HP1 -> HP2
%
% Legacy internal name: maxSpeed_closing_mps
cfg.maxSpeed_closing_mps = 1.0;

% CLOSING: HP2 -> HP3
%
% Legacy internal name: maxSpeed_controlled_mps
cfg.maxSpeed_controlled_mps = 0.3;

% FINAL APPROACH: HP3 -> 5 m
cfg.maxSpeed_final_mps = 0.1;

%% 9. PROPULSION CONSTRAINTS

% Deferred until propulsion preliminary sizing.

cfg.maxPeakAccel_mps2 = inf;
cfg.maxDV_perCycle_mps = inf;

%% 10. HOLDING-POINT DELTA V

% Idealized V-bar stationary holds.
% Additional propulsive allocations currently zero.

cfg.hold1_DVrate_mps_per_min = 0;
cfg.hold2_DVrate_mps_per_min = 0;
cfg.hold3_DVrate_mps_per_min = 0;

%% 11. CLIENT AND ATLANTE OVERRIDES

% No differences assumed yet.
% Can be extended after capture/docking design matures.

cfg.client = struct();
cfg.host = struct();

%% 12. MODEL LIMITATIONS

% - Linearized HCW dynamics
% - Ideal impulsive Homing 1 and Homing 2
% - Continuous quintic closing and final approach
% - No finite-thrust feasibility validation
% - Preliminary authorized KOZ entry corridor
% - No detailed target collision geometry
% - No detailed sensor or navigation model
% - Sampling, not trajectory optimization
% - P95 describes sampled feasible designs
%
% Future:
% - Distance-dependent closing-rate gates
% - Navigation and sensor constraints
% - Propulsion acceleration/thrust constraints
% - Abort and passive safety analysis
% - Objective-based trajectory optimization

end
