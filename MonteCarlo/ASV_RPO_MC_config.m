
function cfg = ASV_RPO_MC_config()
% ASV_RPO_MC_CONFIG
% FABER Taxi - Proximity Operations Monte Carlo
%
% Preliminary literature-informed design bounds.
% These are NOT validated mission requirements.
%
% Reference missions:
% - ESA ATV rendezvous operations
% - Martinez et al. (2023), autonomous RPO GNC
% - Hinkel et al. (2014), NASA rendezvous strategy
%
% All sampled design variables follow uniform distributions.

%% ================================================================
% 1. MONTE CARLO SETTINGS
% ================================================================

% Number of accepted client + host retrieval cycles
% Start small for validation; increase after convergence testing.
cfg.N = 100;

% Reproducibility
cfg.seed = 42;

% Maximum attempts = N * maxAttemptsFactor
cfg.maxAttemptsFactor = 3;

% Samples used to check each HCW trajectory
cfg.nTrajectorySamples = 25;

% Results
cfg.outputFolder = 'ASV_RPO_MC_results';
cfg.makePlots = true;
cfg.saveCSV = true;


%% ================================================================
% 2. OPTIMIZATION SETTINGS
% ================================================================

% Optimization objective:
% 'dv'       = minimum feasible Delta V
% 'tof'      = minimum feasible time of flight
% 'weighted' = weighted Delta V / TOF compromise

cfg.objective = 'dv';

% Available closing / final approach geometries
cfg.methods = {'Vbar'};

% Random candidate designs evaluated per approach
cfg.candidatesPerLeg = 3;

% Weighted objective parameters
cfg.weightDV = 0.5;
cfg.weightTOF = 0.5;

% Normalization values (illustrative)
cfg.refDV_mps = 2.0;
cfg.refTOF_min = 400;


%% ================================================================
% 3. ORBIT PARAMETERS
% ================================================================

% Can be overwritten by MA_deltav_v2.m

cfg.mu_km3s2 = 398600.4418;   % [km^3/s^2]
cfg.Re_km = 6378.1363;       % [km]
cfg.h_km = 457.25;           % [km]


%% ================================================================
% 4. INITIAL RELATIVE STATES
% ================================================================

% LVLH state:
% [x;y;z;vx;vy;vz]
%
% x = radial outward
% y = along-track
% z = cross-track
%
% Position [m], velocity [m/s]

cfg.x0_client = [50;-2000;20;0.01;0.04;-0.01];

cfg.x0_host = [50;-2000;20;0.01;0.04;-0.01];


%% ================================================================
% 5. MANEUVER GEOMETRY [m]
% ================================================================

% First holding point
% Literature reference: SENER 1000 m waypoint
cfg.hp1_m = [800 1200];

% Second holding point
% References: NASA 100 m, SENER R-bar 200 m
%
% Current code requires HP2 inside the 300 m AZ.
cfg.hp2_m = [100 250];

% Final robotic capture standoff
% FABER-specific preliminary assumption
cfg.final_m = [5 5];

% Approach Zone boundary
cfg.AZ_m = 300;

% Keep-Out Zone radius
% Preliminary assumption; refine using target dimensions,
% navigation uncertainties and collision-clearance requirements.
cfg.KOZ_m = 30;

% Final approach corridor half-width
% Preliminary assumption
cfg.corridorHalfWidth_m = 5;

cfg.maxClosingExcursionFactor = 1.50;


%% ================================================================
% 6. MANEUVER TIME BOUNDS [min]
% ================================================================

% Literature-informed exploratory bounds.
% References include ESA ATV operational timelines.
%
% Actual minimum/maximum values are design assumptions,
% not directly published operational limits.

% HOMING
cfg.homing_min = [30 100];

% HOLDING POINT 1
cfg.hold1_min = [10 60];

% CLOSING
cfg.closing_min = [20 130];

% HOLDING POINT 2
cfg.hold2_min = [5 45];

% FINAL APPROACH
cfg.final_min = [5 30];


%% ================================================================
% 7. NUMBER OF MANEUVERS
% ================================================================

% Current Monte Carlo implementation supports:
% - One homing transfer
% - One closing transfer
% - One final approach
%
% ATV multi-maneuver architecture is reserved for future work.
%
% These variables are descriptive only and are NOT currently
% used by ASV_RPO_MonteCarlo.m.

cfg.N_homing = [1 4];
cfg.N_closing = [1 4];


%% ================================================================
% 8. SAFETY AND FEASIBILITY CONSTRAINTS
% ================================================================

% Global relative-speed limit
% Preliminary sensitivity case inspired by NASA speed gates.
%
% NOTE: The published reference uses distance-dependent limits.
% The existing code applies this one value to every trajectory.

cfg.maxSpeed_mps = 1.0;       % [m/s]

% Distance/phase-dependent speed limits [m/s]
cfg.maxSpeed_homing_mps = 1.0;
cfg.maxSpeed_closing_mps = 1.0;
cfg.maxSpeed_final_mps = 1.0;

% Maximum acceleration for continuous final approach
% Unlimited until Taxi thruster performance is available.
cfg.maxPeakAccel_mps2 = inf;

% Maximum Delta V for a complete client + host RPO cycle
% Unlimited until propulsion allocation is established.
cfg.maxDV_perCycle_mps = inf;


%% ================================================================
% 9. HOLDING-POINT CONTROL ALLOCATIONS
% ================================================================

% Additional propulsive allocation per minute of holding.
%
% V-bar stationary holds are idealized as zero-DV.
% R-bar holding control is calculated by the MC function.

cfg.hold1_DVrate_mps_per_min = 0;

cfg.hold2_DVrate_mps_per_min = 0;


%% ================================================================
% 10. MODEL LIMITATIONS
% ================================================================

% - HCW linearized relative-motion dynamics
% - Impulsive homing and closing transfers
% - Ideal continuous final approach
% - V-bar / R-bar option changes closing target geometry
%   and final approach axis
% - Closing is not strictly axis-constrained
% - KOZ screening uses sampled trajectory points
% - No detailed sensor model
% - No navigation uncertainty
% - No finite-thrust model for impulsive transfers
% - No automatic optimization of maneuver count
% - P95 describes feasible sampled designs, not reliability
%
% Future refinements:
% - Distance-dependent relative-speed limits
% - Separate V-bar and R-bar waypoint ranges
% - Multiple homing and closing maneuvers
% - Sensor acquisition and visibility constraints
% - Thruster limits and propellant consumption
% - Navigation uncertainty and robustness analysis

end
