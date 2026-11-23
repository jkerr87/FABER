
function cfg = ASV_RPO_MC_config()
% FABER TAXI - NINE-PHASE RPO MONTE CARLO
%
% 1  Homing 1    -> HP1
% 2  Holding 1
% 3  Homing 2    -> HP2
% 4  Holding 2
% 5  Closing 1   -> HP3, authorized AZ entry
% 6  Holding 3   -> KOZ GO decision
% 7  Closing 2   -> HP4, authorized KOZ entry
% 8  HP4 gate    -> capture-readiness GO decision
% 9  Final       -> target standoff
%
% All decisions GO.
% LVLH states [x;y;z;vx;vy;vz].
% Distances m, speeds m/s, durations min.

%% 1. MONTE CARLO

cfg.N = 100;                    % Fixed number of scenarios
cfg.seed = 42;
cfg.nTrajectorySamples = 301;

cfg.outputFolder = 'FABER_RPO_9phase_results';
cfg.makePlots = true;
cfg.saveCSV = true;

%% 2. METHOD OPTIONS - NEW

cfg.methodsH2 = {'GeneralHCW','Vbar','Rbar'};

cfg.methodsC1 = {'GeneralHCW','Vbar','Rbar'};

cfg.methodsC2 = {'GeneralHCW','Continuous'};

% 3 x 3 x 2 x 2 gate modes = 36 combinations.
% Remove options from these arrays to simplify the study.

% Discrete selection of the already generated candidates.
% 'dv', 'tof', 'weighted', or 'pareto' (all nondominated candidates).
cfg.objective = 'pareto';
cfg.weightDV = 0.5;          % used only when objective='weighted'
cfg.weightTOF = 0.5;
cfg.paretoMinAcceptance = 0.0;  % optional architecture-level filter


%% 3. REFERENCE ORBIT

cfg.mu_km3s2 = 398600.4418;
cfg.Re_km = 6378.1363;

% CHANGED: set reference altitude for this trade study.
% Use 457.25 to reproduce the older orbit assumption.
cfg.h_km = 680.12;
cfg.h_km = 457.25;

%% 4. INITIAL RELATIVE STATES - UNCHANGED

cfg.x0_client = [50;-2000;20;0.01;0.04;-0.01];
cfg.x0_host   = [50;-2000;20;0.01;0.04;-0.01];

cfg.dx0_client = [100;300;100;0.05;0.05;0.05];
cfg.dx0_host   = [100;300;100;0.05;0.05;0.05];

%% 4A. INITIAL ARRIVAL DISTRIBUTION
% 'alongtrack' reproduces the original rectangular perturbation.
% 'shell' samples directions uniformly on a 2 km spherical shell.
% 'volume' samples positions uniformly in a spherical annulus.
cfg.initialPositionMode = 'alongtrack';
cfg.shellRadius_m = [2000 2000];    % shell or volume radii [min max]
cfg.shellVelocityMode = 'perturbed'; % 'perturbed' or 'ball'
cfg.shellMaxRelSpeed_mps = 0.05;    % for 'ball': magnitude bound

%% 5. HOLDING-POINT GEOMETRY - CHANGED

cfg.hp1_m = [800 1200];  % Outside AZ
cfg.hp2_m = [350 450];   % Outside AZ

% NEW: intermediate point inside AZ, outside KOZ
cfg.hp3_m = [100 250];

% NEW: terminal decision gate inside KOZ
cfg.hp4_m = [20 20];

cfg.final_m = [5 5];


%% HP4 GATE COMPARISON - NEW

% FullStop:
%   Closing 2 ends at zero relative velocity.
%   Final approach starts from rest.
%
% FlyThrough:
%   Closing 2 ends with positive along-track velocity.
%   Final approach starts with the same velocity.
%
% Both modes initially have zero dwell at HP4.

cfg.gateModes = {'FullStop','FlyThrough'};

% Exploratory approach velocity at HP4 [m/s].
% Positive y points toward the target from negative V-bar.
cfg.gateSpeed_mps = 0.01;

cfg.hold4_min = [0 0];

%% EVENT DETECTION - NEW

% Number of points used to bracket sphere crossings.
% Each crossing time is subsequently refined by fzero.
cfg.eventScanN = 601;

% Numerical tolerance for interpreting crossing direction
cfg.eventRangeRateTol_mps = 1e-8;

% AZ and KOZ are preliminary spherical boundaries.
% No lateral corridor restriction is imposed.


%% 6. SAFETY ZONES

cfg.AZ_m = 300;  % Preliminary: requires validation
cfg.KOZ_m = 30;  % RPO lecture assumption

% NEW: numerical tolerance ONLY for zone comparisons
cfg.zoneTolerance_m = 1e-6;

% REMOVED:
% cfg.corridorHalfWidth_m
% cfg.maxClosingExcursionFactor
%
% No lateral corridor exclusion is applied.

%% 7. PHASE TIME BOUNDS [min] - CHANGED

cfg.homing1_min = [30 130];
cfg.hold1_min   = [10 60];

cfg.homing2_min = [20 130];
cfg.hold2_min   = [5 45];

% Closing 1: HP2 -> HP3, crosses AZ.
%
% 110 min upper bound is a TRADE-STUDY relaxation
% permitting the one-orbit analytical V-bar transfer.
cfg.closing1_min = [10 110];

% HP3: short hold and KOZ GO decision
cfg.hold3_min = [2 10];

% Closing 2: HP3 -> HP4, crosses KOZ
cfg.closing2_min = [5 60];

% HP4: zero-velocity gate, no dwell by default
cfg.hold4_min = [0 0];

% Final approach: HP4 -> 5 m
cfg.final_min = [5 30];

cfg.minTransferTime_s = 120;

%% 8. PHASE SPEED LIMITS - CHANGED

% Total LVLH relative-speed magnitude.

cfg.maxSpeed_homing1_mps = 1.0;
cfg.maxSpeed_homing2_mps = 1.0;

cfg.maxSpeed_closing1_mps = 0.3;
cfg.maxSpeed_closing2_mps = 0.3;

cfg.maxSpeed_final_mps = 0.1;

%% 9. PROPULSION CONSTRAINTS



%% TAXI LOADED / UNLOADED MASS

% Initial Taxi wet mass before client retrieval [kg].
% This is the reference mission mass before module capture.
cfg.massTaxi_kg = 3503.59;

% Captured payload mass [kg].
% REPLACE with actual transported module + container mass.
cfg.payloadMass_kg = 0;

% Enable mass change after client capture.
cfg.enablePayloadMassChange = true;

% Preliminary propulsion inputs; replace with qualified subsystem values.
cfg.thrust_N = 22;              % effective commanded translational thrust
cfg.Isp_s = 293;               % ILLUSTRATIVE assumption - verify thruster Isp
cfg.nActiveThrusters = 4;       % usable co-aligned thrusters; set actual count
cfg.thrustEfficiency = 1.0;      % factor for usable force, 0 < eta <= 1
% Optional detailed LVLH thruster arrangement (3xN direction columns).
% If provided, thrusterForces_N must contain the matching N rated thrusts.
% Requires Optimization Toolbox (linprog). This assumes directions are
% already expressed in LVLH, not body frame, and ignores attitude torques.
cfg.thrusterDirections_LVLH = [];
cfg.thrusterForces_N = [];
cfg.nThrustSamples = 201;       % directional feasibility for continuous control
cfg.brakingMargin_m = 0;        % extra clearance in radial stopping screen
cfg.maxPropellantCycle_kg = inf; % set finite allocation to enforce fuel cap
cfg.maxBurnDuration_s = inf;   % optional per-impulse burn time limit
cfg.maxBurnFraction = 0.05;    % cumulative impulse durations / phase TOF
cfg.maxPeakAccel_mps2 = inf;   % optional additional limit
cfg.maxDV_perCycle_mps = inf;

% Closing-speed check: 'total', 'radial', or 'both'.
% Radial uses rho-dot = dot(r,v)/norm(r), inward speed = max(0,-rho-dot).
cfg.closingSpeedMode = 'both';

% Braking-distance screening (not a finite-burn propagation).
cfg.checkBrakingDistance = true;
cfg.brakingSafetyFactor = 1.0;


%% PHASE-SPECIFIC TAXI PROPULSION

% Master switch
% false = disable propulsion, fuel, braking and burn-fraction checks
cfg.enablePropulsionConstraints = true;

% 'phase' = phase-specific main/RCS systems
% 'legacy' = existing cfg.thrust_N and cfg.Isp_s
cfg.propulsionMode = 'phase';

% TAXI MAIN - L3Harris R-4D-11
cfg.propulsion.main.thrust_N = 490;
cfg.propulsion.main.nInstalled = 2;
cfg.propulsion.main.nActiveThrusters = 1;
cfg.propulsion.main.Isp_s = 315.5;
cfg.propulsion.main.thrustEfficiency = 1.0;
cfg.propulsion.main.thrusterDirections_LVLH = [];
cfg.propulsion.main.thrusterForces_N = [];

% TAXI RCS - AMPAC-ISP
cfg.propulsion.rcs.thrust_N = 22;
cfg.propulsion.rcs.nInstalled = 20;
cfg.propulsion.rcs.nActiveThrusters = 4;
cfg.propulsion.rcs.Isp_s = 293.0;
cfg.propulsion.rcs.thrustEfficiency = 1.0;
cfg.propulsion.rcs.thrusterDirections_LVLH = [];
cfg.propulsion.rcs.thrusterForces_N = [];

% Nine-phase propulsion assignment
cfg.phasePropulsion = { ...
    'main', ... % 1 Homing 1
    'rcs',  ... % 2 Holding 1
    'rcs',  ... % 3 Homing 2
    'rcs',  ... % 4 Holding 2
    'rcs',  ... % 5 Closing 1
    'rcs',  ... % 6 Holding 3
    'rcs',  ... % 7 Closing 2
    'rcs',  ... % 8 HP4 gate
    'rcs'};     % 9 Final approach


%% 9A. FINITE-THRUST VALIDATION (POST ANALYTICAL MONTE CARLO)
% 'none'     : skip finite-thrust propagation
% 'selected' : scenario-selected solutions, P95 representatives, extreme
%              burn-fraction cases (debugging / preliminary verification)
% 'pareto'   : all scenario-level nondominated analytically feasible options
% 'all'      : all analytically feasible options
cfg.finiteThrustValidation = 'selected';
cfg.finiteMaxCandidates = 100; % Inf for unlimited; subset IS NOT full validation
cfg.finiteOutputSamples = 101;
cfg.finiteRelTol = 1e-8;
cfg.finiteAbsTol = 1e-9;
% Finite trajectory endpoint tolerances, at each holding/final waypoint:
cfg.finitePosTolerance_m = 5;    % preliminary, must justify by mission phase
cfg.finiteVelTolerance_mps = 0.02;
cfg.finiteFinalPosTolerance_m = 0.5;
cfg.finiteFinalVelTolerance_mps = 0.005;
% All idealized finite burn commands are applied in the LVLH direction.
% If body-frame orientation is known, thruster directions should be
% transformed to LVLH as a function of time (NOT implemented here).
% Open-loop pulse timings, without feedback, are checked for drift.
% Each leg starts from the nominal planned waypoint, not the propagated
% previous-leg error. Full chained guidance/hold control not included.

%% 10. HOLDING POINT DELTA V - CHANGED

% Four stationary reference positions.
% Ideal holding cost = 0 until stationkeeping is sized.
%
% HP4 has no dwell by default.
cfg.holdDVrate_mps_per_min = [0 0 0 0];

%% 11. DECISION GATES - NEW

% All GO for this analysis.
cfg.GO_HP1 = true;
cfg.GO_HP2 = true;
cfg.GO_AZ = true;
cfg.GO_KOZ = true;
cfg.GO_capture = true;


%% INITIAL TAXI DEPARTURE FROM ATLANTE

cfg.departure.enabled = true;

% Preliminary geometry: change to match CONOPS
cfg.departure.initialState = ...
    [0;-5;0;0;0;0];

cfg.departure.finalState = ...
    [0;-350;0;0;0;0];

% One initial departure maneuver
cfg.departure.TOF_min = 60;

% Use RCS for preliminary departure
cfg.departure.propulsion = 'rcs';

% Feasibility assumptions
cfg.departure.maxSpeed_mps = 1.0;
cfg.departure.minSeparation_m = 5;

cfg.departure.posTolerance_m = 1;
cfg.departure.velTolerance_mps = 0.01;


%% 12. CONOPS LIMITATIONS

% - Linear HCW circular-orbit dynamics
% - Impulsive HCW / analytical V-bar / analytical R-bar
% - Ideal continuous quintic reference tracking
% - Zero-relative-velocity holding-point endpoints
% - No 5-m lateral corridor exclusion
% - Preliminary AZ and KOZ radii
% - Sampled trajectory feasibility checks
% - All decision gates GO
% - Impulse durations/fuel are screened, but finite-burn trajectories are NOT propagated
% - No passive-safety certification
% - Discrete Pareto screening only; no continuous optimal-control solve
% - P95 conditional on accepted sampled trajectories

end
