% Working script for inspecting and plotting FABER analysis results.
clear; clc; close all;

% out = MA_vehicle_budget_phaseA();
%
% opts = struct();
% opts.savePNG = true;
% opts.makeZoom = false;  % Disable old side-by-side zoom
%
% plots = plot_FABER_P95_RPO_3D(out,opts);
%
% B = print_FABER_budget(out);

% Run ONLY RPO Monte Carlo method comparison; mission budget unaffected.
% opts=struct();
% opts.N=60;                   % Paired candidate draws; raise after initial check
% opts.seed=42;
% opts.h_km=680.12;          % Match current FABER mission target orbit
% opts.expandClosingWindow=true; % Required to allow P and P/2 analytical cases
% opts.runLegacyBaseline=true;
% opts.nSamples=51;           % Raise to 201+ for sampled safety sensitivity
% opts.makePlots=true;         % Exactly 2 method-overlay figures
% opts.saveCSV=true;
% R=FABER_RPO_compare_methods(opts);

% clear; clc; close all;
%
% cfg = ASV_RPO_MC_config();
%
% % Reference orbit
% cfg.h_km = 680.12;
%
% % Monte Carlo
% cfg.N = 250;
% cfg.seed = 42;
% cfg.nTrajectorySamples = 101;
%
% % Enable comparison
% cfg.compareMethods = true;
%
% % Extend the closing window to permit
% % the analytical V-bar maneuver.
% cfg.comparisonExpandClosing = true;
%
% % Output
% cfg.makePlots = true;
% cfg.saveCSV = true;
% cfg.outputFolder = 'FABER_RPO_integrated_results';
%
% % Run
% R = ASV_RPO_MonteCarlo_V2(cfg);
%
% % Results
% disp(R.summary);

%% FABER RPO - THREE-METHOD MONTE CARLO COMPARISON
% Methods:
%   1. General HCW
%   2. Analytical V-bar
%   3. Analytical R-bar
%
% Hard constraints:
%   - KOZ approach corridor
%   - Non-increasing range after KOZ entry
%   - Maximum relative speed
%   - Transfer time limits
%
% Final approach remains continuous quintic HCW.
%
% Requires the updated ASV_RPO_MonteCarlo.m
% and ASV_RPO_MC_config.m
%
% clear; clc; close all;
%
% %% 1. INITIAL CONFIGURATION
%
% cfg = ASV_RPO_MC_config();
%
% cfg.h_km = 680.12;
%
% % Monte Carlo settings
% cfg.N = 250;
% cfg.seed = 42;
% cfg.maxAttemptsFactor = 30;
%
% % Number of samples along each trajectory
% cfg.nTrajectorySamples = 501;
%
% %% 2. ENABLE THREE-METHOD COMPARISON
%
% cfg.compareMethods = true;
%
% % General HCW
% % Analytical V-bar
% % Analytical R-bar
%
% %% 3. RPO GEOMETRY
%
% cfg.AZ_m = 300;
% cfg.KOZ_m = 30;
%
% % Maximum radial/cross-track deviation inside KOZ
% cfg.corridorHalfWidth_m = 5;
%
% %% 4. HARD KOZ CONSTRAINTS
%
% % Reject trajectories violating the V-bar corridor
% cfg.enforceKOZCorridor = true;
%
% % Reject outward range excursions after entering KOZ
% cfg.enforceKOZMonotonic = true;
%
% % Numerical tolerance for sampled range changes
% cfg.maxOutwardRangeStep_m = 1e-3;
%
% %% 5. HARD RELATIVE SPEED CONSTRAINTS
%
% % Homing 1
% cfg.maxSpeed_homing_mps = 1.0;
%
% % Homing 2
% cfg.maxSpeed_closing_mps = 0.3;
%
% % Closing: HP2 to HP3
% cfg.maxSpeed_controlled_mps = 0.3;
%
% % Continuous final approach
% cfg.maxSpeed_final_mps = 0.1;
%
% %% 6. TRANSFER TIME CONSTRAINTS
%
% % Homing 1
% cfg.homing_min = [30 100];
%
% % Holding point 1
% cfg.hold1_min = [10 60];
%
% % Homing 2
% cfg.closing_min = [20 130];
%
% % Holding point 2
% cfg.hold2_min = [5 45];
%
% % Closing: HP2 to HP3
% cfg.controlled_min = [10 70];
%
% % Holding point 3
% cfg.hold3_min = [5 15];
%
% % Continuous final approach
% cfg.final_min = [5 30];
%
% %% 7. CHOOSE CLOSING TIME SCENARIO
%
% % OPTION A: STRICT ORIGINAL CONSTRAINTS
% % V-bar analytical transfer (~98 min) will be infeasible.
%
% cfg.comparisonExpandClosing = false;
% cfg.enforceOriginalClosingTime = true;
%
% % OPTION B: RELAXED TIMING STUDY
% % Uncomment these two lines instead to permit
% % the analytical V-bar transfer.
% %
% % cfg.comparisonExpandClosing = true;
% % cfg.enforceOriginalClosingTime = false;
%
% %% 8. OTHER FEASIBILITY LIMITS
%
% cfg.minTransferTime_s = 120;
%
% % No upper Delta-V rejection threshold for this trade study
% cfg.maxDV_perCycle_mps = inf;
%
% % No finite-thrust acceleration limit until propulsion sizing
% cfg.maxPeakAccel_mps2 = inf;
%
% %% 9. OUTPUT SETTINGS
%
% cfg.makePlots = true;
% cfg.saveCSV = true;
%
% cfg.outputFolder = 'FABER_RPO_strict_results';
%
% %% 10. RUN MONTE CARLO
%
% R = ASV_RPO_MonteCarlo_V2(cfg);
%
% %% 11. SUMMARY RESULTS
%
% fprintf('\n========================================\n');
% fprintf('FABER THREE-METHOD RPO COMPARISON\n');
% fprintf('========================================\n');
%
% disp(R.summary);
%
% %% 12. FEASIBILITY AND REJECTION REASONS
%
% for m = 1:numel(R.methods)
%
%     methodName = R.methods{m};
%
%     nAccepted = sum(R.valid(:,m));
%     nTotal = size(R.valid,1);
%
%     fprintf('\n%s\n',methodName);
%     fprintf('Accepted: %d / %d (%.1f%%)\n', ...
%         nAccepted,nTotal, ...
%         100*nAccepted/nTotal);
%
%     reasons = R.rejectionReason(:,m);
%     reasons = reasons(strlength(reasons)>0);
%
%     if ~isempty(reasons)
%
%         reasonTypes = categorical(reasons);
%
%         fprintf('\nRejection reasons:\n');
%
%         disp(table( ...
%             categories(reasonTypes), ...
%             countcats(reasonTypes), ...
%             'VariableNames',{'Reason','Count'}));
%
%     end
%
% end
%
% %% 13. COMMON FEASIBLE CANDIDATES
%
% commonFeasible = all(R.valid,2);
%
% fprintf('\n========================================\n');
% fprintf('COMMON FEASIBLE CANDIDATES\n');
% fprintf('========================================\n');
%
% fprintf('Feasible for all three methods: %d / %d\n', ...
%     sum(commonFeasible),size(R.valid,1));
%
% %% 14. COMPLETE
%
% fprintf('\nMonte Carlo comparison complete.\n');
%
%

% Load the assembly configuration
run('truss_platform.m');

% Orbital parameters
mu = 398600.4418;  % km^3/s^2
Re = 6378.1363;    % km
h = 680.12;        % km

% Assembly geometry
hostPoint = [0;-250;0];

% Candidate HCW transfer durations
timesMin = 5:5:250;

% Trajectory resolution
Npts = 201;

% Baseline options
opts = struct();

opts.methods = {'GeneralHCW','Vbar','Rbar'};

opts.gateDistance_m = 20;

opts.exitTimesMin = 5:5:30;
opts.finalTimesMin = 5:5:30;

opts.maxExitSpeed_mps = 0.10;
opts.maxTransferSpeed_mps = 0.30;
opts.maxFinalSpeed_mps = 0.10;

opts.maxExitAccel_mps2 = Inf;
opts.maxFinalAccel_mps2 = Inf;

opts.printSummary = false;

% Assembly relocation optimization
opts.objective = 'weighted';

% Equal normalized importance of DV and TOF
opts.weightDV = 0.5;
opts.weightTOF = 0.5;

% Reference scales for objective normalization
opts.dvReference_mps = 0.1;
opts.tofReference_min = 60;

% Run independent sensitivity analysis
S = assembly_relocations_sensitivity_OAT( ...
    mu,Re,h,Npts,timesMin, ...
    hostPoint,modulePoints,opts);
