% Run the nine-phase RPO study from the perturbation-analysis package.
clear; clc; close all;

%% Configuration
cfg = ASV_RPO_MC_config();

cfg.N = 1000;
cfg.seed = 42;
cfg.h_km = 680.12;

%% NEW CONFIGURATION OPTIONS
% Nominal baseline vs directionally uniform 2 km arrival shell.
cfg.initialPositionMode = 'alongtrack'; % 'alongtrack' / 'shell' / 'volume'
cfg.shellRadius_m = [2000 2000];
cfg.shellVelocityMode = 'perturbed'; % or 'ball'

% Propulsion preliminary assumptions: verify values with subsystem design.
cfg.thrust_N = 25;         % N per active thruster (or effective thrust)
cfg.nActiveThrusters = 1;  % REPLACE with usable simultaneous thruster count
cfg.thrustEfficiency = 1; % directional/control derating if applicable
cfg.massTaxi_kg = 3503.59;

cfg.Isp_s = 220;                 % PLACEHOLDER, replace with actual Isp
cfg.maxPropellantCycle_kg = inf; % specify a finite cycle fuel allocation
cfg.maxBurnFraction = 0.05;
cfg.maxBurnDuration_s = inf;
cfg.brakingMargin_m = 0;   % set operational clearance margin

% Closing-speed criterion, using rho-dot = dot(r,v)/norm(r).
cfg.closingSpeedMode = 'both'; % 'total' / 'radial' / 'both'
cfg.checkBrakingDistance = true;

% Candidate selection: 'dv', 'tof', 'weighted', 'pareto'.
cfg.objective = 'pareto';
cfg.weightDV = 0.5;
cfg.weightTOF = 0.5;
cfg.paretoMinAcceptance = 0.70;

% Optional second independent MC run; don't mix shell with nominal P95.
runShellSensitivity = false;

%% Run nine-phase comparison
R = FABER_RPO_9phase_MC(cfg);

%% Display all 36 combinations (default gate modes)
disp(R.summary);

%% Selected candidate per initial-state sample
selectedOK = isfinite(R.selectedArchitecture);
fprintf('Selected feasible trials: %d / %d\n',nnz(selectedOK),cfg.N);
if any(selectedOK)
    fprintf('Selected median DV: %.4f m/s | TOF: %.2f min\n', ...
        median(R.selectedDV_mps(selectedOK)), ...
        median(R.selectedTOF_min(selectedOK)));
end

%% Show best accepted architecture by P95 Delta-V
validRows = find(R.summary.Accepted > 0);

if ~isempty(validRows)

    [~,k] = min( ...
        R.summary.P95CycleDV_mps(validRows));

    best = validRows(k);

    fprintf('\nLowest sampled P95 cycle DV:\n');
    disp(R.summary(best,:));

else
    fprintf('\nNo feasible architectures found.\n');
end

%% Display rejection reasons for each architecture
for m = 1:size(R.valid,2)

    fprintf('\nArchitecture %d: %s / %s / %s\n', ...
        m, ...
        char(R.methods(m,1)), ...
        char(R.methods(m,2)), ...
        char(R.methods(m,3)));

    failures = R.rejectionReason(:,m);

    failures = failures(strlength(failures)>0);

    if ~isempty(failures)

        [names,~,idx] = unique(failures);
        counts = accumarray(idx,1);

        disp(table(names,counts, ...
            'VariableNames',{'Reason','Count'}));
    end
end

%% NEW: Plot all architectures in docked figures
plot_FABER_RPO_P95(R);

%% PHASE-BY-PHASE ANALYSIS
% Architecture 5: Vbar / Vbar / GeneralHCW

m = 5;

ok = R.valid(:,m);

if ~any(ok)
    candidates = find(R.summary.Accepted>0);
    if isempty(candidates)
        warning('No feasible architectures for phase analysis.');
        return;
    end

    [~,j] = max(R.summary.Accepted(candidates));
    m = candidates(j);
    ok = R.valid(:,m);
    fprintf('Architecture 5 unavailable, using architecture %d instead.\n',m);
end

% Dimensions of original arrays:
% Attempts x Architectures x Approach x Phase
%
% Approach 1 = Client
% Approach 2 = ATLANTE

DVclient = reshape( ...
    R.DV_phase_mps(ok,m,1,:),[],9);

DVhost = reshape( ...
    R.DV_phase_mps(ok,m,2,:),[],9);

TOFclient = reshape( ...
    R.TOF_phase_min(ok,m,1,:),[],9);

TOFhost = reshape( ...
    R.TOF_phase_min(ok,m,2,:),[],9);

% Sum both approaches trial by trial
DVphase = DVclient + DVhost;
TOFphase = TOFclient + TOFhost;

phaseNames = string(R.phaseNames(:));

% Empirical nearest-rank P95 for each phase
p95DV = zeros(9,1);
p95TOF = zeros(9,1);

for k = 1:9
    x = sort(DVphase(:,k));
    idx = max(1,ceil(0.95*numel(x)));
    p95DV(k) = x(idx);

    x = sort(TOFphase(:,k));
    idx = max(1,ceil(0.95*numel(x)));
    p95TOF(k) = x(idx);
end

results = table( ...
    phaseNames, ...
    mean(DVphase,1)', ...
    p95DV, ...
    mean(TOFphase,1)', ...
    p95TOF, ...
    'VariableNames',{ ...
    'Phase', ...
    'MeanDV_mps', ...
    'P95DV_mps', ...
    'MeanTOF_min', ...
    'P95TOF_min'});

disp(results);

%% Plot phase Delta-V
figure('WindowStyle','docked','Color','w');
barh(p95DV);
set(gca,'YTick',1:9, ...
    'YTickLabel',phaseNames, ...
    'YDir','reverse');
xlabel('Phase P95 Delta-V [m/s]');
title('FABER - Architecture 5 phase contributions');
grid on;

%% Compare all feasible architectures
P = FABER_RPO_Pareto(R,cfg.paretoMinAcceptance);

%% Optional directional sensitivity on the 2 km shell
if runShellSensitivity
    shellCfg = cfg;
    shellCfg.initialPositionMode = 'shell';
    shellCfg.outputFolder = 'FABER_RPO_9phase_shell_results';
    shellCfg.makePlots = false;
    shellCfg.saveCSV = true;
    R_shell = FABER_RPO_9phase_MC(shellCfg);
    disp(R_shell.summary);
    if any(R_shell.summary.Accepted>0)
        P_shell = FABER_RPO_Pareto(R_shell,cfg.paretoMinAcceptance);
    end

    fprintf('Nominal feasible scenarios: %d/%d\n', ...
        nnz(any(R.valid,2)),cfg.N);
    fprintf('Shell feasible scenarios:   %d/%d\n', ...
        nnz(any(R_shell.valid,2)),shellCfg.N);
end

%% OPTIONAL: FINITE-THRUST PROPAGATION OF SCREENED CANDIDATES
% Executes after analytical MC + discrete Pareto candidate selection.
% Use 'selected' to debug; 'pareto' for full frontier validation.
% 'all' validates all analytically feasible candidates.
if ~strcmpi(cfg.finiteThrustValidation,'none')
    V = FABER_RPO_finite_validation(R,cfg.finiteThrustValidation);
    fprintf('Finite validation: %d / %d passed (of %d analytically feasible).\n', ...
        V.passedCount,V.attemptedCount,V.nominalValidCount);
    if V.passedCount>0
        plot_FABER_RPO_finite(R,V);
    end
end
