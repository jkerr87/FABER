
%% FABER NINE-PHASE RPO MONTE CARLO
% Phase-specific propulsion / legacy / unconstrained comparison
% Discrete Pareto compromise: Delta-V vs TOF

clear; clc; close all;

%% =========================================================
% 1. USER SETTINGS
% ==========================================================

cfg = ASV_RPO_MC_config();

cfg.N = 1000;
cfg.seed = 42;
cfg.h_km = 680.12;

% Initial-condition sampling
cfg.initialPositionMode = 'alongtrack';
cfg.shellRadius_m = [2000 2000];
cfg.shellVelocityMode = 'perturbed';

%% PROPULSION STUDY

% Choose one:
% 'phase'  = actual main/RCS propulsion per phase
% 'legacy' = original 25 N and 220 s Isp
% 'off'    = no propulsion screening

propulsionStudy = 'phase';

switch lower(propulsionStudy)

    case 'phase'
        cfg.enablePropulsionConstraints = true;
        cfg.propulsionMode = 'phase';

    case 'legacy'
        cfg.enablePropulsionConstraints = true;
        cfg.propulsionMode = 'legacy';

        cfg.thrust_N = 25;
        cfg.nActiveThrusters = 1;
        cfg.thrustEfficiency = 1.0;
        cfg.Isp_s = 220;

    case 'off'
        cfg.enablePropulsionConstraints = false;
        cfg.propulsionMode = 'phase';

    otherwise
        error('Unknown propulsion study mode.');

end

% Vehicle mass
cfg.massTaxi_kg = 3503.59;

% Propulsion screening
cfg.maxPropellantCycle_kg = inf;
cfg.maxBurnFraction = 0.05;
cfg.maxBurnDuration_s = inf;

% Closing speed and braking
cfg.closingSpeedMode = 'both';
cfg.checkBrakingDistance = true;
cfg.brakingMargin_m = 0;

%% PARETO CONFIGURATION

cfg.objective = 'pareto';
cfg.weightDV = 0.5;
cfg.weightTOF = 0.5;
cfg.paretoMinAcceptance = 0;

%% OPTIONAL ANALYSES

runShellSensitivity = false;

% 'none', 'selected', 'pareto', 'all'
cfg.finiteThrustValidation = 'none';
cfg.finiteMaxCandidates = 30;

% Finite thrust requires propulsion to be enabled
if ~cfg.enablePropulsionConstraints
    cfg.finiteThrustValidation = 'none';
end

% Plot settings
showPareto = true;
showSelectedScatter = true;
showPhaseResults = true;

%% OUTPUT SETTINGS

cfg.outputFolder = fullfile( ...
    fileparts(mfilename('fullpath')), ...
    ['FABER_RPO_results_' propulsionStudy]);

if ~isfolder(cfg.outputFolder)
    mkdir(cfg.outputFolder);
end

cfg.makePlots = false;
cfg.saveCSV = false;

%% =========================================================
% 2. RUN MONTE CARLO
% ==========================================================

fprintf('\n=======================================\n');
fprintf('FABER NINE-PHASE RPO MONTE CARLO\n');
fprintf('=======================================\n');

fprintf('Propulsion mode: %s\n',propulsionStudy);
fprintf('Trials: %d\n',cfg.N);

R = FABER_RPO_9phase_MC(cfg);

%% =========================================================
% 3. FEASIBILITY
% ==========================================================

feasible = isfinite(R.selectedArchitecture);

fprintf('\n=======================================\n');
fprintf('FEASIBILITY RESULTS\n');
fprintf('=======================================\n');

fprintf('Feasible scenarios: %d / %d\n', ...
    nnz(feasible),cfg.N);

fprintf('Scenario acceptance: %.2f %%\n', ...
    100*nnz(feasible)/cfg.N);

%% =========================================================
% 4. PARETO COMPROMISE SELECTION
% ==========================================================

if ~any(feasible)
    error('No feasible scenarios. Review constraints.');
end

trialIdx = find(feasible);
archIdx = R.selectedArchitecture(feasible);

DV = R.selectedDV_mps(feasible);
TOF = R.selectedTOF_min(feasible);

fuel = nan(numel(trialIdx),1);

for j = 1:numel(trialIdx)
    fuel(j) = R.FuelCycle_kg( ...
        trialIdx(j),archIdx(j));
end

% Nearest-rank statistics
Nfeas = numel(DV);
idx95 = max(1,ceil(0.95*Nfeas));

sortedDV = sort(DV);
sortedTOF = sort(TOF);
sortedFuel = sort(fuel);

fprintf('\n=======================================\n');
fprintf('SELECTED PARETO COMPROMISE\n');
fprintf('=======================================\n');

fprintf('Mean DV:  %.4f m/s\n',mean(DV));
fprintf('P50 DV:   %.4f m/s\n',median(DV));
fprintf('P95 DV:   %.4f m/s\n',sortedDV(idx95));

fprintf('Mean TOF: %.2f min\n',mean(TOF));
fprintf('P50 TOF:  %.2f min\n',median(TOF));
fprintf('P95 TOF:  %.2f min\n',sortedTOF(idx95));

fprintf('P95 fuel: %.4f kg\n',sortedFuel(idx95));

%% =========================================================
% 5. SELECTED ARCHITECTURE FREQUENCIES
% ==========================================================

M = size(R.valid,2);

counts = zeros(M,1);

for m = 1:M
    counts(m) = nnz(R.selectedArchitecture == m);
end

ArchitectureResults = R.summary;

ArchitectureResults.SelectionCount = counts;
ArchitectureResults.SelectionFraction = ...
    counts/Nfeas;

ArchitectureResults = sortrows( ...
    ArchitectureResults,'SelectionCount','descend');

fprintf('\nMOST FREQUENT COMPROMISE ARCHITECTURES\n');

disp(ArchitectureResults( ...
    1:min(10,height(ArchitectureResults)),:));

%% =========================================================
% 6. PHASE STATISTICS OF SELECTED COMPROMISES
% ==========================================================

nPhases = numel(R.phaseNames);

phaseDV = nan(Nfeas,nPhases);
phaseTOF = nan(Nfeas,nPhases);

for j = 1:Nfeas

    i = trialIdx(j);
    m = archIdx(j);

    dv1 = reshape(R.DV_phase_mps(i,m,1,:),1,nPhases);
    dv2 = reshape(R.DV_phase_mps(i,m,2,:),1,nPhases);

    t1 = reshape(R.TOF_phase_min(i,m,1,:),1,nPhases);
    t2 = reshape(R.TOF_phase_min(i,m,2,:),1,nPhases);

    phaseDV(j,:) = dv1+dv2;
    phaseTOF(j,:) = t1+t2;

end

phaseP95DV = nan(nPhases,1);
phaseP95TOF = nan(nPhases,1);

for k = 1:nPhases
    sDV = sort(phaseDV(:,k));
    sTOF = sort(phaseTOF(:,k));

    phaseP95DV(k) = sDV(idx95);
    phaseP95TOF(k) = sTOF(idx95);
end

PhaseResults = table( ...
    string(R.phaseNames(:)), ...
    mean(phaseDV,1)', ...
    phaseP95DV, ...
    mean(phaseTOF,1)', ...
    phaseP95TOF, ...
    'VariableNames',{ ...
    'Phase','MeanDV_mps','P95DV_mps', ...
    'MeanTOF_min','P95TOF_min'});

if showPhaseResults
    disp(PhaseResults);

    figure('Color','w','WindowStyle','docked');
    barh(phaseP95DV);
    set(gca,'YTick',1:nPhases, ...
        'YTickLabel',R.phaseNames,'YDir','reverse');
    xlabel('P95 Delta-V [m/s]');
    title('Selected Pareto Compromises - Phase Delta-V');
    grid on;
end

%% =========================================================
% 7. PARETO PLOTS
% ==========================================================

if showPareto
    P = FABER_RPO_Pareto( ...
        R,cfg.paretoMinAcceptance);
end

if showSelectedScatter

    figure('Color','w','WindowStyle','docked');
    hold on; grid on; box on;

    scatter( ...
        R.CycleTOF_min(R.valid), ...
        R.CycleDV_mps(R.valid), ...
        6,[0.75 0.75 0.75],'filled');

    scatter(TOF,DV,20, ...
        [0.10 0.45 0.85],'filled');

    xlabel('Complete cycle TOF [min]');
    ylabel('Complete cycle Delta-V [m/s]');

    title('FABER - Pareto Compromise Selection');

    legend('All feasible','Selected compromises', ...
        'Location','best');

end

%% =========================================================
% 8. REJECTION REASONS
% ==========================================================

failures = R.rejectionReason(:);
failures = failures(strlength(failures)>0);

if ~isempty(failures)

    [names,~,ii] = unique(failures);

    countsFail = accumarray(ii,1);

    Rejections = table(names,countsFail, ...
        'VariableNames',{'Reason','Count'});

    Rejections = sortrows(Rejections, ...
        'Count','descend');

    fprintf('\nTOP REJECTION REASONS\n');

    disp(Rejections( ...
        1:min(20,height(Rejections)),:));

else
    Rejections = table();
end

%% =========================================================
% 9. FINITE-THRUST VALIDATION
% ==========================================================

if ~strcmpi(cfg.finiteThrustValidation,'none')

    V = FABER_RPO_finite_validation( ...
        R,cfg.finiteThrustValidation);

    fprintf('\nFinite-thrust validation:\n');
    fprintf('Tested: %d\n',V.attemptedCount);
    fprintf('Passed: %d\n',V.passedCount);

    if V.attemptedCount > 0
        plot_FABER_RPO_finite(R,V);
        writetable(V.table,fullfile( ...
            cfg.outputFolder,'finite_validation.csv'));
    end
end

%% =========================================================
% 10. OPTIONAL SPHERICAL-SHELL STUDY
% ==========================================================

if runShellSensitivity

    shellCfg = cfg;

    shellCfg.initialPositionMode = 'shell';
    shellCfg.makePlots = false;
    shellCfg.saveCSV = false;

    R_shell = FABER_RPO_9phase_MC(shellCfg);

    fprintf('\nSHELL SENSITIVITY\n');

    fprintf('Nominal: %d/%d feasible\n', ...
        nnz(feasible),cfg.N);

    fprintf('Shell:   %d/%d feasible\n', ...
        nnz(any(R_shell.valid,2)),shellCfg.N);

end

%% =========================================================
% 11. EXPORT
% ==========================================================

SelectedResults = table( ...
    trialIdx,archIdx,DV,TOF,fuel, ...
    'VariableNames',{ ...
    'Trial','Architecture', ...
    'DV_mps','TOF_min','Fuel_kg'});

writetable(SelectedResults,fullfile( ...
    cfg.outputFolder,'selected_pareto_results.csv'));

writetable(ArchitectureResults,fullfile( ...
    cfg.outputFolder,'architecture_results.csv'));

writetable(PhaseResults,fullfile( ...
    cfg.outputFolder,'phase_results.csv'));

if ~isempty(Rejections)
    writetable(Rejections,fullfile( ...
        cfg.outputFolder,'rejection_reasons.csv'));
end

save(fullfile(cfg.outputFolder,'RPO_results.mat'), ...
    'R','cfg','SelectedResults','PhaseResults', ...
    'ArchitectureResults','-v7.3');

fprintf('\n=======================================\n');
fprintf('FABER RPO ANALYSIS COMPLETE\n');
fprintf('=======================================\n');
