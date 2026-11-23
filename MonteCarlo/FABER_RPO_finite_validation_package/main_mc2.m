
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
%cfg.h_km = 680.12;
cfg.h_km = 457.25;

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
cfg.maxBurnFraction = 0.10;
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
cfg.finiteThrustValidation = 'all';
cfg.finiteMaxCandidates = 72;

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

cfg.makePlots = true;
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
% INITIAL TAXI DEPARTURE FROM ATLANTE
% Separate from the nominal RPO retrieval-cycle budget
% ==========================================================

D = struct();

if cfg.departure.enabled

    D = FABER_Taxi_initial_departure(cfg);

    fprintf('\n====================================\n');
    fprintf('INITIAL TAXI DEPARTURE FROM ATLANTE\n');
    fprintf('====================================\n');

    fprintf('Feasible:         %d\n',D.feasible);
    fprintf('Delta-V:          %.5f m/s\n',D.DV_mps);
    fprintf('TOF:              %.2f min\n',D.TOF_min);
    fprintf('Propellant:       %.5f kg\n',D.propellant_kg);
    fprintf('Active thrust:    %.2f s\n',D.burnTime_s);
    fprintf('Peak thrust:      %.3f N\n',D.maxThrust_N);
    fprintf('Available thrust: %.3f N\n', ...
        D.availableThrust_N);

    if ~D.feasible
        warning('Initial Taxi departure failed one or more checks.');
    end

    figure('Color','w','WindowStyle','docked');

    tiledlayout(3,1);

    nexttile;
    plot(D.time_s/60, ...
        vecnorm(D.thrust_LVLH_N,2,2));
    ylabel('Thrust [N]');
    grid on;

    nexttile;
    plot(D.time_s/60,D.rho_m);
    ylabel('Range [m]');
    grid on;

    nexttile;
    plot(D.time_s/60,D.rhoDot_mps);
    xlabel('Time [min]');
    ylabel('Range rate [m/s]');
    grid on;

end


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


%% FINITE-THRUST PLOTS FOR ALL TESTED ARCHITECTURES

% if V.attemptedCount > 0
% 
%     M = size(R.valid,2);
% 
%     fprintf('\nGENERATING FINITE-THRUST PLOTS\n');
% 
%     for m = 1:M
% 
%         % Find a passing finite-thrust candidate
%         i = find(V.passed(:,m),1,'first');
% 
%         if isempty(i)
%             fprintf('Architecture %d: no passing candidate.\n',m);
%             continue;
%         end
% 
%         fprintf('Plotting trial %d, architecture %d\n',i,m);
% 
%         % Generate plot in MATLAB only
%         plot_FABER_RPO_finite(R,V,i,m);
% 
%     end
% 
% end
% 
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


%% PLOT PREFERRED PARETO ARCHITECTURE

T = R.summary;

m = find( ...
    T.methodH2 == "Vbar" & ...
    T.methodC1 == "Rbar" & ...
    T.methodC2 == "GeneralHCW" & ...
    T.gateMode == "FlyThrough");

assert(isscalar(m),'Architecture not found uniquely.');

% Find a passing finite-thrust trial
i = find(V.passed(:,m),1,'first');

if ~isempty(i)

    plot_FABER_RPO_finite(R,V,i,m);

    fprintf('Preferred architecture: %d, trial %d\n',m,i);

else

    warning(['No passing finite-thrust validation ' ...
        'available for the preferred architecture.']);

end



%% ================================================================
% 10. RPO TRAJECTORY GEOMETRY - REPORT FIGURES
% Preferred architecture: Vbar / Rbar / GeneralHCW / FlyThrough
% Two separate horizontal figures
% ================================================================

% Find preferred architecture
Tarch = R.summary;

mPlot = find( ...
    string(Tarch.methodH2) == "Vbar" & ...
    string(Tarch.methodC1) == "Rbar" & ...
    string(Tarch.methodC2) == "GeneralHCW" & ...
    string(Tarch.gateMode) == "FlyThrough");

assert(isscalar(mPlot), ...
    'Preferred RPO architecture not found uniquely.');

% Select accepted trial closest to median cycle Delta-V
validTrials = find(R.valid(:,mPlot));

assert(~isempty(validTrials), ...
    'No accepted trials for preferred architecture.');

dvValid = R.CycleDV_mps(validTrials,mPlot);
dvMedian = median(dvValid);

[~,idx] = min(abs(dvValid - dvMedian));
iPlot = validTrials(idx);

% Extract actual Client approach trajectories
rec = R.records{iPlot,mPlot};
traces = rec.Client.trace;

% Nine-phase sequence:
% 1 Homing 1
% 2 Holding 1
% 3 Homing 2
% 4 Holding 2
% 5 Closing 1
% 6 Holding 3
% 7 Closing 2
% 8 HP4 gate
% 9 Final approach

phaseIdx = [1 3 5 7 9];

phaseLabels = { ...
    'Homing 1 (HCW)', ...
    'Homing 2 (V-bar)', ...
    'Closing 1 (R-bar)', ...
    'Closing 2 (HCW)', ...
    'Final approach'};

colors = [ ...
    0.13 0.42 0.70;
    0.12 0.60 0.48;
    0.85 0.45 0.15;
    0.58 0.35 0.72;
    0.80 0.20 0.30];

% Safety-zone geometry
AZ = cfg.AZ_m;
KOZ = cfg.KOZ_m;

theta = linspace(0,2*pi,500);

xAZ = AZ*cos(theta);
yAZ = AZ*sin(theta);

xKOZ = KOZ*cos(theta);
yKOZ = KOZ*sin(theta);

% Figure settings
lw = 2.2;
fs = 12;
ms = 7;

%% ================================================================
% FIGURE 1 - COMPLETE RPO TRAJECTORY
% ================================================================

fig1 = figure( ...
    'Color','w', ...
    'Position',[100 100 1100 550], ...
    'Name','FABER - Complete RPO Trajectory');

ax1 = axes(fig1);

hold(ax1,'on');
grid(ax1,'on');
box(ax1,'on');
axis(ax1,'equal');

set(ax1, ...
    'Color','w', ...
    'XColor','k', ...
    'YColor','k', ...
    'GridColor',[0.75 0.75 0.75], ...
    'GridAlpha',0.3, ...
    'FontName','Times New Roman', ...
    'FontSize',fs);

% Safety zones
plot(ax1,xAZ,yAZ,'--', ...
    'Color',[0.4 0.4 0.4], ...
    'LineWidth',1.3, ...
    'DisplayName','Approach Zone');

plot(ax1,xKOZ,yKOZ,':', ...
    'Color',[0.8 0.2 0.2], ...
    'LineWidth',1.5, ...
    'DisplayName','Keep-Out Zone');

% All five maneuver phases
for k = 1:numel(phaseIdx)

    X = traces{phaseIdx(k)};

    if isempty(X)
        continue;
    end

    plot(ax1,X(1,:),X(2,:), ...
        'Color',colors(k,:), ...
        'LineWidth',lw, ...
        'DisplayName',phaseLabels{k});

    % Phase endpoint
    plot(ax1,X(1,end),X(2,end),'o', ...
        'MarkerSize',ms, ...
        'MarkerFaceColor',colors(k,:), ...
        'MarkerEdgeColor','k', ...
        'HandleVisibility','off');

end

% Initial position
X0 = traces{1};

plot(ax1,X0(1,1),X0(2,1),'ks', ...
    'MarkerFaceColor','k', ...
    'MarkerSize',8, ...
    'DisplayName','Initial position');

% Target
plot(ax1,0,0,'p', ...
    'MarkerEdgeColor','k', ...
    'MarkerFaceColor',[1 0.8 0], ...
    'MarkerSize',14, ...
    'DisplayName','Target');

% Holding point labels
for k = 1:4

    X = traces{phaseIdx(k)};

    if isempty(X)
        continue;
    end

    px = X(1,end);
    py = X(2,end);

    text(ax1,px+20,py, ...
        sprintf('HP%d',k), ...
        'FontSize',fs-1, ...
        'Color','k');

end

xlabel(ax1,'Radial position x [m]');
ylabel(ax1,'Along-track position y [m]');

%title(ax1,'Complete RPO trajectory', ...
%    'FontWeight','bold');

title(ax1,'Complete RPO trajectory', ...
    'Color','k', ...
    'FontSize',14, ...
    'FontWeight','bold');

%legend(ax1,'Location','eastoutside');
lgd1 = legend(ax1,'Location','eastoutside');
lgd1.Color = 'w';
lgd1.TextColor = 'k';
lgd1.EdgeColor = [0.7 0.7 0.7];
lgd1.FontSize = 10;

% Save report figure
exportgraphics(fig1, ...
    fullfile(cfg.outputFolder,'RPO_full_trajectory.png'), ...
    'Resolution',300);

exportgraphics(fig1, ...
    fullfile(cfg.outputFolder,'RPO_full_trajectory.pdf'), ...
    'ContentType','vector');


%% ================================================================
% FIGURE 2 - CLOSE-RANGE APPROACH
% ================================================================

fig2 = figure( ...
    'Color','w', ...
    'Position',[150 150 1100 550], ...
    'Name','FABER - Close-Range RPO');

ax2 = axes(fig2);

hold(ax2,'on');
grid(ax2,'on');
box(ax2,'on');
axis(ax2,'equal');

set(ax2, ...
    'Color','w', ...
    'XColor','k', ...
    'YColor','k', ...
    'GridColor',[0.75 0.75 0.75], ...
    'GridAlpha',0.3, ...
    'FontName','Times New Roman', ...
    'FontSize',fs);

% Safety zones
plot(ax2,xAZ,yAZ,'--', ...
    'Color',[0.4 0.4 0.4], ...
    'LineWidth',1.3, ...
    'DisplayName','Approach Zone');

plot(ax2,xKOZ,yKOZ,':', ...
    'Color',[0.8 0.2 0.2], ...
    'LineWidth',1.5, ...
    'DisplayName','Keep-Out Zone');

% Closing 1, Closing 2, Final approach
for k = 3:5

    X = traces{phaseIdx(k)};

    if isempty(X)
        continue;
    end

    plot(ax2,X(1,:),X(2,:), ...
        'Color',colors(k,:), ...
        'LineWidth',lw, ...
        'DisplayName',phaseLabels{k});

    plot(ax2,X(1,end),X(2,end),'o', ...
        'MarkerSize',ms, ...
        'MarkerFaceColor',colors(k,:), ...
        'MarkerEdgeColor','k', ...
        'HandleVisibility','off');

end

% Target
plot(ax2,0,0,'p', ...
    'MarkerEdgeColor','k', ...
    'MarkerFaceColor',[1 0.8 0], ...
    'MarkerSize',14, ...
    'DisplayName','Target');

% Holding points 3 and 4
for k = 3:4

    X = traces{phaseIdx(k)};

    if isempty(X)
        continue;
    end

    px = X(1,end);
    py = X(2,end);

    text(ax2,px+10,py, ...
        sprintf('HP%d',k), ...
        'FontSize',fs-1, ...
        'Color','k');

end

xlabel(ax2,'Radial position x [m]');
ylabel(ax2,'Along-track position y [m]');

%title(ax2,'Close-range RPO trajectory', ...
%    'FontWeight','bold');

title(ax2,'Close-range RPO trajectory', ...
    'Color','k', ...
    'FontSize',14, ...
    'FontWeight','bold');

% Zoom into close-range region
xlim(ax2,[-350 350]);
ylim(ax2,[-360 60]);

%legend(ax2,'Location','eastoutside');
lgd2 = legend(ax2,'Location','eastoutside');
lgd2.Color = 'w';
lgd2.TextColor = 'k';
lgd2.EdgeColor = [0.7 0.7 0.7];
lgd2.FontSize = 10;

% Save report figure
exportgraphics(fig2, ...
    fullfile(cfg.outputFolder,'RPO_close_range.png'), ...
    'Resolution',300);

exportgraphics(fig2, ...
    fullfile(cfg.outputFolder,'RPO_close_range.pdf'), ...
    'ContentType','vector');


%% ================================================================
% PLOT INFORMATION
% ================================================================

fprintf('\n=======================================\n');
fprintf('RPO TRAJECTORY FIGURES GENERATED\n');
fprintf('=======================================\n');

fprintf('Architecture: %d\n',mPlot);
fprintf('Trial: %d\n',iPlot);
fprintf('Cycle Delta-V: %.4f m/s\n', ...
    R.CycleDV_mps(iPlot,mPlot));

fprintf('Figures saved in:\n%s\n',cfg.outputFolder);


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

