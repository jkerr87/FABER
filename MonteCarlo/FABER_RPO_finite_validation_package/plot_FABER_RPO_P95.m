function plot_FABER_RPO_P95(R)
% PLOT_FABER_RPO_P95
%
% FABER nine-phase RPO Monte Carlo visualization.
%
% FIGURE 1:
%   P95 complete-cycle Delta-V comparison
%
% FIGURE 2:
%   P95 complete-cycle TOF comparison
%
% FIGURES 3 onward:
%   One 3D trajectory figure per architecture
%   - Left: Taxi -> Client
%   - Right: Taxi -> ATLANTE
%
% 3D trajectories correspond to actual accepted
% Monte Carlo trials at the nearest-rank P95
% complete-cycle Delta-V.
%
% ALL FIGURES ARE DOCKED IN MATLAB.

%% 1. SETTINGS

T = R.summary;
M = height(T);

folder = R.cfg.outputFolder;

if ~exist(folder,'dir')
    mkdir(folder);
end

% Dock all figures created by this function.
% Does not permanently change MATLAB's figure default.
oldStyle = get(groot,'DefaultFigureWindowStyle');
set(groot,'DefaultFigureWindowStyle','docked');

restoreStyle = onCleanup(@() ...
    set(groot,'DefaultFigureWindowStyle',oldStyle));

% Architecture labels
labels = strings(M,1);

for m = 1:M
    labels(m) = sprintf('%02d: %s / %s / %s', ...
        m, ...
        char(T.methodH2(m)), ...
        char(T.methodC1(m)), ...
        char(T.methodC2(m)));
end

% Determine architecture with lowest P95 DV.
eligible = find(T.Accepted > 0 & ...
    isfinite(T.P95CycleDV_mps));

if isempty(eligible)
    warning('No feasible architectures to plot.');
    return;
end

[~,ix] = min(T.P95CycleDV_mps(eligible));
bestArchitecture = eligible(ix);

fprintf('\n====================================\n');
fprintf('FABER RPO P95 VISUALIZATION\n');
fprintf('====================================\n');
fprintf('Architectures: %d\n',M);
fprintf('Lowest P95 DV architecture: %d\n', ...
    bestArchitecture);

%% ============================================================
% FIGURE 1: P95 DELTA-V
% =============================================================

fig1 = figure( ...
    'Name','01 - P95 Delta-V Comparison', ...
    'NumberTitle','off', ...
    'Color','w', ...
    'WindowStyle','docked');

ax1 = axes('Parent',fig1);

yDV = T.P95CycleDV_mps;

b = bar(ax1,1:M,yDV,0.75);

b.FaceColor = 'flat';
b.CData = repmat([0.30 0.53 0.78],M,1);

% Highlight lowest-P95-DV architecture
b.CData(bestArchitecture,:) = [0.90 0.38 0.20];

grid(ax1,'on');

set(ax1, ...
    'XTick',1:M, ...
    'XTickLabel',labels, ...
    'XTickLabelRotation',45, ...
    'FontSize',9);

xlabel(ax1,'Maneuver architecture');

ylabel(ax1,'P95 complete-cycle Delta-V [m/s]');

title(ax1,'FABER RPO - P95 Delta-V Comparison');

xlim(ax1,[0.3 M+0.7]);

finiteDV = yDV(isfinite(yDV));

if ~isempty(finiteDV)
    ylim(ax1,[0 1.15*max(finiteDV)]);
end

% Values above bars
for m = 1:M

    if isfinite(yDV(m))

        text(ax1,m,yDV(m), ...
            sprintf('%.2f',yDV(m)), ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom', ...
            'FontSize',8);
    end
end

exportgraphics(fig1, ...
    fullfile(folder,'P95_DeltaV_comparison.png'), ...
    'Resolution',200);

%% ============================================================
% FIGURE 2: P95 TIME OF FLIGHT
% =============================================================

fig2 = figure( ...
    'Name','02 - P95 TOF Comparison', ...
    'NumberTitle','off', ...
    'Color','w', ...
    'WindowStyle','docked');

ax2 = axes('Parent',fig2);

yTOF = T.P95CycleTOF_min;

b = bar(ax2,1:M,yTOF,0.75);

b.FaceColor = 'flat';
b.CData = repmat([0.25 0.62 0.52],M,1);

% Highlight lowest-P95-DV architecture
b.CData(bestArchitecture,:) = [0.90 0.38 0.20];

grid(ax2,'on');

set(ax2, ...
    'XTick',1:M, ...
    'XTickLabel',labels, ...
    'XTickLabelRotation',45, ...
    'FontSize',9);

xlabel(ax2,'Maneuver architecture');

ylabel(ax2,'P95 complete-cycle TOF [min]');

title(ax2,'FABER RPO - P95 Time of Flight Comparison');

xlim(ax2,[0.3 M+0.7]);

finiteTOF = yTOF(isfinite(yTOF));

if ~isempty(finiteTOF)
    ylim(ax2,[0 1.15*max(finiteTOF)]);
end

for m = 1:M

    if isfinite(yTOF(m))

        text(ax2,m,yTOF(m), ...
            sprintf('%.0f',yTOF(m)), ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom', ...
            'FontSize',8);
    end
end

exportgraphics(fig2, ...
    fullfile(folder,'P95_TOF_comparison.png'), ...
    'Resolution',200);

%% ============================================================
% FIGURES 3+: ONE 3D FIGURE PER ARCHITECTURE
% =============================================================

phaseNames = R.phaseNames;
colors = lines(9);

for m = 1:M

    %% Skip architectures with zero accepted trials

    acceptedTrials = find(R.valid(:,m));

    if isempty(acceptedTrials)

        fprintf('Architecture %d: no accepted trials.\n',m);
        continue;
    end

    %% Find actual P95 Delta-V trial

    dvAccepted = R.CycleDV_mps(acceptedTrials,m);

    [sortedDV,order] = sort(dvAccepted);

    rank95 = max(1,ceil(0.95*numel(sortedDV)));

    selectedTrial = acceptedTrials(order(rank95));

    rec = R.records{selectedTrial,m};

    if isempty(rec)
        warning('Empty record for architecture %d.',m);
        continue;
    end

    actualDV = R.CycleDV_mps(selectedTrial,m);
    actualTOF = R.CycleTOF_min(selectedTrial,m);

    %% Create one docked figure

    fig = figure( ...
        'Name',sprintf('%02d - %s / %s / %s', ...
            m+2, ...
            char(T.methodH2(m)), ...
            char(T.methodC1(m)), ...
            char(T.methodC2(m))), ...
        'NumberTitle','off', ...
        'Color','w', ...
        'WindowStyle','docked');

    tl = tiledlayout(fig,1,2, ...
        'TileSpacing','compact', ...
        'Padding','compact');

    %% Plot Client and ATLANTE side by side

    for leg = 1:2

        if leg==1

            result = rec.Client;
            draw = rec.designClient;
            destination = 'Taxi to Client';

        else

            result = rec.ATLANTE;
            draw = rec.designATLANTE;
            destination = 'Taxi to ATLANTE';
        end

        ax = nexttile(tl);

        hold(ax,'on');
        grid(ax,'on');
        axis(ax,'equal');

        %% Plot nine phases

        for k = 1:9

            X = result.trace{k};

            if isempty(X)
                continue;
            end

            % Skip a zero-length holding gate when drawing.
            if size(X,2)<2
                continue;
            end

            plot3(ax, ...
                X(1,:),X(2,:),X(3,:), ...
                'LineWidth',2, ...
                'Color',colors(k,:), ...
                'DisplayName',phaseNames{k});
        end

        %% AZ and KOZ

        [sx,sy,sz] = sphere(30);

        AZ = R.cfg.AZ_m;
        KOZ = R.cfg.KOZ_m;

        surf(ax, ...
            AZ*sx,AZ*sy,AZ*sz, ...
            'FaceColor',[0.20 0.50 0.90], ...
            'FaceAlpha',0.025, ...
            'EdgeColor',[0.20 0.50 0.90], ...
            'EdgeAlpha',0.10, ...
            'HandleVisibility','off');

        surf(ax, ...
            KOZ*sx,KOZ*sy,KOZ*sz, ...
            'FaceColor',[0.90 0.25 0.25], ...
            'FaceAlpha',0.08, ...
            'EdgeColor',[0.90 0.25 0.25], ...
            'EdgeAlpha',0.25, ...
            'HandleVisibility','off');

        %% Holding points

        for h = 1:4

            position = [0;-draw.hp(h);0];

            scatter3(ax, ...
                position(1), ...
                position(2), ...
                position(3), ...
                60,'k','filled', ...
                'HandleVisibility','off');

            text(ax, ...
                position(1), ...
                position(2), ...
                position(3), ...
                sprintf('  HP%d',h), ...
                'FontSize',9);
        end

        %% Target

        scatter3(ax,0,0,0, ...
            70,'r','filled', ...
            'HandleVisibility','off');

        text(ax,0,0,0,'  Target', ...
            'FontSize',9);

        %% Final standoff

        scatter3(ax,0,-draw.final,0, ...
            55,'m','filled', ...
            'HandleVisibility','off');

        %% Axes and view

        xlabel(ax,'Radial x [m]');
        ylabel(ax,'Along-track y [m]');
        zlabel(ax,'Cross-track z [m]');

        title(ax,destination);

        view(ax,35,25);

        if leg==1
            % Show one legend per architecture.
            legend(ax,'Location','bestoutside');
        end
    end

    %% Figure title

    title(tl,sprintf( ...
        ['Architecture %02d: %s / %s / %s\n' ...
         'P95-DV trial %d | Cycle DV %.3f m/s | ' ...
         'Cycle TOF %.1f min'], ...
        m, ...
        char(T.methodH2(m)), ...
        char(T.methodC1(m)), ...
        char(T.methodC2(m)), ...
        selectedTrial,actualDV,actualTOF));

    %% Export figure

    filename = sprintf( ...
        'P95_3D_Architecture_%02d.png',m);

    exportgraphics(fig, ...
        fullfile(folder,filename), ...
        'Resolution',200);

    fprintf(['Architecture %02d | P95-DV trial %d | ' ...
             'DV %.3f m/s | TOF %.1f min\n'], ...
        m,selectedTrial,actualDV,actualTOF);

end

%% ============================================================
% COMPLETE
% =============================================================

fprintf('\n====================================\n');
fprintf('P95 FIGURES COMPLETE\n');
fprintf('====================================\n');

fprintf('2 bar graphs + up to %d 3D figures.\n',M);
fprintf('All figures created as docked tabs.\n');
fprintf('Saved to: %s\n',folder);

end
