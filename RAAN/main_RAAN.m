
%% FABER RAAN ACCESSIBILITY - BOTH CLUSTERS
% 30, 90, 180-day screening with 100 m/s DV cap
%
% Requires:
%   FABER_RAAN_Phase_Joint30d.m
%
% Outputs:
%   1. RAAN feasibility envelope plot
%   2. Delta-V versus RAAN separation plot
%   3. Summary CSV
%   4. Complete MATLAB results
%
% -------------------------------------------------------

clearvars; close all; clc;

%% CONFIGURATION

Tdays = [30 90 180];
DVcap = 100;              % [m/s]
phase0 = 20;              % [deg]
raanStep = 1;             % [deg]

% FABER cluster reference orbits
clusterNames = { ...
    'Cluster 1: 457.25 km / 51.02 deg', ...
    'Cluster 2: 680.12 km / 98.12 deg'};

hReference_km = [457.25, 680.12];
iReference_deg = [51.02, 98.12];

% Drift-altitude search range
hDriftCandidates_km = 380:1:700;

% Output directory
outputFolder = 'FABER_RAAN_report_outputs';

if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end

%% RUN ALL SIX CASES

Results = cell(2,numel(Tdays));
summaryRows = cell(0,10);

for c = 1:2

    for j = 1:numel(Tdays)

        cfg = struct();

        cfg.hReference_km = hReference_km(c);
        cfg.iReference_deg = iReference_deg(c);

        cfg.DVcap_mps = DVcap;
        cfg.Tmax_days = Tdays(j);

        cfg.phaseDiff_deg = phase0;
        cfg.raanDiff_deg = -180:raanStep:180;

        cfg.hDriftCandidates_km = hDriftCandidates_km;

        cfg.maxPoweredRAAN_deg = 0.5;

        cfg.makeFigures = false;
        cfg.exportFiles = false;

        fprintf('\n--- %s | Tmax = %d days ---\n', ...
            clusterNames{c},Tdays(j));

        R = FABER_RAAN_Phase_Joint30d(cfg);
        Results{c,j} = R;

        raan = R.cfg.raanDiff_deg(:)';
        feasible = logical(R.feasible(1,:));

        % Overall feasible RAAN range
        if any(feasible)

            allMin = min(raan(feasible));
            allMax = max(raan(feasible));

            maxAbs = max(abs(raan(feasible)));

        else

            allMin = NaN;
            allMax = NaN;
            maxAbs = NaN;

        end

        % Continuous feasible interval containing zero
        cMin = NaN;
        cMax = NaN;

        i0 = find(abs(raan)<1e-10,1);

        if ~isempty(i0) && feasible(i0)

            iLeft = i0;
            iRight = i0;

            while iLeft > 1 && feasible(iLeft-1)
                iLeft = iLeft - 1;
            end

            while iRight < numel(raan) && feasible(iRight+1)
                iRight = iRight + 1;
            end

            cMin = raan(iLeft);
            cMax = raan(iRight);

        end

        coverage = 100*nnz(feasible)/numel(feasible);

        summaryRows(end+1,:) = { ...
            clusterNames{c}, ...
            hReference_km(c), ...
            iReference_deg(c), ...
            Tdays(j), ...
            allMin, ...
            allMax, ...
            cMin, ...
            cMax, ...
            maxAbs, ...
            coverage};

        fprintf('Feasible RAAN range: %.1f to %.1f deg\n', ...
            allMin,allMax);

        fprintf('Continuous envelope: %.1f to %.1f deg\n', ...
            cMin,cMax);

    end
end

%% SUMMARY TABLE

Summary = cell2table(summaryRows, ...
    'VariableNames',{ ...
    'Cluster', ...
    'Altitude_km', ...
    'Inclination_deg', ...
    'Tmax_days', ...
    'MinFeasible_deg', ...
    'MaxFeasible_deg', ...
    'ContinuousMin_deg', ...
    'ContinuousMax_deg', ...
    'MaxAbsRAAN_deg', ...
    'FeasibleGrid_percent'});

fprintf('\nFABER RAAN ACCESSIBILITY SUMMARY\n');
disp(Summary);

% writetable(Summary,fullfile(outputFolder, ...
%     'FABER_RAAN_two_clusters_summary.csv'));
% 
% save(fullfile(outputFolder, ...
%     'FABER_RAAN_two_clusters_results.mat'), ...
%     'Results','Summary','Tdays','DVcap','phase0', ...
%     'hReference_km','iReference_deg');

%% AVERAGE SELECTED TOF FOR EACH RAAN ENVELOPE

AverageTOF_days = nan(2,numel(Tdays));

for c = 1:2
    for j = 1:numel(Tdays)

        R = Results{c,j};
        raan = R.cfg.raanDiff_deg(:)';

        % Continuous envelope from summary
        row = (c-1)*numel(Tdays) + j;
        raanMin = Summary.ContinuousMin_deg(row);
        raanMax = Summary.ContinuousMax_deg(row);

        % Include only feasible RAAN points in the envelope
        mask = R.feasible(1,:) & ...
            raan >= raanMin & ...
            raan <= raanMax;

        if any(mask)
            AverageTOF_days(c,j) = ...
                mean(R.bestTOF_days(1,mask),'omitnan');
        end

        fprintf(['Cluster %d | Tmax = %3d days | ', ...
            'Average selected TOF = %.2f days\n'], ...
            c,Tdays(j),AverageTOF_days(c,j));
    end
end

% Add average selected TOF to your summary table
Summary.AverageSelectedTOF_days = ...
    reshape(AverageTOF_days.',[],1);

disp(Summary);
%% FIGURE 1 - RAAN FEASIBILITY ENVELOPES

figure('Color','w', ...
    'Name','FABER RAAN Servicing Envelopes', ...
    'Position',[80 80 1120 580]);

tiledlayout(2,1, ...
    'TileSpacing','compact', ...
    'Padding','compact');

colors = [ ...
    0.22 0.48 0.70;
    0.92 0.55 0.19;
    0.31 0.62 0.45];

for c = 1:2

    ax = nexttile;
    hold(ax,'on');

    for j = 1:numel(Tdays)

        R = Results{c,j};

        x = R.cfg.raanDiff_deg(:)';
        mask = logical(R.feasible(1,:));

        % Plot only feasible points.
        % This preserves gaps in the RAAN envelope.
        plot(ax,x(mask),j*ones(1,nnz(mask)), ...
            's', ...
            'MarkerFaceColor',colors(j,:), ...
            'MarkerEdgeColor',colors(j,:), ...
            'MarkerSize',5, ...
            'LineStyle','none');

    end

    xline(ax,0,':', ...
        'Color',[0.35 0.35 0.35]);

    yticks(ax,1:3);
    yticklabels(ax,{'30 days','90 days','180 days'});

    ylim(ax,[0.5 3.5]);

    grid(ax,'on');
    box(ax,'off');

    ylabel(ax,'Maximum transfer duration');
    title(ax,clusterNames{c},'FontWeight','bold');

    if c == 2
        xlabel(ax,'Initial relative RAAN, \Delta\Omega [deg]');
    end

end

% Set shared x-axis range from feasible results
allRanges = [];

for c = 1:2
    for j = 1:3

        R = Results{c,j};

        x = R.cfg.raanDiff_deg(:)';
        f = R.feasible(1,:);

        allRanges = [allRanges x(f)]; %#ok<AGROW>

    end
end

if ~isempty(allRanges)

    lim = [floor(min(allRanges))-2, ...
           ceil(max(allRanges))+2];

    axs = findall(gcf,'Type','axes');

    for k = 1:numel(axs)
        xlim(axs(k),lim);
    end

end

sgtitle(sprintf( ...
    'FABER RAAN Accessibility | DV <= %.0f m/s | Initial phase = %.0f deg', ...
    DVcap,phase0));

% exportgraphics(gcf, ...
%     fullfile(outputFolder,'FABER_RAAN_envelopes.png'), ...
%     'Resolution',300);
% 
% exportgraphics(gcf, ...
%     fullfile(outputFolder,'FABER_RAAN_envelopes.pdf'), ...
%     'ContentType','vector');
% 
% savefig(gcf, ...
%     fullfile(outputFolder,'FABER_RAAN_envelopes.fig'));

%% FIGURE 2 - DELTA-V VS RAAN SEPARATION

figure('Color','w', ...
    'Name','FABER Delta-V vs RAAN', ...
    'Position',[100 100 1120 690]);

tiledlayout(2,1, ...
    'TileSpacing','compact', ...
    'Padding','compact');

for c = 1:2

    ax = nexttile;
    hold(ax,'on');

    for j = 1:numel(Tdays)

        R = Results{c,j};

        x = R.cfg.raanDiff_deg(:)';
        y = R.minDV_mps(1,:);

        % NaN values break lines at infeasible points
        plot(ax,x,y,'-o', ...
            'Color',colors(j,:), ...
            'MarkerSize',3, ...
            'LineWidth',1.5, ...
            'DisplayName',sprintf('%d days',Tdays(j)));

    end

    yline(ax,DVcap,'k--','100 m/s cap');

    xline(ax,0,':', ...
        'Color',[0.4 0.4 0.4]);

    ylabel(ax,'Minimum feasible \Delta v [m/s]');

    title(ax,clusterNames{c},'FontWeight','bold');

    grid(ax,'on');
    box(ax,'off');

    ylim(ax,[0 110]);

    if ~isempty(allRanges)
        xlim(ax,lim);
    end

    if c == 1
        legend(ax,'Location','best');
    end

    if c == 2
        xlabel(ax,'Initial relative RAAN, \Delta\Omega [deg]');
    end

end

sgtitle(sprintf( ...
    'Joint J_2 Drift and Phasing Screening | Initial phase = %.0f deg', ...
    phase0));

% exportgraphics(gcf, ...
%     fullfile(outputFolder,'FABER_RAAN_DV_vs_separation.png'), ...
%     'Resolution',300);
% 
% exportgraphics(gcf, ...
%     fullfile(outputFolder,'FABER_RAAN_DV_vs_separation.pdf'), ...
%     'ContentType','vector');

% savefig(gcf, ...
%     fullfile(outputFolder,'FABER_RAAN_DV_vs_separation.fig'));

%% COMPLETION MESSAGE

fprintf('\n========================================\n');
fprintf('FABER RAAN SCREENING COMPLETE\n');
fprintf('========================================\n');
fprintf('Results exported to: %s\n',outputFolder);
