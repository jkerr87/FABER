
%% FABER GMAT ANNUAL ECLIPSE ANALYSIS
% GMAT EclipseLocator reports: ISS-like and SSO-like
% Epoch: 01 Jan 2038 to 01 Jan 2039
%
% IMPORTANT:
% Verify that both GMAT propagations actually cover the entire year.
% The EclipseLocator reports supplied previously contain identical
% events, so the ISS/SSO comparison is not yet validated.

clear; clc; close all;

%% USER SETTINGS

files = {
    'ASV_SAA_Results\FABER_ISS_Eclipse_2038.txt'
    'ASV_SAA_Results\FABER_SSO_Eclipse_2038.txt'
};

orbitNames = {'ISS-like','SSO-like'};
orbitTags = {'ISS','SSO'};

outDir = 'FABER_Eclipse_Results';

analysisStart = datetime(2038,1,1,0,0,0);
analysisEnd   = datetime(2039,1,1,0,0,0);

dt_min = 1;

% Set true ONLY after verifying full-year GMAT propagation.
fullYearVerified = [false false];

if ~exist(outDir,'dir')
    mkdir(outDir);
end

%% READ GMAT ECLIPSE REPORTS

reports = cell(2,1);

for k = 1:2
    reports{k} = readGMATEclipse(files{k});
end

% Detect identical orbit reports
if isequal(reports{1}.Start, reports{2}.Start) && ...
   isequal(reports{1}.End, reports{2}.End) && ...
   isequal(reports{1}.Type, reports{2}.Type)

    warning(['ISS-like and SSO-like eclipse reports have ' ...
             'identical event times. Verify GMAT exports.']);
end

%% ANNUAL TIME GRID

t = (analysisStart:minutes(dt_min): ...
     analysisEnd-minutes(dt_min))';

N = numel(t);
samplesPerDay = round(1440/dt_min);
nDays = N/samplesPerDay;

assert(mod(nDays,1)==0, ...
    'Simulation duration must contain whole days.');

dayDates = (analysisStart:days(1): ...
            analysisEnd-days(1))';

Summary = table();

%% PROCESS EACH ORBIT

for k = 1:2

    E = reports{k};

    %% GROUP INDIVIDUAL ECLIPSE EVENTS

    [eventID,~,group] = unique(E.EventNumber,'stable');

    nEvents = numel(eventID);

    eventStart = NaT(nEvents,1);
    totalDuration_min = zeros(nEvents,1);
    umbraDuration_min = zeros(nEvents,1);

    for j = 1:nEvents

        idx = group == j;

        eventStart(j) = min(E.Start(idx));

        totalDuration_min(j) = ...
            sum(E.Duration_s(idx))/60;

        umbraDuration_min(j) = ...
            sum(E.Duration_s(idx & E.Type=="Umbra"))/60;
    end

    %% CREATE ONE-YEAR ECLIPSE TIME HISTORY

    anyDiff = zeros(N+1,1);
    umbraDiff = zeros(N+1,1);

    for j = 1:height(E)

        % Ignore events completely outside analysis window
        if E.End(j) <= analysisStart || ...
           E.Start(j) >= analysisEnd
            continue
        end

        startTime = max(E.Start(j),analysisStart);
        endTime = min(E.End(j),analysisEnd);

        i1 = floor(seconds(startTime-analysisStart)/ ...
                   (60*dt_min))+1;

        i2 = ceil(seconds(endTime-analysisStart)/ ...
                  (60*dt_min))+1;

        i1 = max(1,min(N+1,i1));
        i2 = max(1,min(N+1,i2));

        anyDiff(i1) = anyDiff(i1)+1;
        anyDiff(i2) = anyDiff(i2)-1;

        if E.Type(j) == "Umbra"
            umbraDiff(i1) = umbraDiff(i1)+1;
            umbraDiff(i2) = umbraDiff(i2)-1;
        end
    end

    isEclipse = cumsum(anyDiff(1:N)) > 0;
    isUmbra = cumsum(umbraDiff(1:N)) > 0;

    isSunlight = ~isEclipse;

    %% DAILY STATISTICS

    eclipseMatrix = reshape(isEclipse, ...
        samplesPerDay,nDays);

    umbraMatrix = reshape(isUmbra, ...
        samplesPerDay,nDays);

    dailyEclipseFraction = ...
        mean(eclipseMatrix,1)'*100;

    dailyUmbraFraction = ...
        mean(umbraMatrix,1)'*100;

    dailySunlightFraction = ...
        100-dailyEclipseFraction;

    dailyEclipseMin = ...
        sum(eclipseMatrix,1)'*dt_min;

    %% SUMMARY STATISTICS

    maxTotal = max(totalDuration_min);
    maxUmbra = max(umbraDuration_min);

    meanTotal = mean(totalDuration_min);

    p95Total = prctile(totalDuration_min,95);

    fprintf('\n====================================\n');
    fprintf('FABER %s ECLIPSE ANALYSIS\n', ...
        orbitNames{k});
    fprintf('====================================\n');

    fprintf('First eclipse event: %s\n', ...
        string(min(E.Start)));

    fprintf('Last eclipse event: %s\n', ...
        string(max(E.End)));

    fprintf('Number of eclipse events: %d\n',nEvents);

    fprintf('Maximum total eclipse: %.2f min\n',maxTotal);
    fprintf('Maximum umbra: %.2f min\n',maxUmbra);

    fprintf('Mean total eclipse: %.2f min\n',meanTotal);
    fprintf('P95 total eclipse: %.2f min\n',p95Total);

    if fullYearVerified(k)
        fprintf('Annual sunlight fraction: %.2f %%\n', ...
            mean(isSunlight)*100);
    else
        warning(['Annual sunlight fraction for %s is ' ...
            'NOT VERIFIED. Confirm GMAT propagation coverage.'], ...
            orbitNames{k});
    end

    %% STORE SUMMARY

    row = table( ...
        string(orbitNames{k}), ...
        nEvents, ...
        maxTotal, ...
        maxUmbra, ...
        meanTotal, ...
        p95Total, ...
        mean(isSunlight)*100, ...
        fullYearVerified(k), ...
        'VariableNames',{ ...
        'Orbit', ...
        'EclipseEvents', ...
        'MaxEclipse_min', ...
        'MaxUmbra_min', ...
        'MeanEclipse_min', ...
        'P95Eclipse_min', ...
        'ComputedSunlightPercent', ...
        'FullYearVerified'});

    Summary = [Summary;row];

    %% SAVE TIME HISTORY

    Daily = table( ...
        dayDates, ...
        dailySunlightFraction, ...
        dailyEclipseFraction, ...
        dailyUmbraFraction, ...
        dailyEclipseMin, ...
        'VariableNames',{ ...
        'Date', ...
        'Sunlight_percent', ...
        'Eclipse_percent', ...
        'Umbra_percent', ...
        'Eclipse_min_per_day'});

    writetable(Daily,fullfile(outDir, ...
        [orbitTags{k} '_Daily_Eclipse.csv']));

    Events = table( ...
        eventStart,totalDuration_min,umbraDuration_min, ...
        'VariableNames',{ ...
        'StartUTC','TotalEclipse_min','Umbra_min'});

    writetable(Events,fullfile(outDir, ...
        [orbitTags{k} '_Eclipse_Events.csv']));

    %% FIGURE 1: ECLIPSE DURATION THROUGHOUT YEAR

    f = figure('Color','w', ...
        'Units','inches', ...
        'Position',[1 1 10 3.8]);

    plot(eventStart,totalDuration_min, ...
        'Color',[0.05 0.35 0.75], ...
        'LineWidth',1.2);

    hold on;

    plot(eventStart,umbraDuration_min, ...
        'Color',[0.85 0.25 0.15], ...
        'LineWidth',1.0);

    grid on;

    xlabel('Date in 2038','Color','k');
    ylabel('Duration [min]','Color','k');

    title([orbitNames{k} ...
        ' - Eclipse Duration Throughout 2038'], ...
        'Color','k');

    legend('Total eclipse','Umbra', ...
        'Location','best', ...
        'Color','w','TextColor','k');

    xlim([analysisStart analysisEnd]);

    set(gca,'Color','w', ...
        'XColor','k','YColor','k', ...
        'FontSize',11);

    exportgraphics(f,fullfile(outDir, ...
        [orbitTags{k} '_Eclipse_Duration.png']), ...
        'Resolution',300);

    %% FIGURE 2: DAILY SUNLIGHT AVAILABILITY

    f = figure('Color','w', ...
        'Units','inches', ...
        'Position',[1 1 10 3.8]);

    plot(dayDates,dailySunlightFraction, ...
        'Color',[0.90 0.55 0.05], ...
        'LineWidth',1.5);

    hold on;

    plot(dayDates,dailyEclipseFraction, ...
        'Color',[0.15 0.35 0.75], ...
        'LineWidth',1.2);

    grid on;

    xlabel('Date in 2038','Color','k');
    ylabel('Fraction of day [%]','Color','k');

    title([orbitNames{k} ...
        ' - Daily Sunlight and Eclipse Fraction'], ...
        'Color','k');

    legend('Sunlight','Eclipse', ...
        'Location','best', ...
        'Color','w','TextColor','k');

    ylim([0 105]);
    xlim([analysisStart analysisEnd]);

    set(gca,'Color','w', ...
        'XColor','k','YColor','k', ...
        'FontSize',11);

    exportgraphics(f,fullfile(outDir, ...
        [orbitTags{k} '_Daily_Sunlight.png']), ...
        'Resolution',300);

    %% FIGURE 3: ECLIPSE BEHAVIOUR OVER ONE YEAR

    % 0 = full sunlight
    % 1 = penumbra
    % 2 = umbra

    state = zeros(N,1);

    state(isEclipse) = 1;
    state(isUmbra) = 2;

    % Plot at 10-minute resolution for manageable graphics.
    plotStep = max(1,round(10/dt_min));

    f = figure('Color','w', ...
        'Units','inches', ...
        'Position',[1 1 11 3.2]);

    stairs(t(1:plotStep:end), ...
        state(1:plotStep:end), ...
        'Color',[0.1 0.25 0.6], ...
        'LineWidth',0.65);

    grid on;

    xlabel('Date in 2038','Color','k');
    ylabel('Solar condition','Color','k');

    title([orbitNames{k} ...
        ' - Annual Eclipse Timeline'], ...
        'Color','k');

    yticks([0 1 2]);
    yticklabels({'Sunlight','Penumbra','Umbra'});

    ylim([-0.2 2.2]);
    xlim([analysisStart analysisEnd]);

    set(gca,'Color','w', ...
        'XColor','k','YColor','k', ...
        'FontSize',10);

    exportgraphics(f,fullfile(outDir, ...
        [orbitTags{k} '_Eclipse_Timeline.png']), ...
        'Resolution',300);

    %% FIGURE 4: ECLIPSE DURATION DISTRIBUTION

    f = figure('Color','w', ...
        'Units','inches', ...
        'Position',[1 1 8 3.8]);

    histogram(totalDuration_min,30, ...
        'FaceColor',[0.15 0.45 0.75], ...
        'EdgeColor','w');

    hold on;

    xline(p95Total,'--r', ...
        sprintf('P95 = %.1f min',p95Total), ...
        'LineWidth',1.4);

    xline(maxTotal,'--k', ...
        sprintf('Max = %.1f min',maxTotal), ...
        'LineWidth',1.4);

    grid on;

    xlabel('Eclipse duration [min]','Color','k');
    ylabel('Number of events','Color','k');

    title([orbitNames{k} ...
        ' - Eclipse Duration Distribution'], ...
        'Color','k');

    set(gca,'Color','w', ...
        'XColor','k','YColor','k', ...
        'FontSize',11);

    exportgraphics(f,fullfile(outDir, ...
        [orbitTags{k} '_Eclipse_Histogram.png']), ...
        'Resolution',300);

end

%% SAVE COMPARISON

disp(Summary);

writetable(Summary, ...
    fullfile(outDir,'FABER_Eclipse_Comparison.csv'));

fprintf('\nAll plots saved to:\n%s\n', ...
    fullfile(pwd,outDir));

%% LOCAL FUNCTION: READ GMAT ECLIPSELOCATOR REPORT

function E = readGMATEclipse(filename)

    if ~isfile(filename)
        error('Cannot find GMAT report: %s',filename);
    end

    lines = readlines(filename);

    % GMAT report fields:
    % Start Time | Stop Time | Duration (s) |
    % Occulting Body | Type | Event Number | Total Duration

    expression = [ ...
        '^\s*(\d{2}\s+\w{3}\s+\d{4}\s+' ...
        '\d{2}:\d{2}:\d{2}\.\d+)\s+' ...
        '(\d{2}\s+\w{3}\s+\d{4}\s+' ...
        '\d{2}:\d{2}:\d{2}\.\d+)\s+' ...
        '([\d.Ee+-]+)\s+' ...
        '(\w+)\s+' ...
        '(Penumbra|Umbra|Antumbra)\s+' ...
        '(\d+)\s+' ...
        '([\d.Ee+-]+)'];

    startList = datetime.empty(0,1);
    endList = datetime.empty(0,1);

    durationList = [];
    typeList = strings(0,1);
    numberList = [];

    for i = 1:numel(lines)

        tokens = regexp(char(lines(i)), ...
            expression,'tokens','once');

        if isempty(tokens)
            continue
        end

        startList(end+1,1) = datetime(tokens{1}, ...
            'InputFormat','dd MMM yyyy HH:mm:ss.SSS', ...
            'Locale','en_US');

        endList(end+1,1) = datetime(tokens{2}, ...
            'InputFormat','dd MMM yyyy HH:mm:ss.SSS', ...
            'Locale','en_US');

        durationList(end+1,1) = str2double(tokens{3});

        typeList(end+1,1) = string(tokens{5});

        numberList(end+1,1) = str2double(tokens{6});
    end

    if isempty(durationList)
        error('No eclipse events parsed from %s',filename);
    end

    E = table( ...
        startList,endList,durationList, ...
        typeList,numberList, ...
        'VariableNames',{ ...
        'Start','End','Duration_s', ...
        'Type','EventNumber'});

end
