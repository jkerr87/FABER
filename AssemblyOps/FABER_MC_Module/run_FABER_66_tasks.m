%% FABER — Simple operation-by-operation module assembly Monte Carlo
% Goal: quantify assembly time for 50 valid modules, including lost elements.
% This is a schedule sensitivity model, NOT a spacecraft reliability model.
% REC-01...REC-15 add a single recoverable delay if triggered.
% REC-18 interrupts the current element attempt and requires a new element.
% REC-16/REC-17 are excluded: a critical spacecraft failure cannot be fixed
% merely by repeating the module attempt. Their risks are assessed separately.
clearvars; clc; close all;

%% User settings
nSim = 10000;                         % Independent campaigns per scenario
nTarget = 50;                        % Successfully installed modules
seed = 20261104;
recoveryCase = 'reference';          % 'low', 'reference', or 'high'
containerMode = 'once';             % 'once' or 'every'
restartAfterLoss_s = 0;             % Unknown reset allowance; zero omits repositioning
maxAttempts = 100;                  % Numerical safety guard, not an FDIR rule
scenarios = {'optimistic','reference','pessimistic'};

%% Read the original V4 tables (operation IDs must remain strings)
folder = fileparts(mfilename('fullpath'));
files = {'operations.csv','failure_assumptions.csv', ...
         'failure_events.csv','recovery_audit.csv'};
for k = 1:numel(files)
    assert(isfile(fullfile(folder,files{k})), ...
        'Missing required CSV: %s',files{k});
end
opts = detectImportOptions(fullfile(folder,files{1}),'TextType','string');
opts = setvartype(opts,'id','string');
O = readtable(fullfile(folder,files{1}),opts);
opts = detectImportOptions(fullfile(folder,files{2}),'TextType','string');
opts = setvartype(opts,'id','string');
F = readtable(fullfile(folder,files{2}),opts);
opts = detectImportOptions(fullfile(folder,files{3}),'TextType','string');
opts = setvartype(opts,{'event_id','recovery_id', ...
    'origin_operation','detection_operation'},'string');
E = readtable(fullfile(folder,files{3}),opts);
opts = detectImportOptions(fullfile(folder,files{4}),'TextType','string');
opts = setvartype(opts,'recovery_id','string');
A = readtable(fullfile(folder,files{4}),opts);
O = O(O.enabled == 1,:);
assert(height(O)==66 && height(E)==49, ...
    'Expected 66 nominal operations and 49 fault opportunities');
assert(all(ismember(E.detection_operation,O.id)), ...
    'Fault detection operation not found');
assert(all(ismember(E.recovery_id,F.id)), ...
    'Fault recovery ID not found');
assert(all(ismember(F.id,A.recovery_id)), ...
    'Recovery duration not found');
assert(ismember(string(recoveryCase),["low","reference","high"]));
assert(ismember(string(containerMode),["once","every"]));

% Map every event to its detection operation and recovery type.
[~, detectionIdx] = ismember(E.detection_operation,O.id);
[~, faultIdx] = ismember(E.recovery_id,F.id);
[~, auditIdx] = ismember(F.id,A.recovery_id);
activeAt = cell(height(O),1);
for e = 1:height(E)
    activeAt{detectionIdx(e)}(end+1) = e;
end
nOpportunity = accumarray(faultIdx,1,[height(F),1]);
recoverySeconds = A.(sprintf('%s_s',recoveryCase));
recoverySeconds = recoverySeconds(auditIdx);
assert(all(isfinite(recoverySeconds(1:15))));

% Build the nominal budget, keeping the original V4 margins.
nominalOperation = O.nominal_s + O.margin_s;
openIdx = find(O.id=="3.04");
assert(numel(openIdx)==1,'Container opening task is missing');
firstBaseline_s = sum(nominalOperation);
subsequentBaseline_s = firstBaseline_s - nominalOperation(openIdx);
assert(firstBaseline_s==5516 && subsequentBaseline_s==5406, ...
    'Unexpected nominal timing baseline');

%% Monte Carlo
rng(seed,'twister');
summary = cell(numel(scenarios),17);
allDelays = zeros(numel(scenarios),15);
for q = 1:numel(scenarios)
    scenario = scenarios{q};
    pCycle = F.(sprintf('p_%s',scenario));
    % Split a per-module REC probability across its mapped opportunities.
    % This assumes equal and independent opportunities; it is a sensitivity rule.
    pEvent = 1-(1-pCycle(faultIdx)).^(1./nOpportunity(faultIdx));

    campaignHours = zeros(nSim,1);
    totalAttempts = zeros(nSim,1);
    lostModules = zeros(nSim,1);
    delaySeconds = zeros(nSim,1);
    eventDelays = zeros(nSim,15);
    completed = false(nSim,1);
    moduleTimes_s = nan(nSim,nTarget); % Elapsed time to install each valid module, including failed attempts.

    for s = 1:nSim
        installed = 0;
        attempts = 0;
        elapsed = 0;
        lastInstalled_s = 0; % Start time of the next successful-module interval.
        while installed < nTarget && attempts < maxAttempts
            attempts = attempts + 1;
            moduleLost = false;
            for oi = 1:height(O)
                % Only the initial attempt pays for opening the container.
                if string(containerMode)=="once" && attempts>1 && oi==openIdx
                    continue
                end
                elapsed = elapsed + nominalOperation(oi);
                eventIndices = activeAt{oi};
                for e = eventIndices
                    r = faultIdx(e);
                    % REC-16 and REC-17 are critical VEHICLE faults, not module losses.
                    if r==16 || r==17
                        continue
                    end
                    if rand >= pEvent(e)
                        continue
                    end
                    if r==18
                        % Count already performed operations, then replace element.
                        moduleLost = true;
                        lostModules(s) = lostModules(s) + 1;
                        elapsed = elapsed + restartAfterLoss_s;
                        break
                    end
                    extra = recoverySeconds(r);
                    % REC-10 at 7.03 excludes future contact/engagement tasks.
                    if r==10 && O.id(oi)=="7.03"
                        if strcmp(recoveryCase,'reference')
                            extra = 35;
                        elseif strcmp(recoveryCase,'low')
                            extra = min(extra,35);
                        end
                    end
                    elapsed = elapsed + extra;
                    delaySeconds(s) = delaySeconds(s) + extra;
                    eventDelays(s,r) = eventDelays(s,r) + extra;
                end
                if moduleLost
                    break
                end
            end
            if ~moduleLost
                installed = installed + 1;
                % Attribute all work since the previous successful installation
                % to this valid module, including any element losses in between.
                moduleTimes_s(s,installed) = elapsed - lastInstalled_s;
                lastInstalled_s = elapsed;
            end
        end
        completed(s) = (installed==nTarget);
        campaignHours(s) = elapsed / 3600;
        totalAttempts(s) = attempts;
    end

    valid = completed;
    assert(any(valid),'No completed campaigns: increase maxAttempts or inspect assumptions');
    baselineHours = (nTarget*firstBaseline_s - ...
        strcmp(containerMode,'once')*(nTarget-1)*nominalOperation(openIdx))/3600;
    % Store the first ten campaign-level statistics
    summary(q,1:10) = {string(scenario),nSim,mean(valid), ...
        mean(totalAttempts),mean(lostModules), ...
        mean(campaignHours(valid)),prctile(campaignHours(valid),95), ...
        mean(delaySeconds(valid))/60/nTarget, ...
        mean(campaignHours(valid))*60/nTarget,baselineHours};
    allDelays(q,:) = mean(eventDelays(valid,:),1)/60/nTarget;
    % Compute module statistics from individual installed-module durations.
    completedModuleMin = moduleTimes_s(valid,:)/60;
    allModuleMin = completedModuleMin(:);
    assert(all(isfinite(allModuleMin)), ...
        'A completed campaign has missing module durations');
    moduleP = prctile(allModuleMin,[50 90 95 99]);
    firstMean = mean(completedModuleMin(:,1));
    laterMean = mean(completedModuleMin(:,2:end),'all');
    % Campaign time must equal the sum of individual module intervals.
    assert(max(abs(sum(completedModuleMin,2)*60 - campaignHours(valid)*3600)) < 1e-6, ...
        'Module durations do not sum to campaign duration');
    summary(q,11:17) = {mean(allModuleMin),moduleP(1),moduleP(2), ...
        moduleP(3),moduleP(4),firstMean,laterMean};
    % Save one row per successfully installed module.
    [simIdx,moduleIdx] = find(~isnan(moduleTimes_s));
    positions = sub2ind(size(moduleTimes_s),simIdx,moduleIdx);
    M = table(simIdx,moduleIdx,moduleTimes_s(positions)/60, ...
        'VariableNames',{'simulation','module_number','duration_min'});
    out = fullfile(folder,'results');
    if ~exist(out,'dir'), mkdir(out); end
    writetable(M,fullfile(out,sprintf('modules_%s.csv',scenario)));
    fprintf('%s: 50 modules %.2f h mean; P95 %.2f h; lost %.3f per campaign\n', ...
        upper(scenario),mean(campaignHours(valid)), ...
        prctile(campaignHours(valid),95),mean(lostModules));
    fprintf('  Each module: mean %.2f min | P50 %.2f | P90 %.2f | P95 %.2f | P99 %.2f min\n', ...
        mean(allModuleMin),moduleP);
    % Save per-campaign data for independent checks.
    T = table((1:nSim)',completed,campaignHours,totalAttempts,lostModules, ...
        delaySeconds/60,'VariableNames',{'simulation','completed', ...
        'campaign_h','attempts','lost_modules','recovery_delay_min'});
    out = fullfile(folder,'results');
    if ~exist(out,'dir'), mkdir(out); end
    writetable(T,fullfile(out,sprintf('campaigns_%s.csv',scenario)));
end

%% Results and plots
S = cell2table(summary,'VariableNames',{'scenario','nSim', ...
    'fraction_completed','mean_attempts','mean_lost_modules', ...
    'mean_campaign_h','campaign_P95_h','mean_recovery_min_per_module', ...
    'mean_time_min_per_module','baseline_campaign_h', ...
    'mean_individual_module_min','module_P50_min','module_P90_min', ...
    'module_P95_min','module_P99_min','first_module_mean_min', ...
    'subsequent_module_mean_min'});
out = fullfile(folder,'results');
writetable(S,fullfile(out,'summary.csv'));
D = array2table(allDelays,'VariableNames', ...
    cellstr(matlab.lang.makeValidName(F.id(1:15))));
D = addvars(D,string(scenarios(:)),'Before',1,'NewVariableNames','scenario');
writetable(D,fullfile(out,'recovery_contributions.csv'));
disp(S)
figure('Color','w'); bar(categorical(S.scenario), ...
    [S.mean_campaign_h S.campaign_P95_h]);
ylabel('Time to install 50 modules [h]');
legend('Mean','P95','Location','northwest');
title('FABER — Assembly time, including replacement attempts');
exportgraphics(gcf,fullfile(out,'campaign_times.png'),'Resolution',180);
figure('Color','w'); bar(categorical(S.scenario),S.mean_lost_modules);
ylabel('Mean modules lost per campaign');
title('FABER — Expected lost modules, not a validated risk rate');
exportgraphics(gcf,fullfile(out,'module_losses.png'),'Resolution',180);

figure('Color','w');
bar(categorical(S.scenario), ...
    [S.mean_individual_module_min S.module_P95_min]);
ylabel('Time per successfully installed module [min]');
legend('Mean','P95','Location','northwest');
title('FABER — Individual module assembly duration');
exportgraphics(gcf,fullfile(out,'module_times.png'),'Resolution',180);
