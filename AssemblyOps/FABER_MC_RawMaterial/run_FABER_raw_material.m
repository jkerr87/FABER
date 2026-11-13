%% FABER — Scenario B: semi-finished raw-material element assembly Monte Carlo
% Goal: quantify assembly time for 50 valid modules, including lost elements.
% This is a schedule sensitivity model, NOT a spacecraft reliability model.
% Raw-material fabrication is represented by one trim and one joining step,
% NOT processing bulk feedstock into a finished part.
% RAW-06 replacement assumes safe removal of the defective joint; this is
% unvalidated and may actually require terminating the campaign.
% Applicable REC events and RAW-01...RAW-05 add a single recoverable delay.
% REC-18 and RAW-06 stop an element attempt and start a replacement.
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

%% Read the Scenario B tables (operation IDs must remain strings)
folder = fileparts(mfilename('fullpath'));
files = {'operations_raw.csv','failure_assumptions_raw.csv', ...
         'failure_events_raw.csv','recovery_audit_raw.csv'};
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
assert(height(O)==70 && height(E)==53, ...
    'Expected 70 nominal operations and 53 fault opportunities');
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
recoverable = ~ismember(F.id,["REC-16","REC-17","REC-18","RAW-06"]);
assert(all(isfinite(recoverySeconds(recoverable))));

% Build the nominal budget, keeping the original V4 margins.
nominalOperation = O.nominal_s + O.margin_s;
openIdx = find(O.id=="3.04");
assert(numel(openIdx)==1,'Container opening task is missing');
firstBaseline_s = sum(nominalOperation);
subsequentBaseline_s = firstBaseline_s - nominalOperation(openIdx);
assert(firstBaseline_s==6156 && subsequentBaseline_s==6046, ...
    'Unexpected nominal timing baseline');

%% Monte Carlo
rng(seed,'twister');
summary = cell(numel(scenarios),17);
recoverableIDs = F.id(recoverable);
allDelays = zeros(numel(scenarios),numel(recoverableIDs));
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
    eventDelays = zeros(nSim,numel(recoverableIDs));
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
                    % Critical vehicle failures cannot be represented by a bounded retry.
                    if ismember(F.id(r),["REC-16","REC-17"])
                        continue
                    end
                    if rand >= pEvent(e)
                        continue
                    end
                    if ismember(F.id(r),["REC-18","RAW-06"])
                        % Count already performed operations, then replace element.
                        moduleLost = true;
                        lostModules(s) = lostModules(s) + 1;
                        elapsed = elapsed + restartAfterLoss_s;
                        break
                    end
                    extra = recoverySeconds(r);
                    % REC-10 at 7.03 excludes future contact/engagement tasks.
                    if F.id(r)=="REC-10" && O.id(oi)=="7.03"
                        if strcmp(recoveryCase,'reference')
                            extra = 35;
                        elseif strcmp(recoveryCase,'low')
                            extra = min(extra,35);
                        end
                    end
                    elapsed = elapsed + extra;
                    delaySeconds(s) = delaySeconds(s) + extra;
                    contributionIdx = find(recoverableIDs == F.id(r),1);
                    eventDelays(s,contributionIdx) = eventDelays(s,contributionIdx) + extra;
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
    % Store campaign-level statistics in the first ten columns.
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
    % Save one row per successfully installed raw-material element.
    [simIdx,moduleIdx] = find(~isnan(moduleTimes_s));
    positions = sub2ind(size(moduleTimes_s),simIdx,moduleIdx);
    M = table(simIdx,moduleIdx,moduleTimes_s(positions)/60, ...
        'VariableNames',{'simulation','element_number','duration_min'});
    out = fullfile(folder,'results');
    if ~exist(out,'dir'), mkdir(out); end
    writetable(M,fullfile(out,sprintf('elements_%s.csv',scenario)));
    fprintf('%s: 50 elements %.2f h mean; P95 %.2f h; lost %.3f per campaign\n', ...
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
    'fraction_completed','mean_attempts','mean_lost_elements', ...
    'mean_campaign_h','campaign_P95_h','mean_recovery_min_per_module', ...
    'mean_time_min_per_module','baseline_campaign_h', ...
    'mean_individual_element_min','element_P50_min','element_P90_min', ...
    'element_P95_min','element_P99_min','first_element_mean_min', ...
    'subsequent_element_mean_min'});
out = fullfile(folder,'results');
writetable(S,fullfile(out,'summary.csv'));
D = array2table(allDelays,'VariableNames', ...
    cellstr(matlab.lang.makeValidName(recoverableIDs)));
D = addvars(D,string(scenarios(:)),'Before',1,'NewVariableNames','scenario');
writetable(D,fullfile(out,'recovery_contributions.csv'));
disp(S)
figure('Color','w'); bar(categorical(S.scenario), ...
    [S.mean_campaign_h S.campaign_P95_h]);
ylabel('Time to install 50 raw-material elements [h]');
legend('Mean','P95','Location','northwest');
title('FABER — Scenario B campaign time, including replacements');
exportgraphics(gcf,fullfile(out,'campaign_times.png'),'Resolution',180);
figure('Color','w'); bar(categorical(S.scenario),S.mean_lost_elements);
ylabel('Mean elements lost per campaign');
title('FABER — Expected lost elements (exploratory)');
exportgraphics(gcf,fullfile(out,'module_losses.png'),'Resolution',180);

figure('Color','w');
bar(categorical(S.scenario), ...
    [S.mean_individual_element_min S.element_P95_min]);
ylabel('Time per successfully installed element [min]');
legend('Mean','P95','Location','northwest');
title('FABER — Scenario B individual element duration');
exportgraphics(gcf,fullfile(out,'module_times.png'),'Resolution',180);
