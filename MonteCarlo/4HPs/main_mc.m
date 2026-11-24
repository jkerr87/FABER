% Run the four-holding-point RPO architecture study.
clear; clc; close all;

%% Configuration
cfg = ASV_RPO_MC_config();

cfg.N = 1000;
cfg.seed = 42;
cfg.h_km = 680.12;

%% Run nine-phase comparison
R = FABER_RPO_9phase_MC(cfg);

%% Display all 18 combinations
disp(R.summary);

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
    error('No feasible trials for architecture %d',m);
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
P = FABER_RPO_Pareto(R);

%% Optional: require at least 70% acceptance
P70 = FABER_RPO_Pareto(R,0.70);
