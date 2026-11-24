function P = FABER_RPO_Pareto(R,minAcceptance)
% FABER_RPO_PARETO
%
% Architecture-level P95 screening and
% scenario-by-scenario discrete Pareto comparison.
%
% Lower TOF and lower Delta-V are preferred.
%
% Usage:
%   P = FABER_RPO_Pareto(R)
%   P = FABER_RPO_Pareto(R,0.70)

if nargin<2
    minAcceptance = 0;
end

T = R.summary;

N = size(R.valid,1);
M = size(R.valid,2);

%% 1. ARCHITECTURE-LEVEL P95 PARETO

x = T.P95CycleTOF_min;
y = T.P95CycleDV_mps;

eligible = T.Acceptance>=minAcceptance & ...
    T.Accepted>0 & isfinite(x) & isfinite(y);

paretoP95 = false(M,1);

for m = find(eligible)'

    dominates = eligible & ...
        x<=x(m)+1e-9 & ...
        y<=y(m)+1e-9 & ...
        (x<x(m)-1e-9 | y<y(m)-1e-9);

    paretoP95(m) = ~any(dominates);
end

%% 2. SCENARIO-BY-SCENARIO PARETO

scenarioPareto = false(N,M);

for i = 1:N

    idx = find(R.valid(i,:));

    for m = idx

        dv = R.CycleDV_mps(i,m);
        tof = R.CycleTOF_min(i,m);

        dominated = false;

        for j = idx

            if j==m
                continue;
            end

            dv2 = R.CycleDV_mps(i,j);
            tof2 = R.CycleTOF_min(i,j);

            if dv2<=dv+1e-9 && ...
               tof2<=tof+1e-9 && ...
               (dv2<dv-1e-9 || tof2<tof-1e-9)

                dominated = true;
                break;
            end
        end

        scenarioPareto(i,m) = ~dominated;
    end
end

ParetoCount = sum(scenarioPareto,1)';
ParetoFractionOfFeasible = nan(M,1);

ok = T.Accepted>0;

ParetoFractionOfFeasible(ok) = ...
    ParetoCount(ok)./T.Accepted(ok);

%% 3. RESULTS

P = struct();

P.P95Mask = paretoP95;
P.ScenarioMask = scenarioPareto;

P.table = T;

P.table.ParetoP95 = paretoP95;
P.table.ParetoCount = ParetoCount;
P.table.ParetoFractionOfFeasible = ...
    ParetoFractionOfFeasible;

disp(P.table(P.table.ParetoP95,:));

%% 4. P95 SCATTER PLOT

figure( ...
    'WindowStyle','docked', ...
    'Color','w', ...
    'Name','FABER Delta-V vs TOF Pareto');

hold on;
grid on;

for m = 1:M

    if ~eligible(m)
        continue;
    end

    if strcmpi(T.gateMode(m),'FullStop')
        marker = 'o';
    else
        marker = 's';
    end

    if paretoP95(m)
        scatter(x(m),y(m),95, ...
            marker,'filled', ...
            'MarkerEdgeColor','k');
    else
        scatter(x(m),y(m),35, ...
            marker, ...
            'MarkerEdgeColor',[0.6 0.6 0.6]);
    end
end

% Sort nondominated P95 points by TOF.
idx = find(paretoP95);

[~,order] = sort(x(idx));
idx = idx(order);

if ~isempty(idx)

    plot(x(idx),y(idx),'k--', ...
        'LineWidth',1.3, ...
        'HandleVisibility','off');

    for k = 1:numel(idx)

        m = idx(k);

        text(x(m),y(m), ...
            sprintf('  %d',m), ...
            'FontWeight','bold');
    end
end

xlabel('P95 complete-cycle TOF [min]');
ylabel('P95 complete-cycle Delta-V [m/s]');

title(sprintf( ...
    'FABER architecture trade study (acceptance >= %.0f%%)', ...
    100*minAcceptance));

%% 5. PRINT SCENARIO-LEVEL FREQUENCIES

fprintf('\nSCENARIO-LEVEL PARETO FREQUENCIES\n');

disp(P.table(:,{ ...
    'methodH2','methodC1','methodC2', ...
    'gateMode','Accepted', ...
    'ParetoCount','ParetoFractionOfFeasible'}));

end
