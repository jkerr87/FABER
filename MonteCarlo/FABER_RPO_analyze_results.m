function A = FABER_RPO_analyze_results(MC)
% Analyse SELECTED FEASIBLE designs from the existing Monte Carlo.
% This does NOT establish raw-design feasibility by count; use
% FABER_RPO_count_feasibility for that separate experiment.
T = MC.trialTable;
D = MC.DV_phase_mps;
F = MC.TOF_phase_min;
assert(size(D,2)==10 && size(F,2)==10,'Expected ten phase columns.');
assert(all(ismember({'ClientN_Homing','ClientN_Closing','HostN_Homing', ...
 'HostN_Closing'},T.Properties.VariableNames)), 'Transfer-count columns missing.');
labels = {'Client Homing','Client Closing','Host Homing','Host Closing'};
countFields = {'ClientN_Homing','ClientN_Closing','HostN_Homing','HostN_Closing'};
phaseCols = [1 3 6 8];
rows = cell(0,7);
for j=1:4
 for k=1:4
  ix = T.(countFields{j})==k;
  x = D(ix,phaseCols(j)); y = F(ix,phaseCols(j));
  if isempty(x)
   rows(end+1,:)={string(labels{j}),k,0,NaN,NaN,NaN,NaN}; %#ok<AGROW>
  else
   rows(end+1,:)={string(labels{j}),k,numel(x),mean(x),pct(x,95),mean(y),pct(y,95)}; %#ok<AGROW>
  end
 end
end
A.countTable = cell2table(rows,'VariableNames', ...
 {'Phase','NTransfers','NSelected','MeanDV_mps','P95DV_mps','MeanTOF_min','P95TOF_min'});
disp(A.countTable);
if ~exist(MC.cfg.outputFolder,'dir'), mkdir(MC.cfg.outputFolder); end
writetable(A.countTable,fullfile(MC.cfg.outputFolder,'RPO_MC_selected_counts.csv'));

% True non-dominated frontier among the SELECTED feasible client+host cycles.
% Sort ascending DV and retain only strictly decreasing TOF.
dv = T.CycleDV_mps(:); tof = T.CycleTOF_min(:);
[~,ord]=sortrows([dv,tof],[1 2]);
keep=false(numel(dv),1); bestTOF=inf;
for j=1:numel(ord)
 i=ord(j);
 if tof(i)<bestTOF-1e-9
  keep(i)=true;bestTOF=tof(i);
 end
end
A.paretoTable = T(keep,:);
A.paretoTable=sortrows(A.paretoTable,'CycleTOF_min');
writetable(A.paretoTable,fullfile(MC.cfg.outputFolder,'RPO_MC_Pareto.csv'));
f=figure('Color','w','Name','FABER RPO Pareto frontier');
scatter(tof/60,dv,18,[0.65 0.65 0.65],'filled');hold on;
plot(A.paretoTable.CycleTOF_min/60,A.paretoTable.CycleDV_mps, ...
 'o-','Color',[0 0.4 0.75],'LineWidth',1.7,'MarkerFaceColor',[0 0.4 0.75]);
[~,iMinDV]=min(dv);[~,iMinTOF]=min(tof);
plot(tof(iMinDV)/60,dv(iMinDV),'kp','MarkerSize',12,'MarkerFaceColor','k');
plot(tof(iMinTOF)/60,dv(iMinTOF),'ks','MarkerSize',9,'LineWidth',2);
xlabel('Client + Host RPO TOF [hours]');ylabel('Client + Host RPO Delta V [m/s]');
title('Pareto frontier of selected feasible retrieval-cycle designs');grid on;
legend('Selected feasible cycles','Non-dominated solutions','Minimum DV','Minimum TOF', ...
 'Location','best');
saveas(f,fullfile(MC.cfg.outputFolder,'RPO_MC_Pareto.png'));

f=figure('Color','w','Name','RPO manoeuvre count vs phase DV');
tiledlayout(2,2);
for j=1:4
 nexttile;
 scatter(T.(countFields{j}),D(:,phaseCols(j)),20,'filled');
 xlabel('Transfers per phase');ylabel('Phase Delta V [m/s]');
 title(labels{j});xlim([0.5 4.5]);xticks(1:4);grid on;
end
saveas(f,fullfile(MC.cfg.outputFolder,'RPO_MC_counts_scatter.png'));
fprintf('\nPareto points: %d of %d selected feasible cycles.\n',height(A.paretoTable),height(T));
fprintf('Selected-count CSV is NOT a feasibility analysis; run FABER_RPO_count_feasibility.\n');
end
function p=pct(x,q)
x=sort(x(:));n=numel(x);if n==1,p=x;return;end
u=1+(n-1)*q/100;a=floor(u);b=ceil(u);p=x(a)+(u-a)*(x(b)-x(a));
end
