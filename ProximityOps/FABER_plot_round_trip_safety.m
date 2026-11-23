function figs = FABER_plot_round_trip_safety(out)
% LVLH ensembles, gate dispersions, passive drift and sensitivity for ASMAD.
T = out.trials;
records = out.records;
cfg = out.options;
cases = out.summary.Case;
colors = [0.13 0.35 0.62; 0.79 0.33 0.12; 0.17 0.56 0.43; ...
    0.59 0.32 0.68; 0.66 0.45 0.04; 0.29 0.52 0.65; ...
    0.68 0.21 0.39; 0.33 0.45 0.22; 0.42 0.42 0.42];
colors = colors(1:numel(cases),:);
figs = gobjects(0,1);
figs(end+1) = figure('Name','Round trip LVLH ensembles','Tag','Round_Trip_LVLH','Color','w','WindowStyle','docked');
theme(figs(end),'light');
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
farOut = nexttile; farBack = nexttile; clientAx = nexttile; hostAx = nexttile;
for ax = [farOut farBack clientAx hostAx]
    hold(ax,'on'); grid(ax,'on'); axis(ax,'equal');
    xlabel(ax,'Along-track'); ylabel(ax,'Radial');
end
title(farOut,'Outbound far field (km)'); title(farBack,'Return far field (km)');
title(clientAx,'Client RPO (m)'); title(hostAx,'Atlante RPO (m)');
circle(farOut,cfg.farProtectedRadius_m/1000,':');
circle(farBack,cfg.farProtectedRadius_m/1000,':');
rpo = ASV_RPO_MC_config();
for key = ["AZ_m","KOZ_m"]
    if isfield(cfg.rpo,key)
        rpo.(key) = cfg.rpo.(key);
    end
end
for ax = [clientAx hostAx]
    circle(ax,rpo.AZ_m,'--');
    circle(ax,rpo.KOZ_m,':');
    plot(ax,0,0,'k+','HandleVisibility','off');
end
if isfield(cfg.rpo,'terminalCorridorHalfWidth_m')
    w = cfg.rpo.terminalCorridorHalfWidth_m;
    if isfinite(w)
        for ax = [clientAx hostAx]
            plot(ax,[-rpo.KOZ_m 0],[w w],'k:','HandleVisibility','off');
            plot(ax,[-rpo.KOZ_m 0],[-w -w],'k:','HandleVisibility','off');
        end
    end
end

figs(end+1) = figure('Name','Range and closing speed','Tag','Round_Trip_Range_and_Closing','Color','w','WindowStyle','docked');
theme(figs(end),'light');
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
farRange = nexttile; passiveRange = nexttile; rpoRange = nexttile; closing = nexttile;
for ax = [farRange passiveRange rpoRange closing]
    hold(ax,'on'); grid(ax,'on');
end
title(farRange,'Far-field separation'); xlabel(farRange,'Hours since departure');
ylabel(farRange,'Range (m)'); set(farRange,'YScale','log');
yline(farRange,cfg.farProtectedRadius_m,'k:','Protected radius');
title(passiveRange,'Passive drift after fault or rejected gate');
xlabel(passiveRange,'Minutes since passive drift began'); ylabel(passiveRange,'Range (m)');
set(passiveRange,'YScale','log');
yline(passiveRange,cfg.farProtectedRadius_m,'k:','Protected radius');
yline(passiveRange,cfg.contactRadius_m,'k--','Contact radius');
title(rpoRange,'RPO range'); xlabel(rpoRange,'Minutes since RPO entry'); ylabel(rpoRange,'Range (m)');
yline(rpoRange,rpo.AZ_m,'k--','AZ'); yline(rpoRange,rpo.KOZ_m,'k:','KOZ');
title(closing,'RPO radial closing speed');
xlabel(closing,'Minutes since RPO entry'); ylabel(closing,'Closing speed (m/s)');
yline(closing,0,'k:');

legendHandles = gobjects(numel(cases),1);
for j = 1:numel(cases)
    ids = find(T.Case==cases(j));
    ids = ids(unique(round(linspace(1,numel(ids),min(numel(ids),cfg.plotLimit)))));
    color = colors(j,:);
    legendHandles(j) = plot(farOut,NaN,NaN,'Color',color,'DisplayName',cases(j));
    for i = ids'
        v = records{i};
        drawArc(farOut,farRange,v.outbound,color,'-');
        drawArc(farBack,farRange,v.returnArc,color,'--');
        drawApproach(clientAx,rpoRange,closing,v.client,color,'-');
        drawApproach(hostAx,rpoRange,closing,v.host,color,'--');
        drawDrift(passiveRange,v.outbound.passiveDrift,color,'-');
        drawDrift(passiveRange,v.captureDrift,color,':');
        if ~isempty(v.returnArc)
            drawDrift(passiveRange,v.returnArc.passiveDrift,color,'--');
        end
    end
end
lg = legend(farOut,legendHandles,'Location','best','Interpreter','none');
lg.FontSize = 8;
subtitle(rpoRange,'Solid client. Dashed Atlante.');
subtitle(farRange,'Solid outbound. Dashed return.');

figs(end+1) = figure('Name','Round trip gates and distributions','Tag','Round_Trip_Dispersions','Color','w','WindowStyle','docked');
theme(figs(end),'light');
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
gateClient = nexttile; gateHost = nexttile;
for ax = [gateClient gateHost]
    hold(ax,'on'); grid(ax,'on');
    xlabel(ax,'Position error (m)'); ylabel(ax,'Velocity error (m/s)');
    xline(ax,cfg.gatePositionTolerance_m,'k:'); yline(ax,cfg.gateVelocityTolerance_mps,'k:');
end
for j = 1:numel(cases)
    keep = T.Case==cases(j);
    scatter(gateClient,T.ClientGatePositionError_m(keep),T.ClientGateVelocityError_mps(keep), ...
        20,colors(j,:),'filled');
    scatter(gateHost,T.AtlanteGatePositionError_m(keep),T.AtlanteGateVelocityError_mps(keep), ...
        20,colors(j,:),'filled');
end
title(gateClient,'Client RPO entry'); title(gateHost,'Atlante RPO entry');
ax = nexttile;
sx = out.summary;
counts = [sx.OutboundGateAccepted,sx.ClientRPOAccepted,sx.Captures, ...
    sx.ReturnGateAccepted,sx.Completed]./sx.Trials;
bar(ax,1:height(sx),counts,'grouped'); ylim(ax,[0 1]); grid(ax,'on');
xticks(ax,1:height(sx)); xticklabels(ax,sx.Case); xtickangle(ax,30);
title(ax,'Fraction reaching each stage');
legend(ax,{'Client gate','Client RPO','Capture','Atlante gate','Complete'},'Location','best');
distribution(nexttile,T.CompleteDV_mps,'Complete round trips','Total ΔV (m/s)');
distribution(nexttile,T.CompleteDuration_hr,'Complete round trips','Elapsed time (h)');
distribution(nexttile,T.PassiveClosestApproach_m,'Screened passive coasts','Closest range (m)');

figs(end+1) = figure('Name','Round trip dispersion sensitivity','Tag','Round_Trip_Sensitivity','Color','w','WindowStyle','docked');
theme(figs(end),'light');
ax = axes(figs(end));
S = out.sensitivity;
inputs = unique(S.Input,'stable');
metrics = unique(S.Metric,'stable');
values = reshape(S.SpearmanRho,numel(inputs),numel(metrics));
im = imagesc(ax,values,[-1 1]);
im.AlphaData = isfinite(values);
ax.Color = [0.91 0.91 0.91];
colorbar(ax);
xticks(ax,1:numel(metrics)); xticklabels(ax,metrics); xtickangle(ax,20);
yticks(ax,1:numel(inputs)); yticklabels(ax,inputs);
title(ax,'Spearman rank correlation for dispersed trials');
subtitle(ax,'Gray cells indicate insufficient variation or sample size.');

figs(end+1) = figure('Name','Lambert time trade and terminal paths', ...
    'Tag','Round_Trip_Transfer_Trade','Color','w','WindowStyle','docked');
theme(figs(end),'light');
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
trade(nexttile,out.nominalOutbound,'Outbound nominal Lambert');
returnTrade = nexttile;
selectedReturn = find(T.Complete,1);
if isempty(selectedReturn)
    selectedReturn = find(cellfun(@(x) ~isempty(x.returnPlan),records),1);
end
if ~isempty(selectedReturn)
    trade(returnTrade,records{selectedReturn}.returnPlan, ...
        sprintf('Return Lambert | trial %d',selectedReturn));
else
    text(returnTrade,0.5,0.5,'No return arc targeted','Units','normalized', ...
        'HorizontalAlignment','center');
end
terminalClient = nexttile; terminalHost = nexttile;
for ax = [terminalClient terminalHost]
    hold(ax,'on'); grid(ax,'on'); axis(ax,'equal'); view(ax,34,24);
    xlabel(ax,'Along-track (m)'); ylabel(ax,'Radial (m)'); zlabel(ax,'Cross-track (m)');
    drawTerminalBoundary(ax,rpo.KOZ_m,cfg);
end
title(terminalClient,'Client terminal approach'); title(terminalHost,'Atlante terminal approach');
for j = 1:numel(cases)
    ids = find(T.Case==cases(j));
    ids = ids(unique(round(linspace(1,numel(ids),min(numel(ids),cfg.plotLimit)))));
    for i = ids'
        terminalPath(terminalClient,records{i}.client,colors(j,:));
        terminalPath(terminalHost,records{i}.host,colors(j,:));
    end
end

if ~isempty(out.rpoSafetyClient) || ~isempty(out.rpoSafetyAtlante)
    figs(end+1) = figure('Name','Conditional RPO fault screening','Tag','RPO_Fault_Screening','Color','w','WindowStyle','docked');
    theme(figs(end),'light');
    tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    results = {out.rpoSafetyClient,out.rpoSafetyAtlante};
    target = ["Client","Atlante"];
    for j=1:2
        r = results{j};
        clearances = []; cam = [];
        if ~isempty(r)
            t = r.trials;
            clearances = t.PassiveMinClearance_m;
            cam = t.CAM_DV_mps(t.CAMRequired & t.CAMFeasible);
        end
        distribution(nexttile,clearances,target(j)+" passive fault drift",'Minimum clearance (m)');
        distribution(nexttile,cam,'Feasible CAM cases that required CAM','CAM ΔV (m/s)');
    end
end
end

function drawArc(ax,rangeAx,arc,color,style)
if isempty(arc), return; end
x = arc.trajectoryLVLH;
plot(ax,x(2,:)/1000,x(1,:)/1000,style,'Color',color,'HandleVisibility','off');
plot(rangeAx,arc.time_s/3600,max(eps,vecnorm(x(1:3,:))),style,'Color',color);
end

function drawDrift(ax,arc,color,style)
if ~isempty(arc)
    plot(ax,arc.time_s/60,max(eps,vecnorm(arc.trajectoryLVLH(1:3,:))),style,'Color',color);
end
end

function drawApproach(orbit,rangeAx,closingAx,approach,color,style)
if isempty(approach), return; end
elapsed = 0;
for k=1:numel(approach.trace)
    x = approach.trace{k};
    duration = approach.tofPhase(k);
    t = elapsed+linspace(0,duration,size(x,2));
    rho = vecnorm(x(1:3,:));
    plot(orbit,x(2,:),x(1,:),style,'Color',color);
    plot(rangeAx,t,rho,style,'Color',color);
    closing = -sum(x(1:3,:).*x(4:6,:),1)./max(rho,eps);
    plot(closingAx,t,closing,style,'Color',color);
    elapsed = elapsed+duration;
end
end

function circle(ax,r,style)
theta = linspace(0,2*pi,145);
plot(ax,r*cos(theta),r*sin(theta),style,'Color',[0.25 0.25 0.25],'HandleVisibility','off');
end

function distribution(ax,values,name,label)
values = values(isfinite(values));
if isempty(values)
    text(ax,0.5,0.5,'No evaluated cases','Units','normalized','HorizontalAlignment','center');
else
    histogram(ax,values,'FaceColor',[0.17 0.43 0.63]);
    subtitle(ax,sprintf('n = %d, P95 = %.3g',numel(values),prctile(values,95)));
end
title(ax,name);
xlabel(ax,label); ylabel(ax,'Trials'); grid(ax,'on');
end

function trade(ax,arc,name)
s = arc.search;
hold(ax,'on'); grid(ax,'on');
palette = [0.2 0.43 0.75; 0.18 0.58 0.40; 0.12 0.12 0.12];
for branch = unique(s.Branch)'
    x = sortrows(s(s.Branch==branch & s.PassesPerigee,:),'Duration_s');
    style = '-o';
    if branch==1, style = '--o'; end
    plot(ax,x.Duration_s/3600,x.DepartureDV_mps,style,'Color',palette(1,:),'MarkerSize',3);
    plot(ax,x.Duration_s/3600,x.ArrivalDV_mps,style,'Color',palette(2,:),'MarkerSize',3);
    plot(ax,x.Duration_s/3600,x.TotalDV_mps,style,'Color',palette(3,:),'MarkerSize',3);
end
selectedDV = norm(arc.departureBurn_mps)+norm(arc.arrivalBurn_mps);
scatter(ax,arc.transferTime_s/3600,selectedDV,45,[0.77 0.17 0.15],'filled');
xlabel(ax,'Transfer duration (h)'); ylabel(ax,'Delta-V (m/s)');
title(ax,name);
subtitle(ax,'Blue departure, green arrival, black total, red selected');
end

function drawTerminalBoundary(ax,koz,cfg)
[sx,sy,sz] = sphere(22);
surf(ax,koz*sx,koz*sy,koz*sz,'FaceColor',[0.7 0.77 0.84], ...
    'FaceAlpha',0.12,'EdgeColor','none');
if isfield(cfg.rpo,'terminalCorridorHalfWidth_m')
    w = cfg.rpo.terminalCorridorHalfWidth_m;
    if isfinite(w)
        [y,z,x] = cylinder(w,24);
        surf(ax,-koz*x,y,z,'FaceColor',[0.8 0.72 0.24], ...
            'FaceAlpha',0.15,'EdgeColor','none');
    end
end
plot3(ax,0,0,0,'k+');
end

function terminalPath(ax,approach,color)
if ~isempty(approach)
    x = approach.trace{end};
    plot3(ax,x(2,:),x(1,:),x(3,:),'-','Color',color);
end
end
