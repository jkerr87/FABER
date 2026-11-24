function [window,tabs] = FABER_environment_plots(out,plan,reference,solar,names,Re)
% Overlay illumination on the reference geometry and show scheduled EPS use.
cfg = out.options;
window = figure('Name','FABER Environment','NumberTitle','off', ...
    'WindowStyle','normal','Visible','off','Position',[100 100 1150 710]);
theme(window,'light');
group = uitabgroup(window,'Units','normalized','Position',[0 0 1 1]);
trajectoryTabs(group,reference,reference,out,solar,cfg,Re,"Reference");
trajectoryTabs(group,plan,reference,out,solar,cfg,Re,"Scheduled");
sun = tabAxes(group,"Sunlight and eclipse","Environment_Sunlight_Eclipse");
angle = tabAxes(group,"Sun aspect angle","Environment_Sun_Aspect");
beta = tabAxes(group,"Orbit beta angle","Environment_Orbit_Beta");

for v=[1 2 3 4]
    x = out.lightingTrace(out.lightingTrace.Vehicle==names(v),:);
    t = hours(x.UTC-cfg.epochUTC);
    plot(sun,t,x.SunFraction,'DisplayName',names(v));
    plot(angle,t,x.SunAspect_deg,'DisplayName',names(v));
    plot(beta,t,x.Beta_deg,'DisplayName',names(v));
end
r = out.referenceLightingTrace(out.referenceLightingTrace.Vehicle=="Taxi",:);
tr = hours(r.UTC-cfg.epochUTC);
plot(sun,tr,r.SunFraction,'--','Color',[0.35 0.35 0.35],'DisplayName','Taxi reference');
plot(angle,tr,r.SunAspect_deg,'--','Color',[0.35 0.35 0.35],'DisplayName','Taxi reference');
ylabel(sun,'Visible solar disk'); ylim(sun,[-0.05 1.05]);
ylabel(angle,'Sun aspect (deg)'); ylim(angle,[0 180]);
ylabel(beta,'Orbit beta (deg)');
for ax = [sun angle beta]
    xlabel(ax,'Hours since UTC epoch');
    legend(ax,'Location','best','NumColumns',2);
end
linkaxes([sun angle beta],'x');

ops = tabAxes(group,"Mission duty cycle","Environment_EPS_Duty_Cycle");
share = tabAxes(group,"Vehicle duty-cycle allocation","Environment_Duty_Cycle_Allocation");
soc = tabAxes(group,"Battery state of charge","Environment_Battery_SoC");
[categories,keys,colors] = timelineCategories(out,names);
edges = out.eps.t(:);
timelinePlot(ops,edges,categories,keys,colors,names,out.openItems,[0 edges(end)]);
for trip = 1:height(out.retrievalSummary)
    r = out.retrievalSummary(trip,:);
    if ~isfinite(r.Start_s), continue; end
    tripEvents = out.events(out.events.Retrieval==trip,:);
    inbound = tripEvents.Activity=="Return Lambert";
    farField = inbound | tripEvents.Activity=="Outbound Lambert";
    start = r.Start_s/60;
    if any(farField), start = max(start,max(tripEvents.End_s(farField))/60-10); end
    stop = min(edges(end),start+cfg.dutyCycleZoom_min);
    if stop<=start, continue; end

    stage = "client RPO";
    if any(inbound), stage = "Atlante RPO and handover"; end
    name = stage;
    tag = "Environment_Duty_Cycle_Zoom";
    if height(out.retrievalSummary)>1
        name = r.Module+" - "+stage;
        tag = tag+"_"+trip;
        if isfinite(r.Handover_s)
            xline(ops,r.Handover_s/3600,':',r.Module+" delivery", ...
                'HandleVisibility','off');
        end
    end
    zoomAx = tabAxes(group,name,tag);
    timelinePlot(zoomAx,edges,categories,keys,colors,names,out.openItems,[start stop]);
    taxiIssues = out.openItems(out.openItems.Vehicle=="Taxi" & ...
        out.openItems.Epoch_s/60>=start & out.openItems.Epoch_s/60<=stop,:);
    if r.Status=="Blocked" && ~isempty(taxiIssues)
        subtitle(zoomAx,taxiIssues.Activity(1)+" - "+taxiIssues.Reason(1),'Interpreter','none');
    end
end
if ~isempty(cfg.assemblyModuleIndices)
    assembly = tabAxes(group,"Kiki and Rango assembly timeline","Environment_Assembly_Timeline");
    issues = out.openItems(startsWith(out.openItems.Vehicle,"Kiki") | ...
        startsWith(out.openItems.Vehicle,"Rango"),:);
    active = startsWith(out.events.Activity,"Assembly ") & ...
        (startsWith(out.events.Vehicle,"Kiki") | startsWith(out.events.Vehicle,"Rango"));
    if any(active)
        first = min(out.events.Start_s(active))/60;
        last = max([out.events.End_s(active);issues.Epoch_s])/60;
    elseif ~isempty(issues)
        first = min(issues.Epoch_s)/60;
        last = max(issues.Epoch_s)/60;
    else
        first = 0;
        last = edges(end);
    end
    limits = [max(0,first-30) min(edges(end),last+30)];
    if limits(2)<=limits(1), limits = [0 edges(end)]; end
    timelinePlot(assembly,edges,categories,keys,colors,names,issues,limits);
    if ~isempty(issues)
        subtitle(assembly,issues.Vehicle(1)+" - "+issues.Activity(1)+" - "+issues.Reason(1), ...
            'Interpreter','none');
    end
end
allocation(share,edges,categories,keys,colors,names);
plot(soc,out.eps.t/60,out.eps.soc,'LineWidth',1.2);
yline(soc,65,'k:','PM4 minimum SoC');
legend(soc,names,'Location','best','NumColumns',3);
ylabel(soc,'Battery SoC at EOL (%)'); xlabel(soc,'Hours since UTC epoch');
linkaxes([ops soc],'x');
tabs = group.Children(:);
group.SelectedTab = tabs(1);
window.Visible = get(groot,'DefaultFigureVisible');
end

function trajectoryTabs(group,replay,reference,out,solar,cfg,Re,kind)
titles = ["Taxi far field (inertial)","Client RPO (LVLH)", ...
    "Atlante RPO (LVLH)","Kiki relocation (Atlante LVLH)"];
tags = ["Taxi_Far_Field","Client_RPO","Atlante_RPO","Kiki_Relocation"];
ax = gobjects(4,1);
for k=1:4
    ax(k) = tabAxes(group,kind+" | "+titles(k),"Environment_"+kind+"_"+tags(k));
    axis(ax(k),'equal');
    view(ax(k),38,24);
    if k==1
        xlabel(ax(k),'X (km)'); ylabel(ax(k),'Y (km)'); zlabel(ax(k),'Z (km)');
    else
        xlabel(ax(k),'Along track (m)'); ylabel(ax(k),'Radial (m)'); zlabel(ax(k),'Cross track (m)');
    end
end
panelTracks(replay,ax,solar,cfg,Re,1:4,'-');
colors = lightColors();
for k=1:4
    if isempty(ax(k).Children) && kind=="Scheduled"
        panelTracks(reference,ax,solar,cfg,Re,k,'--');
        if ~isempty(ax(k).Children)
            subtitle(ax(k),'Reference trajectory (dashed). Scheduled stage blocked or awaiting handover.');
        end
    end
    if isempty(ax(k).Children)
        if k==4 && ~isempty(out.assemblyReferencePaths)
            for j=1:numel(out.assemblyReferencePaths)
                p = out.assemblyReferencePaths{j};
                plot3(ax(k),p(2,:),p(1,:),p(3,:),'--','Color',[0.35 0.35 0.35]);
            end
            subtitle(ax(k),'Phase A geometry only. No scheduled UTC window for lighting.');
        else
            text(ax(k),0.5,0.5,'No Phase A path available for this stage', ...
                'Units','normalized','HorizontalAlignment','center');
        end
        continue;
    end
    keys = gobjects(3,1);
    for i=1:3
        keys(i) = plot3(ax(k),NaN,NaN,NaN,'Color',colors(i,:),'LineWidth',2);
    end
    legend(ax(k),keys,{'Sunlight','Penumbra','Umbra'},'Location','best');
end
end

function ax = tabAxes(group,name,tag)
tab = uitab(group,'Title',name,'Tag',tag);
ax = axes('Parent',tab,'Units','normalized');
hold(ax,'on');
grid(ax,'on');
title(ax,name,'Interpreter','none');
end

function panelTracks(plan,axesList,solar,cfg,Re,include,style)
colors = lightColors();
for i=1:numel(plan.paths)
    p = plan.paths{i};
    if p.stop<=p.start, continue; end
    if p.vehicle=="Taxi"
        switch p.leg
            case "Far field", k=1;
            case "Client", k=2;
            case "Atlante", k=3;
            otherwise, continue;
        end
    elseif startsWith(p.vehicle,"Kiki") && startsWith(p.leg,"Assembly")
        k=4;
    else
        continue;
    end
    if ~ismember(k,include), continue; end
    t = linspace(p.start,p.stop,max(2,ceil((p.stop-p.start)/cfg.sampleStep_s)+1));
    x = p.state(t);
    f = FABER_environment_light(solar,(t(1:end-1)+t(2:end))/2, ...
        p.state((t(1:end-1)+t(2:end))/2),cfg.aspectAxisLVLH,Re);
    state = ones(size(f));
    state(f<1-1e-8)=2;
    state(f<=1e-9)=3;
    if k==1
        xyz = x(1:3,:)/1000;
    else
        chief = plan.atlante(t);
        if k==2, chief=p.clientChief(t); end
        r = chief(1:3,:);
        radial = r./vecnorm(r);
        normal = cross(r,chief(4:6,:),1);
        normal = normal./vecnorm(normal);
        transverse = cross(normal,radial,1);
        delta = x(1:3,:)-r;
        xyz = [sum(transverse.*delta,1);sum(radial.*delta,1);sum(normal.*delta,1)];
    end
    starts = [1 find(diff(state)~=0)+1];
    ends = [starts(2:end) numel(t)];
    for j=1:numel(starts)
        q = starts(j):ends(j);
        plot3(axesList(k),xyz(1,q),xyz(2,q),xyz(3,q), ...
            'Color',colors(state(starts(j)),:),'LineStyle',style,'LineWidth',1.6);
    end
end
end

function [category,keys,colors] = timelineCategories(out,names)
% Assign each interval exactly once per vehicle using the executed timeline.
keys = ["Standby","Far-field transit","Sunlit RPO","Eclipse hold","Other hold", ...
    "Loading / handover","Capture / installation","Charging","Refueling","Open state"];
colors = [0.82 0.84 0.86;0.19 0.43 0.67;0.89 0.57 0.10;0.24 0.26 0.48; ...
    0.43 0.50 0.67;0.55 0.36 0.67;0.15 0.56 0.54;0.29 0.62 0.36; ...
    0.78 0.42 0.25;0.70 0.29 0.37];
category = ones(size(out.eps.mode));
mid = (out.eps.t(1:end-1)+out.eps.t(2:end))/2;
for i=1:height(out.events)
    e = out.events(i,:);
    v = find(names==e.Vehicle,1);
    k = mid>=e.Start_s/60 & mid<e.End_s/60;
    name = lower(e.Activity);
    code = 1;
    if e.SunRequired, code=3; end
    if any(e.Mode==["T1","T4"]), code=2; end
    if any(e.Mode==["T3","R4"]), code=7; end
    if any(e.Mode==["T6","K3","K5","R2"]) || ...
            contains(name,"deploy rango") || contains(name,"recover rango")
        code=6;
    end
    if any(e.Mode==["K6","K7","R5"]) || contains(name,"top-up"), code=8; end
    if any(e.Mode==["T9","T9b","A4","K8"]), code=9; end
    if startsWith(name,"attached wait") || startsWith(name,"eclipse wait") || ...
            contains(name,"hold") || name=="attached capture dwell"
        code=5;
    end
    if startsWith(name,"post-stop state open"), code=10; end
    category(k,v)=code;
    if code==5
        light = out.lightingTrace.SunFraction(out.lightingTrace.Vehicle==names(v));
        category(k & light<1-1e-8,v)=4;
    end
end
end

function timelinePlot(ax,edges,category,keys,colors,names,issues,window)
for v=1:numel(names)
    modes = category(:,v);
    if isempty(modes), continue; end
    starts = [1;find(diff(modes)~=0)+1];
    stops = [starts(2:end);numel(edges)];
    for j=1:numel(starts)
        a = max(window(1),edges(starts(j)));
        b = min(window(2),edges(stops(j)));
        if b>a
            patch(ax,[a b b a]/60,v+[-.32 -.32 .32 .32], ...
                colors(modes(starts(j)),:),'EdgeColor','none','HandleVisibility','off');
        end
    end
end
handle = gobjects(numel(keys)+1,1);
for k=1:numel(keys)
    handle(k) = patch(ax,nan,nan,colors(k,:),'EdgeColor','none');
end
for j=1:height(issues)
    v = find(names==issues.Vehicle(j),1);
    t = issues.Epoch_s(j)/60;
    if ~isempty(v) && t>=window(1) && t<=window(2)
        plot(ax,t/60,v,'kd','MarkerFaceColor','w','HandleVisibility','off');
    end
end
handle(end) = plot(ax,nan,nan,'kd','MarkerFaceColor','w');
legend(ax,handle,[keys,"Blocked or open item"],'Location','southoutside','NumColumns',3);
set(ax,'YDir','reverse'); yticks(ax,1:numel(names)); yticklabels(ax,names);
ylim(ax,[0.5 numel(names)+0.5]);
if window(2)>window(1), xlim(ax,window/60); end
xlabel(ax,'Hours since UTC epoch');
end

function allocation(ax,edges,category,keys,colors,names)
elapsed = diff(edges);
duration = zeros(numel(names),numel(keys));
if isempty(elapsed), return; end
for v=1:numel(names)
    for k=1:numel(keys)
        duration(v,k)=100*sum(elapsed(category(:,v)==k))/sum(elapsed);
    end
end
bars = barh(ax,1:numel(names),duration,'stacked','EdgeColor','none');
for k=1:numel(keys), bars(k).FaceColor=colors(k,:); end
legend(ax,bars,keys,'Location','southoutside','NumColumns',3);
set(ax,'YDir','reverse'); yticks(ax,1:numel(names)); yticklabels(ax,names);
xlim(ax,[0 100]); xlabel(ax,'Share of elapsed schedule (%)');
end

function color = lightColors()
color = [0.91 0.58 0.10; 0.59 0.42 0.72; 0.16 0.31 0.52];
end
