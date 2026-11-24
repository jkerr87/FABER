function plot_FABER_mission_MC(R,B,architecture)
% Plot one actual accepted paired cycle near the conditional cycle DV P95.
valid = R.valid(:,architecture);
dv = R.CycleDV_mps(:,architecture);
figure('Name','RPO Monte Carlo screening','Color','w');
tiledlayout(1,2);
nexttile;
bar(categorical(["Accepted","Rejected"]),[nnz(valid),nnz(~valid)]);
ylabel('Trial count');
title('Paired RPO cycle screening');
nexttile;
if any(valid)
    scatter(R.CycleTOF_min(valid,architecture),dv(valid),18,'filled');
end

xlabel('Close-field cycle time [min]');
ylabel('Close-field cycle delta-V [m/s]');
title(sprintf('Feasible campaign P95 known subtotal: %.3g m/s',B.P95KnownRawDV_mps));
grid on;
if ~any(valid)
    return;
end

% Select an existing trial so its outbound and return traces stay paired.
accepted = find(valid);
[~,nearest] = min(abs(dv(valid)-prctile(dv(valid),95)));
trial = R.records{accepted(nearest),architecture};
approaches = {trial.Client,trial.ATLANTE};
names = {'Client','Atlante'};
figure('Name','RPO radial and along-track projection','Color','w');
tiledlayout(1,2);
phases = [1 3 5 7 9];
angle = linspace(0,2*pi,361);
for leg = 1:2
    ax = nexttile;
    hold(ax,'on');
    plot(ax,R.cfg.AZ_m*cos(angle),R.cfg.AZ_m*sin(angle), ...
        'k--','DisplayName','AZ');
    plot(ax,R.cfg.KOZ_m*cos(angle),R.cfg.KOZ_m*sin(angle), ...
        'k:','DisplayName','KOZ');
    for phase = phases
        X = approaches{leg}.trace{phase};
        plot(ax,X(1,:),X(2,:),'DisplayName',R.phaseNames{phase});
    end

    for k = 1:4
        X = approaches{leg}.trace{2*k};
        plot(ax,X(1,1),X(2,1),'ko','HandleVisibility','off');
        text(ax,X(1,1),X(2,1),sprintf(' HP%d',k));
    end

    plot(ax,0,0,'kp','DisplayName','Target');
    axis(ax,'equal');
    grid(ax,'on');
    xlabel(ax,'Radial [m]');
    ylabel(ax,'Along-track [m]');
    title(ax,names{leg});
    legend(ax,'Location','best');
end

figure('Name','Representative RPO cycle','Color','w');
tiledlayout(3,2);
for leg = 1:2
    legData = approaches{leg};
    ax = nexttile(leg);
    hold(ax,'on');
    rangeAx = nexttile(leg+2);
    hold(rangeAx,'on');
    speedAx = nexttile(leg+4);
    hold(speedAx,'on');
    elapsed = 0;
    colors = lines(numel(legData.trace));
    for phase = 1:numel(legData.trace)
        X = legData.trace{phase};
        t = elapsed+linspace(0,legData.tofPhase(phase),size(X,2));
        elapsed = t(end);
        r = X(1:3,:);
        range = vecnorm(r);
        closing = -sum(r.*X(4:6,:),1)./max(range,eps);
        plot3(ax,r(1,:),r(2,:),r(3,:),'Color',colors(phase,:));
        plot(rangeAx,t,range,'Color',colors(phase,:));
        plot(speedAx,t,closing,'Color',colors(phase,:));
    end

    view(ax,3);
    axis(ax,'equal');
    grid(ax,'on');
    title(ax,names{leg});
    xlabel(ax,'Radial [m]');
    ylabel(ax,'Along-track [m]');
    zlabel(ax,'Cross-track [m]');
    yline(rangeAx,R.cfg.AZ_m,':','AZ');
    yline(rangeAx,R.cfg.KOZ_m,':','KOZ');
    ylabel(rangeAx,'Range [m]');
    xlabel(rangeAx,'Approach time [min]');
    ylabel(speedAx,'Radial closing speed [m/s]');
    xlabel(speedAx,'Approach time [min]');
    grid(rangeAx,'on');
    grid(speedAx,'on');
end
end
