function plot_FABER_RPO_finite(R,V,trial,architecture)
% Compare finite-thrust profile and propagated LVLH states with ideal path.
if nargin<3
    [trials,arch]=find(V.passed,1,'first');
    if isempty(trials)
        [trials,arch]=find(V.attempted,1,'first');
    end

    if isempty(trials), warning('No finite validation records.'); return; end
    trial=trials; architecture=arch;
end

D=V.details{trial,architecture};
if isempty(D), warning('No results for selected case.'); return; end
rec=R.records{trial,architecture};

fig=figure('Color','w','Name',sprintf('Finite thrust trial %d arch %d',trial,architecture));
tl=tiledlayout(fig,3,1,'TileSpacing','compact');
for approach=1:2
    if approach==1, s=rec.Client; legName='Client';
    else, s=rec.ATLANTE; legName='ATLANTE'; end
    offset=0;
    for k=1:5
        P=D{approach,k};
        if isempty(P), continue; end
        tmin=(P.time_s+offset)/60;
        nexttile(tl,1); hold on; grid on;
        plot(tmin,vecnorm(P.thrust_LVLH_N,2,2),'DisplayName',sprintf('%s P%d',legName,k));
        nexttile(tl,2); hold on; grid on;
        plot(tmin,P.rho_m,'DisplayName',sprintf('%s P%d',legName,k));
        nexttile(tl,3); hold on; grid on;
        plot(tmin,P.rhoDot_mps,'DisplayName',sprintf('%s P%d',legName,k));
        offset=offset+s.transferDetails{k}.TOF_s;
    end
end

nexttile(tl,1); ylabel('Commanded thrust [N]');
nexttile(tl,2); ylabel('Relative range [m]');
nexttile(tl,3); ylabel('Radial rate [m/s]'); xlabel('Transfer-only elapsed time [min]');
title(tl,sprintf('Finite-thrust simulation | Trial %d | Architecture %d',trial,architecture));
end
