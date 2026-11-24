function [metrics,cycles,histogram] = eps_battery_lifetime(c,h,r,label,fleetPass)
% Energy-based cycle analysis of an actual constrained dispatch.
% No Signal Processing Toolbox required. No life prediction without test data.
% EFC and DoD use energy SoC, not measured cell Ah SoC.
p=c.life; names=c.names'; n=numel(names); dt=r.dt_h(:); hours=sum(dt);
assert(hours>0 && all(dt>0));
assert(all(isfinite(r.soc),'all') && all(r.soc>=-1e-7 & r.soc<=100+1e-7,'all'));
missionHours=p.years*p.daysPerYear*24;
scale=missionHours*p.journeyFraction/hours; % exposure rate, not chronological replay
unmodelledHours=missionHours*(1-p.journeyFraction);
cycles=table(); histogram=table(); metrics=table();
for v=1:n
    x=r.soc(:,v)/100;
    cap=r.EOL(v); bol=h.Battery_kWh(v);
    % Battery terminal energy proxy uses existing aggregate path efficiency.
    dis=sum(r.batOut(:,v).*dt)/c.batteryDischargeEfficiency;
    ch=sum(r.batIn(:,v).*dt)*c.batteryChargeEfficiency;
    assert(abs(ch-dis-(r.E(end,v)-r.E(1,v)))<1e-7,'Battery throughput ledger mismatch.');
    rf=eps_rainflow(x);
    count=rf(:,1); depth=rf(:,2); meanSoC=rf(:,3);
    efcRF=sum(count.*depth);
    % Range-weighted rainflow equals half total variation, including half cycles.
    assert(abs(efcRF-(ch+dis)/(2*cap))<1e-7,'Rainflow throughput mismatch.');
    maxDepth=max([0;depth]); weightedDepth=NaN;
    if sum(count)>0, weightedDepth=sum(count.*depth)/sum(count); end
    high=hoursAbove(x,dt,p.highSoC);
    meanX=sum((x(1:end-1)+x(2:end))/2.*dt)/hours;
    over=sum(count(depth>p.referenceDoD(v)+1e-9));
    % Exposure extrapolation is deliberately withheld for failed fleet traces.
    % A clipped trace can underestimate required lifetime throughput.
    sc=scale; if ~fleetPass, sc=NaN; end
    cycleProjection=sum(count)*sc;
    refCount=NaN; if v<=2, refCount=60000; end
    z=table(names(v),label,c.batteryEOL,c.maxDoD,hours,bol,cap,100*min(x), ...
        100*maxDepth,100*weightedDepth,sum(count),over,100*p.referenceDoD(v), ...
        refCount,dis,ch,dis/cap,dis/bol,efcRF,100*meanX,100*high/hours, ...
        100*(x(end)-x(1)),max(r.batOut(:,v))/c.batteryDischargeEfficiency/cap, ...
        max(r.batIn(:,v))*c.batteryChargeEfficiency/cap, ...
        cycleProjection,over*sc,dis/bol*sc,dis/cap*sc,high*sc,unmodelledHours,fleetPass, ...
        'VariableNames',{'Vehicle','Case','Retention','DoDLimit','Trace_h','AvailableBOL_kWh', ...
        'AvailableAtRetention_kWh','MinSoC_pct','MaxCycleDoD_pct','WeightedMeanCycleDoD_pct', ...
        'RainflowCycleCount','CyclesAboveReferenceDoD','SupplierReferenceDoD_pct', ...
        'SupplierReferenceCycleCount','StoredDischarged_kWh','StoredCharged_kWh', ...
        'EFC_DischargeAtRetention','EFC_DischargeBOL','EFC_RainflowAtRetention', ...
        'TimeWeightedMeanSoC_pct','TimeAboveHighSoC_pct','EndMinusStartSoC_pp', ...
        'MaxStoredDischargeRate_per_h','MaxStoredChargeRate_per_h', ...
        'ProjectedCycles','ProjectedCyclesAboveReferenceDoD','ProjectedEFC_BOL', ...
        'ProjectedEFC_AtRetention','ProjectedHighSoC_h','UnmodelledMission_h','FleetPass'});
    if isempty(metrics), metrics=z; else, metrics=[metrics;z]; end %#ok<AGROW>
    m=numel(count);
    cy=table(repmat(names(v),m,1),repmat(label,m,1),repmat(c.batteryEOL,m,1), ...
        count,100*depth,100*meanSoC,r.t(rf(:,4)),r.t(rf(:,5)),count*sc, ...
        'VariableNames',{'Vehicle','Case','Retention','Count','DoD_pct','MeanSoC_pct', ...
        'Start_min','End_min','ProjectedCount'});
    if isempty(cycles), cycles=cy; else, cycles=[cycles;cy]; end %#ok<AGROW>
    edges=p.binEdges_pct(:); lo=edges(1:end-1); hi=edges(2:end); w=zeros(numel(lo),1);
    for j=1:numel(lo)
        inside=100*depth>=lo(j) & 100*depth<hi(j);
        if j==numel(lo), inside=100*depth>=lo(j) & 100*depth<=hi(j)+1e-7; end
        w(j)=sum(count(inside));
    end
    assert(abs(sum(w)-sum(count))<1e-8,'Cycle histogram lost cycles.');
    m=numel(w);
    hiTable=table(repmat(names(v),m,1),repmat(label,m,1),repmat(c.batteryEOL,m,1), ...
        lo,hi,w,w*sc,'VariableNames',{'Vehicle','Case','Retention','DoDLower_pct','DoDUpper_pct','Count','ProjectedCount'});
    if isempty(histogram), histogram=hiTable; else, histogram=[histogram;hiTable]; end %#ok<AGROW>
end
end

function h=hoursAbove(x,dt,threshold)
% Exact crossing fraction for the piecewise-linear energy state.
a=x(1:end-1); b=x(2:end); f=zeros(size(dt));
f(a>=threshold & b>=threshold)=1;
cross=(a<threshold & b>threshold) | (b<threshold & a>threshold);
f(cross)=(max(a(cross),b(cross))-threshold)./abs(a(cross)-b(cross));
h=sum(f.*dt);
end
