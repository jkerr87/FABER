function out = FABER_RPO_compare_methods(opts)
% FABER_RPO_COMPARE_METHODS
% Standalone RPO-only, paired Monte Carlo comparison:
%   'Vbar'       analytical_manos (one orbital period, tangential burns)
%   'Rbar'       analytical_manos (half period, radial burns)
%   'GeneralHCW' fixed-time two-impulse HCW transfer
%
% H1: general HCW for all (initial acquisition state is not a stationary
%     V-bar point, so analytical_manos cannot be applied there).
% H2 and Closing: method-specific transfer between stationary V-bar HPs.
% Final: rpo_continuous(...,'Vbar',...) for ALL methods.
% This is a method comparison, NOT an optimizer or a certified safe CONOPS.
%
% Needs on MATLAB path:
% ASV_RPO_MC_config.m, analytical_manos.m, rpo_continuous.m
% Optionally ASV_RPO_MonteCarlo.m for the unchanged legacy baseline.
%
% EXAMPLE:
%   opts=struct('N',60,'runLegacyBaseline',true,'makePlots',true);
%   R=FABER_RPO_compare_methods(opts);
%
% Important: the baseline cfg.controlled_min=[10 45] min EXCLUDES both
% analytical transfers at ~680 km (Rbar ~49 min, Vbar ~98 min). By default
% this comparison EXTENDS the HP2->HP3 window for ALL methods to permit the
% comparisons; the original time requirement is recorded separately.

if nargin<1 || isempty(opts), opts=struct(); end
if ~isfield(opts,'N'), opts.N=60; end
if ~isfield(opts,'seed'), opts.seed=42; end
if ~isfield(opts,'h_km'), opts.h_km=680.12; end
if ~isfield(opts,'expandClosingWindow'), opts.expandClosingWindow=true; end
if ~isfield(opts,'runLegacyBaseline'), opts.runLegacyBaseline=true; end
if ~isfield(opts,'makePlots'), opts.makePlots=true; end
if ~isfield(opts,'saveCSV'), opts.saveCSV=true; end
if ~isfield(opts,'folder'), opts.folder='FABER_RPO_method_comparison'; end
if ~isfield(opts,'nSamples'), opts.nSamples=51; end
if ~isfield(opts,'homing2SpeedOverride_mps'), opts.homing2SpeedOverride_mps=[]; end
assert(opts.N>=1 && opts.N==floor(opts.N));

cfg=ASV_RPO_MC_config();
cfg.h_km=opts.h_km;
cfg.nTrajectorySamples=max(21,opts.nSamples);
cfg.makePlots=false; cfg.saveCSV=false;
originalControlled=cfg.controlled_min;
if ~isempty(opts.homing2SpeedOverride_mps)
    assert(opts.homing2SpeedOverride_mps>0);
    cfg.maxSpeed_closing_mps=opts.homing2SpeedOverride_mps;
    fprintf('WARNING: exploratory Homing 2 speed override = %.3f m/s (baseline = 0.3).\n', ...
        opts.homing2SpeedOverride_mps);
end

n=sqrt(cfg.mu_km3s2/(cfg.Re_km+cfg.h_km)^3);
P=2*pi/n;
TAnalytical_min=[P/60,P/120]; % Vbar, Rbar
if opts.expandClosingWindow
    cfg.controlled_min(2)=max(cfg.controlled_min(2),ceil(P/60)+1);
end

fprintf('\nReference orbit: %.2f km | period %.2f min\n',cfg.h_km,P/60);
fprintf('Analytical H2/closing durations: Vbar %.2f min, Rbar %.2f min\n',TAnalytical_min);
fprintf('Original controlled closing window: [%.1f, %.1f] min\n',originalControlled);
fprintf('Study controlled closing window: [%.1f, %.1f] min\n',cfg.controlled_min);
if ~opts.expandClosingWindow
    fprintf('NOTE: no full Vbar/Rbar trials will pass the original 45-min closing limit.\n');
else
    fprintf('CAUTION: study time window is relaxed versus baseline for every method.\n');
end

fprintf('Corridor/excursion checks remain DIAGNOSTIC ONLY, as in your MC baseline.\n');

names={'Vbar','Rbar','GeneralHCW'};
N=opts.N;
records=cell(N,3);
DV=nan(N,3,2); TOF=nan(N,3,2);
reason=strings(N,3);
strictTime=false(N,3);
corridorOK=false(N,3);

rng(opts.seed,'twister');
for i=1:N
    % SAME random geometry, initial states and general-HCW durations for
    % all methods in this attempt (paired candidate comparison).
    d=sampleOneDesign(cfg,n);
    xC=cfg.x0_client+(2*rand(6,1)-1).*cfg.dx0_client;
    xH=cfg.x0_host+(2*rand(6,1)-1).*cfg.dx0_host;
    for m=1:3
        [a,whyA]=simulateOne(cfg,n,xC,d,names{m},originalControlled);
        if ~isempty(whyA), reason(i,m)=string(['Client: ' whyA]); continue; end
        [b,whyB]=simulateOne(cfg,n,xH,d,names{m},originalControlled);
        if ~isempty(whyB), reason(i,m)=string(['Host: ' whyB]); continue; end
        if a.DV+b.DV>cfg.maxDV_perCycle_mps
            reason(i,m)="Cycle DV cap"; continue;
        end

        records{i,m}=struct('client',a,'host',b,'draw',d);
        DV(i,m,1)=a.DV;DV(i,m,2)=b.DV;
        TOF(i,m,1)=a.TOF_min;TOF(i,m,2)=b.TOF_min;
        strictTime(i,m)=a.strictTime && b.strictTime;
        corridorOK(i,m)=a.corridorOK && b.corridorOK;
    end
end

cycleDV=sum(DV,3); cycleTOF=sum(TOF,3);
valid=isfinite(cycleDV);
common=all(valid,2);
result=cell(3,1);
for m=1:3
    v=cycleDV(valid(:,m),m);t=cycleTOF(valid(:,m),m);
    nOk=numel(v);
    if nOk>0
        [sorted,ix]=sort(v);
        ii=find(valid(:,m));selected=ii(ix(max(1,ceil(.95*nOk))));
        medDV=median(v); p95DV=sorted(max(1,ceil(.95*nOk)));
        medTOF=median(t);p95TOF=empP95(t);
    else
        selected=NaN; medDV=NaN;p95DV=NaN;medTOF=NaN;p95TOF=NaN;
    end

    result{m}=struct('method',names{m},'selectedTrial',selected, ...
        'accepted',nOk,'P95DV_mps',p95DV,'P95TOF_min',p95TOF);
    fprintf('%-11s %3d / %3d accepted | P95 DV %.4f m/s | P95 TOF %.1f min | original times satisfied %d\n', ...
        names{m},nOk,N,p95DV,p95TOF,sum(strictTime(:,m)));
    summaryFields(m,:)={string(names{m}),nOk,N,nOk/N,medDV,p95DV,medTOF,p95TOF, ...
        sum(strictTime(:,m)),sum(corridorOK(:,m))}; %#ok<AGROW>
end

summary=cell2table(summaryFields, ...
    'VariableNames',{'Method','Accepted','Attempts','AcceptedFraction', ...
    'MedianCycleDV_mps','P95CycleDV_mps','MedianCycleTOF_min', ...
    'P95CycleTOF_min','CountWithinOriginalClosingTime','CountCorridorDiagnosticPass'});
% Convert table numeric columns (mixed cell->table produces cell variables).
for k = 2:width(summary)
    vn = summary.Properties.VariableNames{k};

    if iscell(summary.(vn))
        summary.(vn) = cell2mat(summary.(vn));
    end
end

fprintf('\nCommon feasible candidates for all 3 methods: %d / %d\n',sum(common),N);
paired=table(string(names(:)),nan(3,1),nan(3,1), ...
    'VariableNames',{'Method','MedianCommonDV_mps','P95CommonDV_mps'});
for m=1:3
    v=cycleDV(common,m);
    if ~isempty(v), paired.MedianCommonDV_mps(m)=median(v);paired.P95CommonDV_mps(m)=empP95(v);end
end

disp(summary);disp(paired);

legacy=[];
if opts.runLegacyBaseline
    fprintf('\nRunning UNCHANGED legacy ASV_RPO_MonteCarlo for context ...\n');
    bcfg=ASV_RPO_MC_config();bcfg.h_km=opts.h_km;
    bcfg.N=opts.N;bcfg.seed=opts.seed;
    bcfg.makePlots=false;bcfg.saveCSV=false;
    legacy=ASV_RPO_MonteCarlo(bcfg);
    fprintf('Legacy P95 COMPLETE cycle DV: %.4f m/s (different sampling/segmentation).\n', ...
        legacy.P95_DV_cycle_mps);
end

out=struct('cfg',cfg,'originalControlled_min',originalControlled, ...
    'methods',{names},'summary',summary,'pairedSummary',paired, ...
    'records',{records},'DV_mps',DV,'TOF_min',TOF, ...
    'valid',valid,'commonFeasible',common,'strictTime',strictTime, ...
    'corridorDiagnosticPass',corridorOK,'rejectionReason',reason, ...
    'selected',{result},'legacy',legacy);

if opts.saveCSV || opts.makePlots
    if ~exist(opts.folder,'dir'),mkdir(opts.folder);end
end

if opts.saveCSV
    writetable(summary,fullfile(opts.folder,'RPO_methods_summary.csv'));
    writetable(paired,fullfile(opts.folder,'RPO_methods_paired_summary.csv'));
end

if opts.makePlots
    plotMethodApproaches(out,opts);
end
end

function d=sampleOneDesign(c,n)
d.hp1=draw(c.hp1_m);d.hp2=draw(c.hp2_m);
d.hp3=draw(c.hp3_m);d.fin=draw(c.final_m);
d.Th=draw(c.homing_min)*60;
d.Tc=draw(c.closing_min)*60;
d.Tcontrolled=draw(c.controlled_min)*60;
d.h1=draw(c.hold1_min);d.h2=draw(c.hold2_min);
d.h3=draw(c.hold3_min);d.Tfinal=draw(c.final_min)*60;
d.period_s=2*pi/n;
end

function [r,why]=simulateOne(c,n,x0,d,method,originalControlled)
r=[];why='';
if ~(norm(x0(1:3))>d.hp1 && d.hp1>d.hp2 && d.hp2>c.AZ_m && ...
        c.AZ_m>c.KOZ_m && d.hp3<c.KOZ_m && d.hp3>d.fin && d.fin>0)
    why='geometry';return;
end

hp1=[0;-d.hp1;0;0;0;0]; hp2=[0;-d.hp2;0;0;0;0];
hp3=[0;-d.hp3;0;0;0;0]; endX=[0;-d.fin;0;0;0;0];
trace=cell(1,7); phaseDV=zeros(1,7);phaseT=zeros(1,7);

% Homing 1 is general HCW for ALL methods, because x0 is not on V-bar.
[s1,e]=generalHCW(x0,hp1,d.Th,n,c.nTrajectorySamples);
if ~isempty(e),why=['H1 ' e];return;end
if s1.maxSpeed>c.maxSpeed_homing_mps || s1.minRange<c.KOZ_m
    why='H1 speed/KOZ';return;
end

trace{1}=s1.trajectory(1:3,:);phaseDV(1)=s1.DV;phaseT(1)=d.Th/60;
trace{2}=repmat(hp1(1:3),1,2);phaseDV(2)=d.h1*c.hold1_DVrate_mps_per_min;phaseT(2)=d.h1;

% Homing 2: matched stationary V-bar endpoints.
[s2,e]=methodTransfer(hp1,hp2,d.Tc,method,n,c.nTrajectorySamples);
if ~isempty(e),why=['H2 ' e];return;end
if s2.maxSpeed>c.maxSpeed_closing_mps || s2.minRange<c.KOZ_m || ...
        s2.TOF_s/60<c.closing_min(1) || s2.TOF_s/60>c.closing_min(2)
    why='H2 speed/time/KOZ';return;
end

trace{3}=s2.trajectory(1:3,:);phaseDV(3)=s2.DV;phaseT(3)=s2.TOF_s/60;
trace{4}=repmat(hp2(1:3),1,2);phaseDV(4)=d.h2*c.hold2_DVrate_mps_per_min;phaseT(4)=d.h2;

% Closing: analytical fixed-time Vbar/Rbar, or sampled-time general HCW.
[s3,e]=methodTransfer(hp2,hp3,d.Tcontrolled,method,n,c.nTrajectorySamples);
if ~isempty(e),why=['Closing ' e];return;end
if s3.maxSpeed>c.maxSpeed_controlled_mps || s3.minRange<d.fin || ...
        s3.TOF_s/60<c.controlled_min(1) || s3.TOF_s/60>c.controlled_min(2)
    why='Closing speed/time/range';return;
end

trace{5}=s3.trajectory(1:3,:);phaseDV(5)=s3.DV;phaseT(5)=s3.TOF_s/60;
% Baseline's corridor and inward-range tests are commented out. Keep them
% diagnostic for fair comparison; NOT a safety certification.
xyz=s3.trajectory(1:3,:); rr=vecnorm(xyz,2,1);
inside=rr<c.KOZ_m;
corridorOK=~any(inside) || all(hypot(xyz(1,inside),xyz(3,inside))<=c.corridorHalfWidth_m ...
    & xyz(2,inside)<0);
if any(inside)
    ix=find(inside,1,'first');
    corridorOK=corridorOK && all(diff(rr(ix:end))<=1e-3);
end

trace{6}=repmat(hp3(1:3),1,2);phaseDV(6)=d.h3*c.hold3_DVrate_mps_per_min;phaseT(6)=d.h3;
% EXACT uploaded final-control routine, with quintic endpoint conditions.
f=rpo_continuous(hp3,endX,n,d.Tfinal,'Vbar',max(101,c.nTrajectorySamples));
if ~f.feasible || f.maxSpeed_mps>c.maxSpeed_final_mps || ...
        f.peakAcceleration_mps2>c.maxPeakAccel_mps2
    why='Final speed/acceleration';return;
end

trace{7}=f.trajectory(1:3,:);phaseDV(7)=f.DV_mps;phaseT(7)=f.TOF_min;

r=struct('DV',sum(phaseDV),'TOF_min',sum(phaseT), ...
    'phaseDV',phaseDV,'phaseTOF',phaseT,'trace',{trace}, ...
    'corridorOK',corridorOK, ...
    'strictTime',s3.TOF_s/60>=originalControlled(1) ...
        && s3.TOF_s/60<=originalControlled(2));
end

function [s,why]=methodTransfer(x0,xf,T,method,n,N)
s=[];why='';
if strcmpi(method,'GeneralHCW')
    [s,why]=generalHCW(x0,xf,T,n,N);
else
    % Current uploaded analytical_manos.m: stop->stop on V-bar ONLY.
    a=analytical_manos(x0(1:3),xf(1:3),n,method,N,false);
    if ~a.feasible,why='analytical endpoint';return;end
    s=struct('DV',a.DV_mps,'TOF_s',a.TOF_s, ...
        'trajectory',a.trajectory,'minRange',a.minRange_m, ...
        'maxSpeed',a.maxSpeed_mps);
end
end

function [s,why]=generalHCW(x0,xf,T,n,N)
s=[];why='';
if T<=0,why='time';return;end
P=hcwPhi(n,T);prv=P(1:3,4:6);
if rcond(prv)<1e-10,why='singular Phi_rv';return;end
v=prv\(xf(1:3)-P(1:3,1:3)*x0(1:3));
dv1=v-x0(4:6);
xp=[x0(1:3);v];
t=linspace(0,T,N);traj=hcwTraj(n,t,xp);
dv2=xf(4:6)-traj(4:6,end);
if norm(traj(1:3,end)-xf(1:3))>1e-5
    why='terminal position';return;
end

s=struct('DV',norm(dv1)+norm(dv2),'TOF_s',T, ...
    'trajectory',traj,'minRange',min(vecnorm(traj(1:3,:),2,1)), ...
    'maxSpeed',max(vecnorm(traj(4:6,:),2,1)));
end

function P=hcwPhi(n,t)
c=cos(n*t);s=sin(n*t);nt=n*t;
P=[4-3*c,0,0,s/n,2*(1-c)/n,0; ...
   6*(s-nt),1,0,-2*(1-c)/n,(4*s-3*nt)/n,0; ...
   0,0,c,0,0,s/n; ...
   3*n*s,0,0,c,2*s,0; ...
   6*n*(c-1),0,0,-2*s,4*c-3,0; ...
   0,0,-n*s,0,0,c];
end

function traj=hcwTraj(n,t,X)
nt=n*t;c=cos(nt);s=sin(nt);
x=X(1);y=X(2);z=X(3);vx=X(4);vy=X(5);vz=X(6);
traj=[(4-3*c)*x+s/n*vx+2*(1-c)/n*vy; ...
    6*(s-nt)*x+y-2*(1-c)/n*vx+(4*s-3*nt)/n*vy; ...
    c*z+s/n*vz; ...
    3*n*s*x+c*vx+2*s*vy; ...
    6*n*(c-1)*x-2*s*vx+(4*c-3)*vy; ...
    -n*s*z+c*vz];
end

function v=draw(b)
v=b(1)+(b(2)-b(1))*rand;
end

function p=empP95(v)
v=sort(v(:));if isempty(v),p=NaN;else,p=v(max(1,ceil(.95*numel(v))));end
end

function plotMethodApproaches(out,opts)
% Exactly two figures: one to client and one to ATLANTE.
colors=lines(3); [sx,sy,sz]=sphere(20);
for leg=1:2
    f=figure('Color','w','Position',[110 110 1120 800]);ax=axes('Parent',f);
    hold(ax,'on');grid(ax,'on');axis(ax,'equal');
    for m=1:3
        sel=out.selected{m}.selectedTrial;
        if isnan(sel),continue;end
        rec=out.records{sel,m};
        if leg==1,p=rec.client;else,p=rec.host;end
        for k=1:7
            xyz=p.trace{k};
            if k==1,visibility='on';else,visibility='off';end
            plot3(ax,xyz(1,:),xyz(2,:),xyz(3,:), ...
                'Color',colors(m,:),'LineWidth',1.8, ...
                'DisplayName',out.methods{m},'HandleVisibility',visibility);
        end
    end

    A=out.cfg.AZ_m; K=out.cfg.KOZ_m;
    surf(ax,A*sx,A*sy,A*sz,'FaceColor','none','EdgeColor',[.4 .6 1], ...
        'EdgeAlpha',.10,'DisplayName','AZ');
    surf(ax,K*sx,K*sy,K*sz,'FaceColor','none','EdgeColor',[1 .3 .3], ...
        'EdgeAlpha',.25,'DisplayName','KOZ');
    xlabel(ax,'Radial [m]');ylabel(ax,'Along-track [m]');zlabel(ax,'Cross-track [m]');
    if leg==1,legName='Client';else,legName='ATLANTE';end
    title(ax,['Selected method-specific P95 RPO: Taxi to ' legName]);
    view(ax,35,25);legend(ax,'Location','eastoutside');
    if opts.saveCSV
        exportgraphics(f,fullfile(opts.folder,['RPO_method_comparison_' legName '.png']), ...
            'Resolution',200);
    end
end
end
