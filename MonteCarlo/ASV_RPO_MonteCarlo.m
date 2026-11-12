function out = ASV_RPO_MonteCarlo(cfg)
% ASV_RPO_MONTECARLO  FABER RPO design-space Monte Carlo (client + host).
%
% Usage:
%   cfg = ASV_RPO_MC_config();
%   out = ASV_RPO_MonteCarlo(cfg);
%
% Phase sequence per approach:
%   Homing -> Holding 1 -> Closing -> Holding 2 -> Final approach.
%
% Homing and closing: unconstrained 3D two-impulse HCW. Rbar/Vbar
% designate the final hold axis and terminal continuous approach axis.
% Holding points: ideal stationary HCW equilibrium points at [0;-d;0].
% Final: continuous quintic trajectory along selected Vbar or Rbar axis.
% KOZ checked for homing/closing, not forbidden during authorized final
% approach (which is checked against the prescribed axis corridor).
%
% P95 is the empirical 95th percentile of ACCEPTED uniformly sampled
% DESIGN choices; it is NOT a statistical 95% confidence / reliability
% bound unless physically justified input probability models are supplied.
%
% Code is self-contained: does not require proximity_operations_v2.
% Requires MATLAB functions table, prctile alternatives are provided below.

validateCfg(cfg);
rng(cfg.seed, 'twister');
n = sqrt(cfg.mu_km3s2 / (cfg.Re_km+cfg.h_km)^3);

phaseNames = {'Homing','Holding 1','Closing','Holding 2','Final approach'};
labels = [strcat('Client - ',phaseNames),strcat('Host - ',phaseNames)];
maxAttempts = max(cfg.N,ceil(cfg.N*cfg.maxAttemptsFactor));

DV_phase = nan(cfg.N,10);
TOF_phase = nan(cfg.N,10);

% Seven stored values per direction:
% HP1, HP2, Homing TOF, Closing TOF, Final TOF,
% Number of homing transfers, Number of closing transfers
settings = nan(cfg.N,14);

dv_client = nan(cfg.N,1); dv_host = nan(cfg.N,1);
methodClient = strings(cfg.N,1); methodHost = strings(cfg.N,1);
tof_client = nan(cfg.N,1); tof_host = nan(cfg.N,1);
peakSpeed = nan(cfg.N,1);


fail = struct( ...
    'geometry',0, ...
    'singular',0, ...
    'KOZ',0, ...
    'excursion',0, ...
    'speed_homing',0, ...
    'speed_closing',0, ...
    'speed_final',0, ...
    'acceleration',0, ...
    'dv_limit',0, ...
    'numeric',0);


accepted = 0; attempts = 0;

while accepted < cfg.N && attempts < maxAttempts
    attempts = attempts+1;
    [c,reasonC] = optimizeLeg(cfg,cfg.x0_client,n);
    if ~isempty(reasonC)
        fail.(reasonC)=fail.(reasonC)+1;
        continue;
    end
    [h,reasonH] = optimizeLeg(cfg,cfg.x0_host,n);
    if ~isempty(reasonH)
        fail.(reasonH)=fail.(reasonH)+1;
        continue;
    end
    if c.dvTotal+h.dvTotal > cfg.maxDV_perCycle_mps
        fail.dv_limit=fail.dv_limit+1;
        continue;
    end
    accepted=accepted+1;
    DV_phase(accepted,:)=[c.dvPhase h.dvPhase];
    TOF_phase(accepted,:)=[c.tofPhase h.tofPhase];
    settings(accepted,:)=[c.draw h.draw];
    dv_client(accepted)=c.dvTotal;
    dv_host(accepted)=h.dvTotal;
    tof_client(accepted)=c.tofTotal;
    tof_host(accepted)=h.tofTotal;
    peakSpeed(accepted)=max(c.peakSpeed,h.peakSpeed);
    methodClient(accepted)=string(c.method);
    methodHost(accepted)=string(h.method);
end


if accepted == 0

    fprintf('\n========================================\n');
    fprintf('MONTE CARLO FEASIBILITY DIAGNOSTICS\n');
    fprintf('========================================\n');

    fprintf('Attempted cycles: %d\n',attempts);
    fprintf('Accepted cycles:  %d\n',accepted);

    fprintf('\nRejection counts:\n');
    disp(fail);

    error(['No feasible trials after %d attempts. ', ...
        'Check rejection counts above.'],attempts);

end

if accepted<cfg.N
    warning('Only %d/%d requested samples accepted after %d attempts.',...
        accepted,cfg.N,attempts);
end

DV_phase=DV_phase(1:accepted,:);
TOF_phase=TOF_phase(1:accepted,:);
settings=settings(1:accepted,:);
dv_client=dv_client(1:accepted);dv_host=dv_host(1:accepted);
tof_client=tof_client(1:accepted);tof_host=tof_host(1:accepted);
peakSpeed=peakSpeed(1:accepted);
methodClient=methodClient(1:accepted);methodHost=methodHost(1:accepted);
dv_cycle=dv_client+dv_host;
tof_cycle=tof_client+tof_host;

% Phase-level percentiles are NOT additive to cycle P95.
phaseMeanDV=mean(DV_phase,1)';phaseP95DV=colPercentile(DV_phase,95)';
phaseMeanTOF=mean(TOF_phase,1)';phaseP95TOF=colPercentile(TOF_phase,95)';
phaseTable=table(string(labels(:)),phaseMeanDV,phaseP95DV,...
    phaseMeanTOF,phaseP95TOF,'VariableNames',...
    {'Phase','MeanDV_mps','P95DV_mps','MeanTOF_min','P95TOF_min'});

summaryLabels = ["Client RPO";"Host RPO";"Complete retrieval RPO"];
DVs=[dv_client dv_host dv_cycle]; Ts=[tof_client tof_host tof_cycle];
summaryTable=table(summaryLabels,mean(DVs,1)',colPercentile(DVs,50)',...
    colPercentile(DVs,95)',max(DVs,[],1)',mean(Ts,1)',...
    colPercentile(Ts,50)',colPercentile(Ts,95)',max(Ts,[],1)',...
    'VariableNames',{'Scope','MeanDV_mps','P50DV_mps','P95DV_mps',...
    'MaxDV_mps','MeanTOF_min','P50TOF_min','P95TOF_min','MaxTOF_min'});

trialTable=table((1:accepted)',dv_client,dv_host,dv_cycle,...
    tof_client,tof_host,tof_cycle,peakSpeed,methodClient,methodHost,...
    'VariableNames',{'Trial','ClientDV_mps','HostDV_mps','CycleDV_mps',...
    'ClientTOF_min','HostTOF_min','CycleTOF_min','PeakSpeed_mps',...
    'ClientMethod','HostMethod'});

drawNames = { ...
    'ClientHP1_m', ...
    'ClientHP2_m', ...
    'ClientHoming_min', ...
    'ClientClosing_min', ...
    'ClientFinal_min', ...
    'ClientN_Homing', ...
    'ClientN_Closing', ...
    'HostHP1_m', ...
    'HostHP2_m', ...
    'HostHoming_min', ...
    'HostClosing_min', ...
    'HostFinal_min', ...
    'HostN_Homing', ...
    'HostN_Closing'};

trialTable=[trialTable array2table(settings,'VariableNames',drawNames)];

fprintf('\n===============================================\n');
fprintf('FABER RPO MONTE CARLO: %d/%d accepted (%.1f%%)\n',...
    accepted,attempts,100*accepted/attempts);
fprintf('===============================================\n');
disp(summaryTable);disp(phaseTable);
fprintf('Feasibility rejection counts:\n');disp(fail);
[~,idxMinDV]=min(dv_cycle);[~,idxMinTOF]=min(tof_cycle);
fprintf('Minimum-DV feasible cycle: %.5f m/s, %.2f min (trial %d)\n',...
    dv_cycle(idxMinDV),tof_cycle(idxMinDV),idxMinDV);
fprintf('Minimum-TOF feasible cycle: %.5f m/s, %.2f min (trial %d)\n',...
    dv_cycle(idxMinTOF),tof_cycle(idxMinTOF),idxMinTOF);
fprintf('Objective: %s, candidates per approach: %d\n', ...
    cfg.objective,cfg.candidatesPerLeg);
fprintf('Selected client approaches: Vbar %d, Rbar %d\n', ...
    nnz(methodClient=="Vbar"),nnz(methodClient=="Rbar"));
fprintf('Selected host approaches: Vbar %d, Rbar %d\n', ...
    nnz(methodHost=="Vbar"),nnz(methodHost=="Rbar"));
fprintf('NOTE: P95 describes selected sampled designs, NOT reliability.\n');

out=struct();out.cfg=cfg;out.accepted=accepted;out.attempts=attempts;
out.failCounts=fail;out.summaryTable=summaryTable;out.phaseTable=phaseTable;
out.trialTable=trialTable;out.DV_phase_mps=DV_phase;
out.TOF_phase_min=TOF_phase;out.bestDV_trial=idxMinDV;
out.bestTOF_trial=idxMinTOF;
out.P95_DV_cycle_mps=percentile(dv_cycle,95);
out.P95_TOF_cycle_min=percentile(tof_cycle,95);
out.selectedClientMethod=methodClient;out.selectedHostMethod=methodHost;
out.jointP95Exceedance=mean(dv_cycle<=out.P95_DV_cycle_mps & ...
                                 tof_cycle<=out.P95_TOF_cycle_min);

if cfg.saveCSV || cfg.makePlots
    if ~exist(cfg.outputFolder,'dir'),mkdir(cfg.outputFolder);end
end
if cfg.saveCSV
    writetable(trialTable,fullfile(cfg.outputFolder,'RPO_MC_trials.csv'));
    writetable(summaryTable,fullfile(cfg.outputFolder,'RPO_MC_summary.csv'));
    writetable(phaseTable,fullfile(cfg.outputFolder,'RPO_MC_phases.csv'));
end
if cfg.makePlots
    makePlots(cfg,dv_cycle,tof_cycle,DV_phase,TOF_phase,labels,...
        out.P95_DV_cycle_mps,out.P95_TOF_cycle_min,idxMinDV,idxMinTOF);
end
end

function [best,reason] = optimizeLeg(cfg,x0,n)
% Evaluate random candidate times/geometry for every permitted method.
% Select the feasible minimum of the requested objective. This is a
% random-search approximation, not proof of a global optimum.
best=[];reason='geometry';bestScore=inf;
for c=1:cfg.candidatesPerLeg
    % Same sampled geometries/times are offered to both methods, making
    % the route comparison paired within each candidate group.
    drawVars=sampleDesign(cfg);
    for j=1:numel(cfg.methods)
        method=char(cfg.methods{j});
        [trial,why]=simulateLeg(cfg,x0,n,drawVars,method);
        if ~isempty(why),reason=why;continue;end
        % 'weighted' uses explicit reference scales for comparability.
        switch lower(cfg.objective)
            case 'dv'
                score=trial.dvTotal;
            case 'tof'
                score=trial.tofTotal;
            case 'weighted'
                score=cfg.weightDV*(trial.dvTotal/cfg.refDV_mps) + ...
                    cfg.weightTOF*(trial.tofTotal/cfg.refTOF_min);
            otherwise
                error('Invalid objective.');
        end
        if score<bestScore
            best=trial;bestScore=score;reason='';
        end
    end
end
end



function d = sampleDesign(cfg)

% Holding-point geometry
d.hp1 = draw(cfg.hp1_m);
d.hp2 = draw(cfg.hp2_m);
d.final_m = draw(cfg.final_m);

% Total phase durations [s]
d.Th = 60*draw(cfg.homing_min);
d.Tc = 60*draw(cfg.closing_min);
d.Tfinal = 60*draw(cfg.final_min);

% Holding durations [min]
d.T_hold1 = draw(cfg.hold1_min);
d.T_hold2 = draw(cfg.hold2_min);

% Random number of transfers (inclusive)
d.Nh = randi(cfg.N_homing);
d.Nc = randi(cfg.N_closing);

% Minimum duration per individual transfer [s]
minLeg_s = 120;


% Random intermediate waypoint locations
d.homingFrac = sort([rand(1,d.Nh-1),1]);
d.closingFrac = sort([rand(1,d.Nc-1),1]);

% Distance fractions for individual transfers
homingSegments = diff([0 d.homingFrac]);
closingSegments = diff([0 d.closingFrac]);

% Random time-allocation variation (+/- 20%)
homingVariation = 0.8 + 0.4*rand(1,d.Nh);
closingVariation = 0.8 + 0.4*rand(1,d.Nc);

% Weighted allocation: longer segments receive more time
homingWeights = homingSegments .* homingVariation;
closingWeights = closingSegments .* closingVariation;

% Minimum duration per transfer
minLeg_s = 120;

% Allocate the remaining time proportionally
d.homingLegTimes_s = minLeg_s + ...
    (d.Th-d.Nh*minLeg_s)* ...
    homingWeights/sum(homingWeights);

d.closingLegTimes_s = minLeg_s + ...
    (d.Tc-d.Nc*minLeg_s)* ...
    closingWeights/sum(closingWeights);


end


function Tlegs = randomLegTimes(Ttotal,N,minLeg_s)
% RANDOMLEGTIMES
% Randomly divides total phase duration among N transfers.
% Enforces a minimum duration for every transfer.
%
% Inputs:
%   Ttotal   - total phase duration [s]
%   N        - number of transfers
%   minLeg_s - minimum duration per transfer [s]
%
% Output:
%   Tlegs    - 1xN transfer durations [s]

if Ttotal < N*minLeg_s
    error('Total phase time is too short for selected N.');
end

% Generate random positive weights
weights = -log(max(rand(1,N),realmin));

% Normalize weights
weights = weights/sum(weights);

% Time available after minimum allocation
Tremaining = Ttotal-N*minLeg_s;

% Allocate remaining time randomly
Tlegs = minLeg_s + Tremaining*weights;

end


function [res,fail] = simulateLeg(cfg,x0,n,d,method)
res=[];fail='';
hp1=d.hp1;hp2=d.hp2;Th=d.Th;Tc=d.Tc;Tfinal=d.Tfinal;
T_hold1=d.T_hold1;T_hold2=d.T_hold2;final_m=d.final_m;
if ~(hp1>hp2 && hp2>cfg.KOZ_m && hp2<norm(x0(1:3)) && ...
     final_m>0 && final_m<cfg.KOZ_m && cfg.AZ_m>cfg.KOZ_m && ...
     cfg.AZ_m<hp1 && cfg.AZ_m>hp2 && cfg.corridorHalfWidth_m>0)
    fail='geometry';return;
end


%% RANDOMIZED MULTIPLE HOMING TRANSFERS

Nh = d.Nh;

x_start = x0;
x_initial = x0;
x_final = [0;-hp1;0;0;0;0];

dh = 0;
spH = 0;
minH = inf;
maxH = 0;

for k = 1:Nh

    % Randomly selected cumulative distance fraction
    alpha = d.homingFrac(k);

    % Intermediate waypoint position
    r_target = x_initial(1:3) + ...
        alpha*(x_final(1:3)-x_initial(1:3));

    % Zero relative velocity at each waypoint
    x_target = [r_target;zeros(3,1)];

    % Randomly allocated transfer duration
    T_leg = d.homingLegTimes_s(k);

    % Two-impulse HCW maneuver
    [dv_k,sp_k,min_k,max_k,bad] = hcwImpulse( ...
        x_start,x_target,T_leg,n,cfg.nTrajectorySamples);

    if bad
        fail = 'singular';
        return;
    end

    dh = dh + dv_k;

    spH = max(spH,sp_k);
    minH = min(minH,min_k);
    maxH = max(maxH,max_k);

    x_start = x_target;

end



% Closing route choice: Vbar terminates at negative along-track HP2;
% Rbar approaches +radial HP2, requiring a Vbar-to-Rbar transition.
if strcmpi(method,'Vbar')
    targetHP2=[0;-hp2;0;0;0;0];
    rStart=[0;-hp2;0];rEnd=[0;-final_m;0];
    dvRadialHold2=0;
elseif strcmpi(method,'Rbar')
    targetHP2=[hp2;0;0;0;0;0];
    rStart=[hp2;0;0];rEnd=[final_m;0;0];
    % Stationary radial offset requires x-axis acceleration 3*n^2*x.
    dvRadialHold2=3*n^2*hp2*(T_hold2*60);
else
    error('Unknown approach method: %s',method);
end


%% RANDOMIZED MULTIPLE CLOSING TRANSFERS

Nc = d.Nc;

x_start = [0;-hp1;0;0;0;0];
x_initial = x_start;
x_final = targetHP2;

dcl = 0;
spC = 0;
minC = inf;
maxC = 0;

for k = 1:Nc

    % Randomly selected cumulative distance fraction
    alpha = d.closingFrac(k);

    % Intermediate waypoint position
    r_target = x_initial(1:3) + ...
        alpha*(x_final(1:3)-x_initial(1:3));

    % Zero relative velocity at each waypoint
    x_target = [r_target;zeros(3,1)];

    % Randomly allocated transfer duration
    T_leg = d.closingLegTimes_s(k);

    % Two-impulse HCW maneuver
    [dv_k,sp_k,min_k,max_k,bad] = hcwImpulse( ...
        x_start,x_target,T_leg,n,cfg.nTrajectorySamples);

    if bad
        fail = 'singular';
        return;
    end

    dcl = dcl + dv_k;

    spC = max(spC,sp_k);
    minC = min(minC,min_k);
    maxC = max(maxC,max_k);

    x_start = x_target;

end

% Continuous quintic with exact HCW feedforward along selected axis.
[dFinal,spF,accF] = finalQuinticAxis(rStart,rEnd,Tfinal,n,...
    cfg.nTrajectorySamples);
peak=max([spH spC spF]);
if ~all(isfinite([dh dcl dFinal peak minH minC maxC accF]))
    fail='numeric';return;
end


%% GEOMETRY / SAFETY CHECKS

% 1. Actual KOZ penetration
% Homing and closing must remain outside the KOZ.
minRange = min(minH,minC);

if minRange < cfg.KOZ_m
    fail = 'KOZ';
    return;
end

% 2. Excessive outward excursion during closing
% Separate from KOZ penetration.
if maxC > cfg.maxClosingExcursionFactor*hp1
    fail = 'excursion';
    return;
end



%% PHASE-SPECIFIC SPEED CONSTRAINTS

if spH > cfg.maxSpeed_homing_mps
    fail = 'speed_homing';
    return;
end

if spC > cfg.maxSpeed_closing_mps
    fail = 'speed_closing';
    return;
end

if spF > cfg.maxSpeed_final_mps
    fail = 'speed_final';
    return;
end

if accF > cfg.maxPeakAccel_mps2
    fail='acceleration';return;
end
% Closing may cross AZ without stopping there; excessive outward excursions
% beyond the first hold distance are screened approximately.
% Final segment follows either axis, with zero orthogonal components;
% final standoff inside KOZ is permitted, not a capture collision check.
res.dvPhase=[dh,T_hold1*cfg.hold1_DVrate_mps_per_min,...
    dcl,T_hold2*cfg.hold2_DVrate_mps_per_min+dvRadialHold2,dFinal];
res.tofPhase=[Th/60,T_hold1,Tc/60,T_hold2,Tfinal/60];
res.dvTotal=sum(res.dvPhase);res.tofTotal=sum(res.tofPhase);
res.peakSpeed=peak;res.method=method;
% Save geometry, times and maneuver counts
res.draw = [ ...
    hp1 hp2 Th/60 Tc/60 Tfinal/60 Nh Nc];
end

function [dv,maxSp,minRange,maxRange,bad] = hcwImpulse(x0,xf,T,n,N)
% Analytical HCW state-transition matrix, two unconstrained 3D burns.
[Prr,Prv,Pvr,Pvv]=hcwPhi(n,T);
bad=rcond(Prv)<1e-11;
if bad,dv=NaN;maxSp=NaN;minRange=NaN;maxRange=NaN;return;end
v1=Prv\(xf(1:3)-Prr*x0(1:3));
b1=v1-x0(4:6);
vPre=Pvr*x0(1:3)+Pvv*v1;
b2=xf(4:6)-vPre;
dv=norm(b1)+norm(b2);
maxSp=0;minRange=inf;maxRange=0;
for k=1:N
    tau=T*(k-1)/(N-1);
    [A,B,C,D]=hcwPhi(n,tau);
    r=A*x0(1:3)+B*v1;
    v=C*x0(1:3)+D*v1;
    maxSp=max(maxSp,norm(v));
    minRange=min(minRange,norm(r));
    maxRange=max(maxRange,norm(r));
end
end

function [A,B,C,D] = hcwPhi(n,t)
nt=n*t;c=cos(nt);s=sin(nt);
A=[4-3*c,0,0;6*(s-nt),1,0;0,0,c];
B=[s,2*(1-c),0;2*(c-1),4*s-3*nt,0;0,0,s]/n;
C=[3*n*s,0,0;6*n*(c-1),0,0;0,0,-n*s];
D=[c,2*s,0;-2*s,4*c-3,0;0,0,c];
end

function [dv,maxSpeed,peakAccel]=finalQuinticAxis(rStart,rEnd,T,n,N)
% Continuous exact-feedforward HCW quintic in LVLH; starts/ends at rest.
t=linspace(0,T,N);s=t/T;
p=10*s.^3-15*s.^4+6*s.^5;
pd=(30*s.^2-60*s.^3+30*s.^4)/T;
pdd=(60*s-180*s.^2+120*s.^3)/T^2;
dr=rEnd-rStart;
r=rStart+dr*p;
v=dr*pd;
a=dr*pdd;
u=[a(1,:)-3*n*n*r(1,:)-2*n*v(2,:);
   a(2,:)+2*n*v(1,:);
   a(3,:)+n*n*r(3,:)];
aMag=sqrt(sum(u.^2,1));
dv=trapz(t,aMag);
maxSpeed=max(sqrt(sum(v.^2,1)));peakAccel=max(aMag);
end

function x=draw(bounds)
x=bounds(1)+(bounds(2)-bounds(1))*rand();
end

function p=percentile(x,q)
% MATLAB prctile-style linear interpolation without Statistics Toolbox.
x=sort(x(:));n=numel(x);
if n==1,p=x;return;end
r=1+(n-1)*q/100;lo=floor(r);hi=ceil(r);
p=x(lo)+(r-lo)*(x(hi)-x(lo));
end

function vals=colPercentile(X,q)
vals=zeros(1,size(X,2));
for j=1:size(X,2),vals(j)=percentile(X(:,j),q);end
end

function validateCfg(cfg)

% Validate number of transfers
assert(all(cfg.N_homing >= 1) && ...
    all(cfg.N_homing == floor(cfg.N_homing)) && ...
    cfg.N_homing(2) >= cfg.N_homing(1), ...
    'Invalid homing maneuver-count bounds.');

assert(all(cfg.N_closing >= 1) && ...
    all(cfg.N_closing == floor(cfg.N_closing)) && ...
    cfg.N_closing(2) >= cfg.N_closing(1), ...
    'Invalid closing maneuver-count bounds.');

assert(isstruct(cfg),'Input must be a configuration struct.');
assert(cfg.N>=1 && cfg.N==floor(cfg.N),'N must be a positive integer.');
assert(cfg.nTrajectorySamples>=3,'Need >=3 trajectory samples.');
assert(cfg.maxAttemptsFactor>=1);
assert(cfg.h_km>100 && cfg.Re_km>0 && cfg.mu_km3s2>0);
assert(all(size(cfg.x0_client)==[6 1]) && all(size(cfg.x0_host)==[6 1]));
keys={'hp1_m','hp2_m','final_m','homing_min','closing_min',...
 'final_min','hold1_min','hold2_min'};
for k=1:numel(keys)
 b=cfg.(keys{k});
 assert(isnumeric(b) && numel(b)==2 && all(isfinite(b)) && ...
        b(1)>0 && b(2)>=b(1),...
        'Invalid bounds for %s (need positive [min max]).',keys{k});
end
assert(cfg.KOZ_m>0 && cfg.AZ_m>cfg.KOZ_m);
assert(cfg.maxSpeed_mps>0 && cfg.maxPeakAccel_mps2>0);
assert(cfg.candidatesPerLeg>=1 && cfg.candidatesPerLeg==floor(cfg.candidatesPerLeg));
assert(any(strcmpi(cfg.objective,{'dv','tof','weighted'})));
assert(iscell(cfg.methods) && ~isempty(cfg.methods));
assert(all(ismember(lower(string(cfg.methods)),["vbar","rbar"])));
assert(cfg.refDV_mps>0 && cfg.refTOF_min>0);
assert(cfg.weightDV>=0 && cfg.weightTOF>=0 && cfg.weightDV+cfg.weightTOF>0);
end

function makePlots(cfg,dv,tof,DVphase,TOFphase,labels,p95dv,p95tof,iDV,iTOF)
f1=figure('Color','w','Name','RPO Monte Carlo - Delta V');
histogram(dv,40);hold on;grid on;xline(p95dv,'--','P95','LineWidth',1.6);
xlabel('Client + Host RPO Delta V [m/s]');ylabel('Accepted designs');
title(sprintf('RPO Delta V (N=%d, P95=%.3f m/s)',numel(dv),p95dv));
saveas(f1,fullfile(cfg.outputFolder,'RPO_MC_DV_histogram.png'));

f2=figure('Color','w','Name','RPO Monte Carlo - TOF');
histogram(tof/60,40);hold on;grid on;xline(p95tof/60,'--','P95','LineWidth',1.6);
xlabel('Client + Host RPO TOF [hours]');ylabel('Accepted designs');
title(sprintf('RPO TOF (P95=%.2f h)',p95tof/60));
saveas(f2,fullfile(cfg.outputFolder,'RPO_MC_TOF_histogram.png'));

f3=figure('Color','w','Name','RPO Monte Carlo - TOF vs Delta V');
scatter(tof/60,dv,11,'filled');hold on;grid on;
plot(tof(iDV)/60,dv(iDV),'kp','MarkerSize',13,'MarkerFaceColor','k');
plot(tof(iTOF)/60,dv(iTOF),'ks','MarkerSize',9,'LineWidth',2);
xline(p95tof/60,'--','P95 TOF');yline(p95dv,'--','P95 DV');
xlabel('Complete retrieval RPO TOF [hours]');ylabel('Complete retrieval RPO Delta V [m/s]');
legend('Feasible trials','Min DV','Min TOF','Location','best');
title('TOF-DV trade space (design-sampling, not probability of failure)');
saveas(f3,fullfile(cfg.outputFolder,'RPO_MC_TOF_vs_DV.png'));

f4=figure('Color','w','Name','RPO Monte Carlo - Phase P95 DV');
barh(colPercentile(DVphase,95));grid on;
set(gca,'YTick',1:numel(labels),'YTickLabel',labels,'YDir','reverse');
xlabel('Per-phase P95 Delta V [m/s]');title('Individual phase P95 (do not sum)');
saveas(f4,fullfile(cfg.outputFolder,'RPO_MC_Phase_P95_DV.png'));

f5=figure('Color','w','Name','RPO Monte Carlo - Phase P95 TOF');
barh(colPercentile(TOFphase,95));grid on;
set(gca,'YTick',1:numel(labels),'YTickLabel',labels,'YDir','reverse');
xlabel('Per-phase P95 TOF [min]');title('Individual phase P95 (do not sum)');
saveas(f5,fullfile(cfg.outputFolder,'RPO_MC_Phase_P95_TOF.png'));
end
