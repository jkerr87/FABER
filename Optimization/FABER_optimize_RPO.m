function R = FABER_optimize_RPO(userCfg)
% FABER_OPTIMIZE_RPO  Joint homing+closing nominal HCW optimization.
% One client OR host approach per run. Homing and closing each have
% 1..Nmax HCW coast arcs through V-bar fly-through waypoints.
% Intermediate arc boundaries use velocity-correction impulses (no stop).
% Only HP1 and HP2 require zero terminal velocities.
% No orbit-period timing assumption. Holds and quintic V-bar final approach
% are included in DV and TOF; their durations are fixed for this starter.
% Requires Optimization Toolbox (fmincon).
% Run: R=FABER_optimize_RPO();
% cfg=struct('x0',[50;-2000;20;0.01;0.04;-0.01]); R=FABER_optimize_RPO(cfg);

cfg=defaults();
if nargin>0 && ~isempty(userCfg)
    f=fieldnames(userCfg);
    for j=1:numel(f), cfg.(f{j})=userCfg.(f{j}); end
end
assert(exist('fmincon','file')==2,'fmincon requires Optimization Toolbox.');
assert(numel(cfg.x0)==6,'x0 must have 6 elements.');
assert(all(diff(cfg.tofCaps_min)>0),'TOF caps must be increasing.');
assert(cfg.Nmax>=1 && cfg.Nmax<=4 && cfg.minLeg_min>0);
assert(cfg.minWaypointGap>0 && cfg.Nmax*cfg.minWaypointGap<1);
assert(cfg.hp1_m(1)>cfg.hp2_m(2) && cfg.hp2_m(1)>cfg.KOZ_m);
assert(cfg.finalStandoff_m>0 && cfg.finalStandoff_m<cfg.KOZ_m);
assert(cfg.validationSamples>=5 && cfg.nSamples>=5);
assert(cfg.maxSpeed_homing_mps>0 && cfg.maxSpeed_closing_mps>0 && cfg.maxSpeed_final_mps>0);
n=sqrt(cfg.mu_km3s2/(cfg.Re_km+cfg.h_km)^3);
rng(cfg.seed,'twister');
opts=optimoptions('fmincon','Algorithm','sqp','Display','off', ...
    'MaxIterations',cfg.maxIterations,'MaxFunctionEvaluations',cfg.maxFunctionEvaluations, ...
    'ConstraintTolerance',1e-6,'StepTolerance',1e-9);

% z = [HP1;HP2; homing internal fractions; closing internal fractions;
%      homing arc durations (min); closing arc durations (min)]
rows=struct('Nh',{},'Nc',{},'TOFcap_min',{},'DV_mps',{},'TOF_min',{}, ...
    'HomingDV_mps',{},'ClosingDV_mps',{},'FinalDV_mps',{}, ...
    'HP1_m',{},'HP2_m',{},'PeakHoming_mps',{},'PeakClosing_mps',{}, ...
    'MinRange_m',{},'Verified',{},'z',{});
for Nh=1:cfg.Nmax
    for Nc=1:cfg.Nmax
        for ic=1:numel(cfg.tofCaps_min)
            cap=cfg.tofCaps_min(ic);
            fixedTime=cfg.hold1_min+cfg.hold2_min+cfg.final_min;
            if cap<fixedTime+(Nh+Nc)*cfg.minLeg_min,continue;end
            [lb,ub,A,b,idx]=makeBounds(Nh,Nc,cap,cfg);
            best=[];
            for s=1:cfg.nStarts
                z0=initialGuess(Nh,Nc,cap,cfg,idx);
                try
                    [z,~,flag]=fmincon(@(zz)objective(zz,Nh,Nc,n,cfg,idx), ...
                        z0,A,b,[],[],lb,ub, ...
                        @(zz)nonlinear(zz,Nh,Nc,n,cfg,idx),opts);
                catch ME
                    warning('Nh=%d Nc=%d cap=%g: %s',Nh,Nc,cap,ME.message);
                    continue;
                end
                test=evaluate(z,Nh,Nc,n,cfg,idx,cfg.validationSamples,true);
                if flag<=0 || ~isFeasible(test,cap,cfg),continue;end
                if isempty(best) || test.dv<best.test.dv
                    best=struct('z',z,'test',test);
                end
            end
            if ~isempty(best)
                t=best.test;
                rows(end+1)=struct('Nh',Nh,'Nc',Nc,'TOFcap_min',cap, ...
                    'DV_mps',t.dv,'TOF_min',t.tof,'HomingDV_mps',t.hdv, ...
                    'ClosingDV_mps',t.cdv,'FinalDV_mps',t.fdv, ...
                    'HP1_m',best.z(1),'HP2_m',best.z(2), ...
                    'PeakHoming_mps',t.spH,'PeakClosing_mps',t.spC, ...
                    'MinRange_m',min(t.minH,t.minC), ...
                    'Verified',true,'z',best.z); %#ok<AGROW>
            else
                fprintf('No verified design Nh=%d Nc=%d cap=%g min\n',Nh,Nc,cap);
            end
        end
    end
end
R=struct('cfg',cfg,'solutions',table(),'pareto',table(), ...
    'raw',rows,'bestHomingArcs',table(),'bestClosingArcs',table());
if isempty(rows)
    warning('No verified designs found; revisit constraints and optimizer starts.');
    return;
end
S=struct2table(rmfield(rows,'z'));
mask=paretoMask(S.TOF_min,S.DV_mps);
P=sortrows(S(mask,:),'TOF_min');
R.solutions=S;R.pareto=P;
fprintf('\nJOINT RPO VERIFIED DESIGNS\n');disp(S);
fprintf('\nJOINT RPO PARETO FRONTIER\n');disp(P);
[~,ib]=min([rows.DV_mps]);
[H,C]=arcTables(rows(ib).z,rows(ib).Nh,rows(ib).Nc,n,cfg);
R.bestHomingArcs=H;R.bestClosingArcs=C;
fprintf('\nLOWEST-DV HOMING ARCS\n');disp(H);
fprintf('\nLOWEST-DV CLOSING ARCS\n');disp(C);
if cfg.makePlots
    figure('Color','w','Name','FABER joint RPO Pareto');
    scatter(S.TOF_min,S.DV_mps,38,S.Nh+S.Nc,'filled');hold on;grid on;
    plot(P.TOF_min,P.DV_mps,'k-o','LineWidth',1.7);
    xlabel('Complete approach TOF [min]');ylabel('Complete approach Delta V [m/s]');
    title('Joint homing + closing HCW optimization (single approach)');
    cb=colorbar;ylabel(cb,'Homing + closing arc count');
    legend('Verified candidates','Sampled Pareto','Location','best');
end
if cfg.saveCSV
    writetable(S,'FABER_joint_RPO_solutions.csv');
    writetable(P,'FABER_joint_RPO_pareto.csv');
    writetable(H,'FABER_joint_RPO_best_homing.csv');
    writetable(C,'FABER_joint_RPO_best_closing.csv');
end
end

function cfg=defaults()
cfg.mu_km3s2=398600.4418;cfg.Re_km=6371;cfg.h_km=457.25;
cfg.x0=[50;-2000;20;0.01;0.04;-0.01];
cfg.hp1_m=[800 1200];cfg.hp2_m=[100 250];
cfg.finalStandoff_m=5;
cfg.Nmax=4;cfg.tofCaps_min=[170 210 260];
cfg.homing_min=[30 100];cfg.closing_min=[20 130];
cfg.hold1_min=30;cfg.hold2_min=10;cfg.final_min=10;
cfg.minLeg_min=5;cfg.minWaypointGap=0.05;
cfg.maxSpeed_homing_mps=1.0;cfg.maxSpeed_closing_mps=1.0;
cfg.maxSpeed_final_mps=1.0;cfg.KOZ_m=30;
cfg.homingExcursionFactor=1.20;cfg.closingExcursionFactor=1.20;
cfg.maxPeakAccel_mps2=inf;
cfg.hold1_DVrate_mps_per_min=0;cfg.hold2_DVrate_mps_per_min=0;
cfg.nStarts=3;cfg.seed=42;cfg.nSamples=21;cfg.validationSamples=201;
cfg.maxIterations=130;cfg.maxFunctionEvaluations=2300;
cfg.terminalTol_m=1e-2;cfg.terminalVelTol_mps=1e-5;
cfg.saveCSV=true;cfg.makePlots=true;
end

function [lb,ub,A,b,idx]=makeBounds(Nh,Nc,cap,cfg)
idx.h=(3:2+Nh-1); % homing internal fractions
idx.c=(3+Nh-1:2+Nh+Nc-2); % closing internal fractions
% More robust explicit cursor, including empty fraction sets
p=3;idx.h=p:(p+Nh-2);p=p+Nh-1;
idx.c=p:(p+Nc-2);p=p+Nc-1;
idx.th=p:(p+Nh-1);p=p+Nh;
idx.tc=p:(p+Nc-1);
ng=cfg.minWaypointGap;
lb=[cfg.hp1_m(1);cfg.hp2_m(1);repmat(ng,Nh+Nc-2,1); ...
    repmat(cfg.minLeg_min,Nh+Nc,1)];
ub=[cfg.hp1_m(2);cfg.hp2_m(2);repmat(1-ng,Nh+Nc-2,1); ...
    repmat(max(cfg.homing_min(2),cfg.closing_min(2)),Nh+Nc,1)];
K=numel(lb);A=zeros(0,K);b=zeros(0,1);
% Enforce fractions in increasing order with minimum separation
for inds={idx.h,idx.c}
    ii=inds{1};
    for k=2:numel(ii)
        row=zeros(1,K);row(ii(k-1))=1;row(ii(k))=-1;
        A(end+1,:)=row;b(end+1,1)=-ng; %#ok<AGROW>
    end
end
% Phase lower/upper total times
[A,b]=addSum(A,b,idx.th,cfg.homing_min(2),K,1);
[A,b]=addSum(A,b,idx.th,cfg.homing_min(1),K,-1);
[A,b]=addSum(A,b,idx.tc,cfg.closing_min(2),K,1);
[A,b]=addSum(A,b,idx.tc,cfg.closing_min(1),K,-1);
[A,b]=addSum(A,b,[idx.th idx.tc], ...
    cap-(cfg.hold1_min+cfg.hold2_min+cfg.final_min),K,1);
end
function [A,b]=addSum(A,b,inds,val,K,sgn)
row=zeros(1,K);row(inds)=sgn;
A(end+1,:)=row;b(end+1,1)=sgn*val;
end
function z=initialGuess(Nh,Nc,cap,cfg,idx)
fh=drawFrac(Nh,cfg.minWaypointGap);
fc=drawFrac(Nc,cfg.minWaypointGap);
% Draw feasible phase durations by distributing slack from minima
low=[max(cfg.homing_min(1),Nh*cfg.minLeg_min), ...
     max(cfg.closing_min(1),Nc*cfg.minLeg_min)];
fixed=cfg.hold1_min+cfg.hold2_min+cfg.final_min;
avail=max(0,cap-fixed-sum(low));
headroom=[cfg.homing_min(2)-low(1),cfg.closing_min(2)-low(2)];
extra=min(headroom,avail*rand(1,2));
if sum(extra)>avail && sum(extra)>0,extra=extra*(avail/sum(extra));end
Th=low(1)+extra(1);Tc=low(2)+extra(2);
th=splitTimes(Th,Nh,cfg.minLeg_min);
tc=splitTimes(Tc,Nc,cfg.minLeg_min);
z=[mean(cfg.hp1_m);mean(cfg.hp2_m);fh;fc;th;tc];
assert(numel(z)==max(idx.tc));
end
function f=drawFrac(N,gap)
if N==1,f=zeros(0,1);return;end
w=rand(N,1);w=w/sum(w);
steps=gap+(1-N*gap)*w;
f=cumsum(steps(1:end-1));
end
function t=splitTimes(T,N,minLeg)
w=rand(N,1);w=w/sum(w);
t=minLeg+(T-N*minLeg)*w;
end
function j=objective(z,Nh,Nc,n,cfg,idx)
r=evaluate(z,Nh,Nc,n,cfg,idx,cfg.nSamples,false);
if r.valid,j=r.dv;else,j=1e5;end
end
function [c,ceq]=nonlinear(z,Nh,Nc,n,cfg,idx)
r=evaluate(z,Nh,Nc,n,cfg,idx,cfg.nSamples,false);
if ~r.valid,c=ones(7,1)*1e3;ceq=[];return;end
accelViolation=-1;
if isfinite(cfg.maxPeakAccel_mps2)
    accelViolation=r.accF-cfg.maxPeakAccel_mps2;
end
c=[r.spH-cfg.maxSpeed_homing_mps; ...
   r.spC-cfg.maxSpeed_closing_mps; ...
   r.spF-cfg.maxSpeed_final_mps; ...
   cfg.KOZ_m-min(r.minH,r.minC); ...
   r.maxH-cfg.homingExcursionFactor*norm(cfg.x0(1:3)); ...
   r.maxC-cfg.closingExcursionFactor*z(1); ...
   accelViolation];
ceq=[];
end
function pass=isFeasible(r,cap,cfg)
pass=r.valid && r.tof<=cap+1e-5 && ...
    r.spH<=cfg.maxSpeed_homing_mps+1e-5 && ...
    r.spC<=cfg.maxSpeed_closing_mps+1e-5 && ...
    r.spF<=cfg.maxSpeed_final_mps+1e-5 && ...
    min(r.minH,r.minC)>=cfg.KOZ_m-1e-4 && ...
    r.maxH<=cfg.homingExcursionFactor*norm(cfg.x0(1:3))+1e-4 && ...
    r.maxC<=cfg.closingExcursionFactor*r.hp1+1e-4 && ...
    r.accF<=cfg.maxPeakAccel_mps2+1e-5 && ...
    r.posError<=cfg.terminalTol_m && ...
    r.velError<=cfg.terminalVelTol_mps;
end
function r=evaluate(z,Nh,Nc,n,cfg,idx,Nsamples,odeCheck)
r=struct('valid',false,'dv',inf,'tof',inf, ...
    'hdv',inf,'cdv',inf,'fdv',inf,'spH',inf,'spC',inf,'spF',inf, ...
    'minH',0,'minC',0,'maxH',inf,'maxC',inf,'accF',inf, ...
    'hp1',z(1),'posError',inf,'velError',inf);
if any(~isfinite(z)) || any(z([idx.th idx.tc])<=0),return;end
hp1=z(1);hp2=z(2);
if hp1<=hp2 || hp2<=cfg.KOZ_m || hp2<=cfg.finalStandoff_m,return;end
fh=[z(idx.h);1];fc=[z(idx.c);1];
if any(diff([0;fh])<=0) || any(diff([0;fc])<=0),return;end
xh=cfg.x0(:);
x1=[0;-hp1;0;0;0;0];x2=[0;-hp2;0;0;0;0];
r.hdv=0;r.cdv=0;r.spH=0;r.spC=0;
r.minH=inf;r.minC=inf;r.maxH=0;r.maxC=0;
r.posError=0;r.velError=0;
% Each phase consists of fly-through HCW arcs. Only the phase endpoint
% is braked to zero relative velocity. One correction burn at each
% intermediate waypoint; no artificial arrival/departure burn pair.
[h,validH]=flyPhase(cfg.x0(:),x1,fh,z(idx.th),n,Nsamples,odeCheck);
if ~validH,return;end
[c,validC]=flyPhase(x1,x2,fc,z(idx.tc),n,Nsamples,odeCheck);
if ~validC,return;end
r.hdv=h.dv;r.cdv=c.dv;
r.spH=h.sp;r.spC=c.sp;
r.minH=h.minR;r.minC=c.minR;
r.maxH=h.maxR;r.maxC=c.maxR;
r.posError=max(h.pe,c.pe);
r.velError=max(h.ve,c.ve);
% V-bar final continuous quintic: start from stationary HP2, end at
% stationary capture standoff, with exact HCW feedforward acceleration.
[r.fdv,r.spF,r.accF]=finalQuintic(hp2,cfg.finalStandoff_m, ...
    cfg.final_min*60,n,Nsamples);
r.dv=r.hdv+r.cdv+r.fdv+ ...
    cfg.hold1_DVrate_mps_per_min*cfg.hold1_min + ...
    cfg.hold2_DVrate_mps_per_min*cfg.hold2_min;
r.tof=sum(z([idx.th idx.tc]))+cfg.hold1_min+cfg.hold2_min+cfg.final_min;
r.valid=all(isfinite([r.dv r.tof r.spH r.spC r.spF r.minH r.minC]));
end
function [p,ok]=flyPhase(xStart,xTarget,f,tMinutes,n,Nsamples,odeCheck)
% Reproduce each fly-through waypoint with a velocity impulse.
% The only arrival stop occurs at the FINAL holding point.
p=struct('dv',0,'sp',0,'minR',inf,'maxR',0,'pe',0,'ve',0);
ok=false;
rStart=xStart(1:3);
vBefore=xStart(4:6);
Tsum=numel(tMinutes);
for k=1:Tsum
    rk=xStart(1:3)+(xTarget(1:3)-xStart(1:3))*f(k);
    Tk=60*tMinutes(k);
    [A,B,C,D]=phi(n,Tk);
    if rcond(B)<1e-10,return;end
    vPlus=B\(rk-A*rStart);
    if any(~isfinite(vPlus)),return;end
    correction=vPlus-vBefore;
    p.dv=p.dv+norm(correction);
    vMinus=C*rStart+D*vPlus;
    % Verify coast with state-transition matrices and optionally ode45.
    for j=1:Nsamples
        tau=Tk*(j-1)/(Nsamples-1);
        [At,Bt,Ct,Dt]=phi(n,tau);
        rr=At*rStart+Bt*vPlus;
        vv=Ct*rStart+Dt*vPlus;
        p.minR=min(p.minR,norm(rr));
        p.maxR=max(p.maxR,norm(rr));
        p.sp=max(p.sp,norm(vv));
    end
    posErr=norm(A*rStart+B*vPlus-rk);
    velErr=0;
    if odeCheck
        options=odeset('RelTol',1e-10,'AbsTol',1e-11);
        [~,X]=ode45(@(~,q)[q(4:6);3*n^2*q(1)+2*n*q(5); ...
              -2*n*q(4);-n^2*q(3)],linspace(0,Tk,50), ...
              [rStart;vPlus],options);
        posErr=norm(X(end,1:3)'-rk);
        velErr=norm(X(end,4:6)'-vMinus);
        p.minR=min(p.minR,min(vecnorm(X(:,1:3),2,2)));
        p.maxR=max(p.maxR,max(vecnorm(X(:,1:3),2,2)));
        p.sp=max(p.sp,max(vecnorm(X(:,4:6),2,2)));
    end
    p.pe=max(p.pe,posErr);
    p.ve=max(p.ve,velErr);
    rStart=rk;
    vBefore=vMinus; % fly through: do NOT brake here
end
% Braking impulse exists ONLY at the designated final holding point.
p.dv=p.dv+norm(xTarget(4:6)-vBefore);
ok=all(isfinite([p.dv p.sp p.minR p.maxR p.pe p.ve]));
end
function [dv,sp,acc]=finalQuintic(hp2,standoff,T,n,N)
t=linspace(0,T,N);s=t/T;
p=10*s.^3-15*s.^4+6*s.^5;
pd=(30*s.^2-60*s.^3+30*s.^4)/T;
pdd=(60*s-180*s.^2+120*s.^3)/T^2;
y0=-hp2;dy=hp2-standoff;
y=y0+dy*p;vy=dy*pd;ay=dy*pdd;
% x=0 and z=0. HCW feedforward radial: -2*n*yDot
ux=-2*n*vy;uy=ay;
accmag=sqrt(ux.^2+uy.^2);
dv=trapz(t,accmag);sp=max(abs(vy));acc=max(accmag);
end
function [A,B,C,D]=phi(n,T)
nt=n*T;c=cos(nt);s=sin(nt);
A=[4-3*c,0,0;6*(s-nt),1,0;0,0,c];
B=[s,2*(1-c),0;2*(c-1),4*s-3*nt,0;0,0,s]/n;
C=[3*n*s,0,0;6*n*(c-1),0,0;0,0,-n*s];
D=[c,2*s,0;-2*s,4*c-3,0;0,0,c];
end
function mask=paretoMask(t,d)
[~,o]=sortrows([t(:) d(:)],[1 2]);mask=false(numel(t),1);best=inf;
for k=1:numel(o)
    i=o(k);if d(i)<best-1e-9,mask(i)=true;best=d(i);end
end
end
function [H,C]=arcTables(z,Nh,Nc,n,cfg)
p=3;fh=[z(p:p+Nh-2);1];p=p+Nh-1;
fc=[z(p:p+Nc-2);1];p=p+Nc-1;
th=z(p:p+Nh-1);p=p+Nh;tc=z(p:p+Nc-1);
H=makeTable(cfg.x0(:),[0;-z(1);0;0;0;0],fh,th,n,'Homing');
C=makeTable([0;-z(1);0;0;0;0],[0;-z(2);0;0;0;0],fc,tc,n,'Closing');
end
function T=makeTable(x0,xtarget,frac,times,n,label)
% BurnStart at each arc start is the velocity correction needed to
% continue the sequence. BurnEnd is zero except at the final hold.
N=numel(times);dat=zeros(N,17);
r0=x0(1:3);rStart=r0;vBefore=x0(4:6);
for k=1:N
    rEnd=r0+(xtarget(1:3)-r0)*frac(k);
    [A,B,C,D]=phi(n,times(k)*60);
    if rcond(B)<1e-10
        error('Cannot export near-singular HCW arc %d',k);
    end
    vPlus=B\(rEnd-A*rStart);
    burnStart=vPlus-vBefore;
    vMinus=C*rStart+D*vPlus;
    burnEnd=zeros(3,1);
    if k==N
        burnEnd=xtarget(4:6)-vMinus;
    end
    dat(k,:)=[k rStart' rEnd' times(k) burnStart' burnEnd' vMinus'];
    rStart=rEnd;
    vBefore=vMinus;
end
T=array2table(dat,'VariableNames',{'Arc','Xstart_m','Ystart_m','Zstart_m', ...
    'Xend_m','Yend_m','Zend_m','TOF_min','BurnStartX_mps', ...
    'BurnStartY_mps','BurnStartZ_mps','BurnEndX_mps','BurnEndY_mps', ...
    'BurnEndZ_mps','FlyThroughVx_mps','FlyThroughVy_mps','FlyThroughVz_mps'});
T.Phase=repmat(string(label),N,1);
T.BurnStartMagnitude_mps=sqrt(T.BurnStartX_mps.^2+T.BurnStartY_mps.^2+T.BurnStartZ_mps.^2);
T.BurnEndMagnitude_mps=sqrt(T.BurnEndX_mps.^2+T.BurnEndY_mps.^2+T.BurnEndZ_mps.^2);
end
