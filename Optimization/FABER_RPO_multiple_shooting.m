function R = FABER_RPO_multiple_shooting(user)
% FABER_RPO_MULTIPLE_SHOOTING
% Joint homing + closing: direct multiple shooting using exact HCW STM.
% Fly-through velocities are explicit variables; burn corrections are
% velocity jumps. Only HP1 and HP2 have zero terminal relative velocity.
% Stationary V-bar waypoint architecture (NOT pure V-bar burns).
% TOF per arc is free within phase limits, unrelated to orbital fractions.
% No perturbations, finite thrust, active abort, or certified safety.
% Requires Optimization Toolbox / fmincon.
%
% R = FABER_RPO_multiple_shooting();
% R = FABER_RPO_multiple_shooting(struct('Nmax',2,'nStarts',2));
% Result includes all verified candidates, Pareto, and individual burns.

c = defaults();
if nargin>0 && ~isempty(user)
    f=fieldnames(user);
    for j=1:numel(f),c.(f{j})=user.(f{j});end
end
assert(exist('fmincon','file')==2,'Optimization Toolbox fmincon required.');
assert(numel(c.x0)==6 && c.Nmax>=1 && c.Nmax<=4);
assert(c.hp1_m(1)>c.hp2_m(2) && c.hp2_m(1)>c.KOZ_m);
assert(c.gap>0 && c.gap*c.Nmax<1);
assert(c.minLeg_min*c.Nmax<=min(c.homing_min(2),c.closing_min(2)));
assert(c.finalStandoff_m>0 && c.finalStandoff_m<c.KOZ_m);
n=sqrt(c.mu_km3s2/(c.Re_km+c.h_km)^3);
rng(c.seed,'twister');
opt=optimoptions('fmincon','Algorithm','sqp','Display','off', ...
    'MaxIterations',c.maxIter,'MaxFunctionEvaluations',c.maxEvals, ...
    'ConstraintTolerance',1e-6,'StepTolerance',1e-9, ...
    'OptimalityTolerance',1e-6);
S=struct('Nh',{},'Nc',{},'TOFcap_min',{},'DV_mps',{},'TOF_min',{}, ...
    'HomingDV_mps',{},'ClosingDV_mps',{},'FinalDV_mps',{}, ...
    'HP1_m',{},'HP2_m',{},'MaxSpeed_mps',{},'MinRange_m',{}, ...
    'MaxTerminalError_m',{},'MaxODEError_m',{},'z',{});

for nh=1:c.Nmax
 for nc=1:c.Nmax
  for cap=c.tofCaps_min(:)'
   [lb,ub,A,b,idx,possible]=bounds(nh,nc,cap,c);
   if ~possible,continue;end
   winner=[];
   for trial=1:c.nStarts
    z0=guess(nh,nc,cap,c,idx,n);
    if isempty(z0),continue;end
    try
     [z,~,flag]=fmincon(@(v)cost(v,nh,nc,idx,c,n),z0, ...
         A,b,[],[],lb,ub,@(v)constraints(v,nh,nc,idx,c,n),opt);
    catch ME
     warning('Nh=%d Nc=%d cap=%g: %s',nh,nc,cap,ME.message);
     continue;
    end
    if flag<=0,continue;end
    q=simulate(z,nh,nc,idx,c,n,c.validationSamples,true);
    if ~accepted(q,c,cap),continue;end
    if isempty(winner) || q.dv<winner.q.dv
        winner=struct('q',q,'z',z);
    end
   end
   if ~isempty(winner)
    q=winner.q;z=winner.z;
    S(end+1)=struct('Nh',nh,'Nc',nc,'TOFcap_min',cap, ...
        'DV_mps',q.dv,'TOF_min',q.tof,'HomingDV_mps',q.h.dv, ...
        'ClosingDV_mps',q.cl.dv,'FinalDV_mps',q.finalDV, ...
        'HP1_m',z(1),'HP2_m',z(2), ...
        'MaxSpeed_mps',max([q.h.maxSpeed,q.cl.maxSpeed,q.finalSpeed]), ...
        'MinRange_m',min([q.h.minRange,q.cl.minRange]), ...
        'MaxTerminalError_m',q.maxPosError, ...
        'MaxODEError_m',q.maxODEError,'z',z); %#ok<AGROW>
   else
    fprintf('No verified candidate Nh=%d Nc=%d cap=%.1f min\n',nh,nc,cap);
   end
  end
 end
end
R=struct('cfg',c,'solutions',table(),'pareto',table(), ...
    'bestHomingArcs',table(),'bestClosingArcs',table(),'raw',S);
if isempty(S),warning('No verified solutions. Increase starts or relax constraints.');return;end
R.solutions=struct2table(rmfield(S,'z'));
mask=pareto(R.solutions.TOF_min,R.solutions.DV_mps);
R.pareto=sortrows(R.solutions(mask,:),'TOF_min');
[~,j]=min([S.DV_mps]);best=S(j);
[~,~,~,~,idx]=bounds(best.Nh,best.Nc,best.TOFcap_min,c);
q=simulate(best.z,best.Nh,best.Nc,idx,c,n,c.validationSamples,true);
R.bestHomingArcs=arctable(q.h,'Homing');
R.bestClosingArcs=arctable(q.cl,'Closing');
R.best=struct('design',best,'validation',q);
disp(R.solutions);disp(R.pareto);
if c.plot
 figure('Color','w','Name','FABER multiple shooting Pareto');
 scatter(R.solutions.TOF_min,R.solutions.DV_mps,35, ...
    R.solutions.Nh+R.solutions.Nc,'filled');hold on;grid on;
 plot(R.pareto.TOF_min,R.pareto.DV_mps,'k-o','LineWidth',1.4);
 xlabel('RPO time [min]');ylabel('Delta V [m/s]');
 title('Joint homing/closing: direct multiple shooting (nominal)');
 colorbar;
end
if c.saveCSV
 writetable(R.solutions,'FABER_MS_solutions.csv');
 writetable(R.pareto,'FABER_MS_pareto.csv');
 writetable(R.bestHomingArcs,'FABER_MS_best_homing.csv');
 writetable(R.bestClosingArcs,'FABER_MS_best_closing.csv');
end
end

function c=defaults()
c.mu_km3s2=398600.4418;c.Re_km=6371;c.h_km=457.25;
c.x0=[50;-2000;20;.01;.04;-.01];
c.Nmax=2;c.nStarts=2;c.seed=42;c.tofCaps_min=[170 210 260];
c.hp1_m=[800 1200];c.hp2_m=[100 250];c.finalStandoff_m=5;
c.homing_min=[30 100];c.closing_min=[20 130];
c.hold1_min=30;c.hold2_min=10;c.final_min=10;
c.minLeg_min=5;c.gap=0.05;
c.maxSpeed_homing_mps=1;c.maxSpeed_closing_mps=1;
c.maxSpeed_final_mps=1;c.KOZ_m=30;
c.homingExcursionFactor=1.20;c.closingExcursionFactor=1.20;
c.maxPeakAccel_mps2=inf;
c.hold1_DVrate_mps_per_min=0;c.hold2_DVrate_mps_per_min=0;
c.nSamples=21;c.validationSamples=201;
c.maxIter=200;c.maxEvals=15000;
c.terminalTol_m=0.01;c.terminalVelTol_mps=1e-5;
c.odeTol_m=0.01;c.plot=true;c.saveCSV=true;
end

function [lb,ub,A,b,I,possible]=bounds(nh,nc,cap,c)
% z = [hp1 hp2 fractionsH fractionsC timesH timesC
%       post-impulse velocitiesH (m/s) post-impulse velocitiesC (m/s)]
p=3;I.fh=p:p+nh-2;p=p+nh-1;
I.fc=p:p+nc-2;p=p+nc-1;
I.th=p:p+nh-1;p=p+nh;
I.tc=p:p+nc-1;p=p+nc;
I.vh=p:p+3*nh-1;p=p+3*nh;
I.vc=p:p+3*nc-1;p=p+3*nc;
K=p-1;
lb=[c.hp1_m(1);c.hp2_m(1);repmat(c.gap,nh+nc-2,1); ...
    repmat(c.minLeg_min,nh+nc,1);-5*ones(3*(nh+nc),1)];
ub=[c.hp1_m(2);c.hp2_m(2);repmat(1-c.gap,nh+nc-2,1); ...
    repmat(max(c.homing_min(2),c.closing_min(2)),nh+nc,1); ...
    5*ones(3*(nh+nc),1)];
A=zeros(0,K);b=zeros(0,1);
for ind={I.fh,I.fc}
 ids=ind{1};
 for k=2:numel(ids)
  row=zeros(1,K);row(ids(k-1))=1;row(ids(k))=-1;
  A(end+1,:)=row;b(end+1,1)=-c.gap; %#ok<AGROW>
 end
end
[A,b]=sumRow(A,b,I.th,c.homing_min(2),1,K);
[A,b]=sumRow(A,b,I.th,c.homing_min(1),-1,K);
[A,b]=sumRow(A,b,I.tc,c.closing_min(2),1,K);
[A,b]=sumRow(A,b,I.tc,c.closing_min(1),-1,K);
fixed=c.hold1_min+c.hold2_min+c.final_min;
[A,b]=sumRow(A,b,[I.th I.tc],cap-fixed,1,K);
possible=(cap>=fixed+max(nh*c.minLeg_min,c.homing_min(1)) ...
    +max(nc*c.minLeg_min,c.closing_min(1)));
end
function [A,b]=sumRow(A,b,ids,val,sign,K)
r=zeros(1,K);r(ids)=sign;A(end+1,:)=r;b(end+1,1)=sign*val;
end
function z=guess(nh,nc,cap,c,I,n)
fixed=c.hold1_min+c.hold2_min+c.final_min;
lo=[max(c.homing_min(1),nh*c.minLeg_min), ...
    max(c.closing_min(1),nc*c.minLeg_min)];
remaining=cap-fixed-sum(lo);
if remaining<-1e-8,z=[];return;end
fraction=rand();
extraH=min([fraction*remaining,c.homing_min(2)-lo(1)]);
extraC=min([(remaining-extraH)*rand(),c.closing_min(2)-lo(2)]);
th=split(lo(1)+extraH,nh,c.minLeg_min);
tc=split(lo(2)+extraC,nc,c.minLeg_min);
fh=fractions(nh,c.gap);fc=fractions(nc,c.gap);
hp1=c.hp1_m(1)+rand()*diff(c.hp1_m);
hp2=c.hp2_m(1)+rand()*diff(c.hp2_m);
Vh=seedVel(c.x0(:),[0;-hp1;0;0;0;0],[fh;1],th,n);
Vc=seedVel([0;-hp1;0;0;0;0],[0;-hp2;0;0;0;0],[fc;1],tc,n);
if isempty(Vh)||isempty(Vc),z=[];return;end
z=[hp1;hp2;fh;fc;th;tc;Vh(:);Vc(:)];
assert(numel(z)==max(I.vc));
end
function x=fractions(N,g)
if N==1,x=zeros(0,1);return;end
w=rand(N,1);w=w/sum(w);w=g+(1-N*g)*w;
x=cumsum(w(1:end-1));
end
function t=split(T,N,m)
w=rand(N,1);w=w/sum(w);t=m+(T-N*m)*w;
end
function V=seedVel(x0,xf,f,t,n)
V=zeros(3,numel(t));r=x0(1:3);start=x0(1:3);
for k=1:numel(t)
 dest=x0(1:3)+f(k)*(xf(1:3)-x0(1:3));
 [A,B]=stm(n,t(k)*60);
 if rcond(B)<1e-9,V=[];return;end
 V(:,k)=B\(dest-A*r);r=dest;
end
end

function f=cost(z,nh,nc,I,c,n)
q=simulate(z,nh,nc,I,c,n,3,false);
if q.valid,f=q.dv;else,f=1e5;end
end
function [g,h]=constraints(z,nh,nc,I,c,n)
q=simulate(z,nh,nc,I,c,n,c.nSamples,false);
if ~q.valid,g=ones(7,1)*1000;h=ones(3*(nh+nc),1)*1000;return;end
g=[q.h.maxSpeed-c.maxSpeed_homing_mps; ...
   q.cl.maxSpeed-c.maxSpeed_closing_mps; ...
   q.finalSpeed-c.maxSpeed_final_mps; ...
   c.KOZ_m-min(q.h.minRange,q.cl.minRange); ...
   q.h.maxRange-c.homingExcursionFactor*norm(c.x0(1:3)); ...
   q.cl.maxRange-c.closingExcursionFactor*z(1); ...
   q.finalAccel-c.maxPeakAccel_mps2];
if isinf(c.maxPeakAccel_mps2),g(end)=-1;end
% Position defects are in km to balance gradient scaling against burn m/s.
h=q.defects(:)/1000;
end

function q=simulate(z,nh,nc,I,c,n,N,odeCheck)
q=struct('valid',false,'dv',inf,'tof',inf,'h',[], 'cl',[], ...
    'finalDV',inf,'finalSpeed',inf,'finalAccel',inf, ...
    'maxPosError',inf,'maxODEError',inf,'defects',[]);
if any(~isfinite(z)) || any(z([I.th I.tc])<=0),return;end
hp1=z(1);hp2=z(2);
if hp1<=hp2 || hp2<=c.KOZ_m || hp2<=c.finalStandoff_m,return;end
fh=[z(I.fh);1];fc=[z(I.fc);1];
if any(diff([0;fh])<=0)||any(diff([0;fc])<=0),return;end
vh=reshape(z(I.vh),3,nh);vc=reshape(z(I.vc),3,nc);
q.h=phase(c.x0(:),[0;-hp1;0;0;0;0],fh,z(I.th),vh,n,N,odeCheck);
q.cl=phase([0;-hp1;0;0;0;0],[0;-hp2;0;0;0;0], ...
    fc,z(I.tc),vc,n,N,odeCheck);
[q.finalDV,q.finalSpeed,q.finalAccel]=finalQuintic( ...
    hp2,c.finalStandoff_m,c.final_min*60,n,max(N,21));
q.dv=q.h.dv+q.cl.dv+q.finalDV+ ...
    c.hold1_min*c.hold1_DVrate_mps_per_min+ ...
    c.hold2_min*c.hold2_DVrate_mps_per_min;
q.tof=sum(z([I.th I.tc]))+c.hold1_min+c.hold2_min+c.final_min;
q.defects=[q.h.defects;q.cl.defects];
q.maxPosError=max(abs(q.defects));
q.maxODEError=max(q.h.maxODEError,q.cl.maxODEError);
q.valid=all(isfinite([q.dv q.tof q.maxPosError q.maxODEError ...
    q.h.maxSpeed q.cl.maxSpeed q.h.minRange q.cl.minRange]));
end
function p=phase(x0,xf,f,t,vPlus,n,N,odeCheck)
K=numel(t);
p=struct('dv',0,'maxSpeed',0,'minRange',inf,'maxRange',0, ...
  'defects',zeros(3*K,1),'maxODEError',0, ...
  'rStart',zeros(3,K),'rEnd',zeros(3,K), ...
  'vBefore',zeros(3,K),'vPlus',vPlus,'vMinus',zeros(3,K), ...
  'burns',zeros(3,K),'lastBrake',zeros(3,1),'times',t);
rStart=x0(1:3);vBefore=x0(4:6);r0=x0(1:3);
for k=1:K
 rEnd=r0+f(k)*(xf(1:3)-r0);
 [A,B,C,D]=stm(n,t(k)*60);
 v=vPlus(:,k);
 rActual=A*rStart+B*v;
 vMinus=C*rStart+D*v;
 burn=v-vBefore;
 p.dv=p.dv+norm(burn);
 p.defects(3*k-2:3*k)=rActual-rEnd;
 p.rStart(:,k)=rStart;p.rEnd(:,k)=rEnd;
 p.vBefore(:,k)=vBefore;p.vMinus(:,k)=vMinus;p.burns(:,k)=burn;
 for j=1:N
  [At,Bt,Ct,Dt]=stm(n,t(k)*60*(j-1)/(N-1));
  rr=At*rStart+Bt*v;vv=Ct*rStart+Dt*v;
  p.maxSpeed=max(p.maxSpeed,norm(vv));
  p.minRange=min(p.minRange,norm(rr));
  p.maxRange=max(p.maxRange,norm(rr));
 end
 if odeCheck
  opt=odeset('RelTol',1e-10,'AbsTol',1e-11);
  [~,Y]=ode45(@(~,x)hcwODE(x,n),linspace(0,60*t(k),60),[rStart;v],opt);
  p.maxODEError=max(p.maxODEError,norm(Y(end,1:3)'-rActual));
  p.maxSpeed=max(p.maxSpeed,max(vecnorm(Y(:,4:6),2,2)));
  p.minRange=min(p.minRange,min(vecnorm(Y(:,1:3),2,2)));
  p.maxRange=max(p.maxRange,max(vecnorm(Y(:,1:3),2,2)));
 end
 % Next arc begins at prescribed fly-through position; position mismatch
 % appears explicitly in equality constraints and acceptance checks.
 rStart=rEnd;vBefore=vMinus;
end
p.lastBrake=xf(4:6)-vBefore;
p.dv=p.dv+norm(p.lastBrake);
end
function dx=hcwODE(x,n)
dx=[x(4:6);3*n*n*x(1)+2*n*x(5);-2*n*x(4);-n*n*x(3)];
end
function [A,B,C,D]=stm(n,T)
a=n*T;s=sin(a);c=cos(a);
A=[4-3*c,0,0;6*(s-a),1,0;0,0,c];
B=[s,2*(1-c),0;2*(c-1),4*s-3*a,0;0,0,s]/n;
C=[3*n*s,0,0;6*n*(c-1),0,0;0,0,-n*s];
D=[c,2*s,0;-2*s,4*c-3,0;0,0,c];
end
function [dv,sp,ac]=finalQuintic(d0,d1,T,n,N)
t=linspace(0,T,N);s=t/T;dy=d0-d1;
pd=(30*s.^2-60*s.^3+30*s.^4)/T;
pdd=(60*s-180*s.^2+120*s.^3)/T^2;
vy=dy*pd;uy=dy*pdd;ux=-2*n*vy;
a=sqrt(ux.^2+uy.^2);
dv=trapz(t,a);sp=max(abs(vy));ac=max(a);
end
function ok=accepted(q,c,cap)
ok=q.valid && q.maxPosError<=c.terminalTol_m && ...
 q.maxODEError<=c.odeTol_m && q.tof<=cap+1e-5 && ...
 q.h.maxSpeed<=c.maxSpeed_homing_mps+1e-5 && ...
 q.cl.maxSpeed<=c.maxSpeed_closing_mps+1e-5 && ...
 q.finalSpeed<=c.maxSpeed_final_mps+1e-5 && ...
 min(q.h.minRange,q.cl.minRange)>=c.KOZ_m-1e-4 && ...
 q.h.maxRange<=c.homingExcursionFactor*norm(c.x0(1:3))+1e-4 && ...
 q.cl.maxRange<=c.closingExcursionFactor*q.cl.rStart(2,1)*(-1)+1e-4 && ...
 q.finalAccel<=c.maxPeakAccel_mps2+1e-5;
end
function T=arctable(p,name)
K=numel(p.times);
arc=(1:K)';
T=table(repmat(string(name),K,1),arc,p.times(:), ...
 p.rStart(1,:)',p.rStart(2,:)',p.rStart(3,:)', ...
 p.rEnd(1,:)',p.rEnd(2,:)',p.rEnd(3,:)', ...
 p.burns(1,:)',p.burns(2,:)',p.burns(3,:)', ...
 vecnorm(p.burns,2,1)',p.vMinus(1,:)',p.vMinus(2,:)', ...
 p.vMinus(3,:)', ...
 'VariableNames',{'Phase','Arc','TOF_min','Xstart_m','Ystart_m','Zstart_m', ...
 'Xend_m','Yend_m','Zend_m','BurnX_mps','BurnY_mps','BurnZ_mps', ...
 'BurnMagnitude_mps','ArrivalVx_mps','ArrivalVy_mps','ArrivalVz_mps'});
T.TerminalBrake_mps=zeros(K,1);T.TerminalBrake_mps(end)=norm(p.lastBrake);
end
function m=pareto(t,d)
[~,order]=sortrows([t(:),d(:)],[1,2]);m=false(numel(t),1);best=inf;
for k=1:numel(order)
 j=order(k);if d(j)<best-1e-9,m(j)=true;best=d(j);end
end
end
