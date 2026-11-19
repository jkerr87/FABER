function R = FABER_safety_MC_core(events,cfg)
% FABER_SAFETY_MC_CORE  Conditional failure-injection / finite-burn CAM screening.
% Event fields: xPre [6x1] (m,m/s), dvNom [3x1] (m/s or NaN for cutoff),
% phase, vehicle, index. Optional: hazardsCenters_m [Nh x 3], mass_kg.
% Event failures are sampled CONDITIONALLY, not at real occurrence rates.
% Candidate CAM uses constant FORCE in a fixed LVLH direction, propagated
% through the HCW equations DURING the finite burn. Not a flight controller.
assert(~isempty(events),'No failure events supplied');
assert(isfinite(cfg.camForce_N)&&cfg.camForce_N>0);
assert(isfinite(cfg.vehicleMass_kg)&&cfg.vehicleMass_kg>0);
assert(cfg.screen.collisionRadius_m>0);
assert(cfg.screen.protectedRadius_m>cfg.screen.collisionRadius_m);
assert(cfg.screen.dt_s>0&&cfg.screen.passiveHorizon_s>0);
assert(cfg.detectionDelay_s>=0);
rng(cfg.seed,'twister');
N=cfg.N;
Trial=(1:N)';EventIndex=zeros(N,1);Module=zeros(N,1);
Vehicle=strings(N,1);Phase=strings(N,1);Case=strings(N,1);
PassiveMinClearance_m=nan(N,1);PassiveBreach=false(N,1);
FirstCollision_min=nan(N,1);CAMFeasible=false(N,1);
CAM_DV_mps=nan(N,1);CAMBurn_min=nan(N,1);
TimeToSafeShell_min=nan(N,1);CAMTotalTime_min=nan(N,1);
PostCAMMinClearance_m=nan(N,1);Status=strings(N,1);
isBurn=arrayfun(@(e)all(isfinite(e.dvNom(:))),events);
for k=1:N
    % Alternate burn and cutoff events when available.
    pool=find(isBurn==mod(k,2));
    if isempty(pool),pool=1:numel(events);end
    ei=pool(randi(numel(pool)));e=events(ei);
    EventIndex(k)=ei;Vehicle(k)=string(e.vehicle);
    Phase(k)=string(e.phase);Module(k)=e.index;
    hazards=cfg.hazards;
    if isfield(e,'hazardsCenters_m') && ~isempty(e.hazardsCenters_m)
        hazards.centers_m=e.hazardsCenters_m;
        hazards.radii_m=zeros(size(hazards.centers_m,1),1);
    end
    if isfield(e,'hazardsSegmentsA_m')
        hazards.segmentsA_m=e.hazardsSegmentsA_m;
        hazards.segmentsB_m=e.hazardsSegmentsB_m;
        hazards.segmentRadii_m=e.hazardsSegmentRadii_m;
    end
    mass=cfg.vehicleMass_kg;
    if isfield(e,'mass_kg') && isfinite(e.mass_kg)&&e.mass_kg>0
        mass=e.mass_kg;
    end
    x=e.xPre(:);dv=e.dvNom(:);
    x(1:3)=x(1:3)+cfg.sigmaPos_m*randn(3,1);
    x(4:6)=x(4:6)+cfg.sigmaVel_mps*randn(3,1);
    if all(isfinite(dv))
        z=rand();
        if z<cfg.caseWeights(1)
            Case(k)="Missed burn";actual=zeros(3,1);
        elseif z<sum(cfg.caseWeights(1:2))
            Case(k)="Thrust magnitude error";
            actual=dv*(1+cfg.sigmaThrustFrac*randn());
        else
            Case(k)="Pointing error";
            actual=tilt(dv,cfg.sigmaPointing_rad*randn());
        end
        x(4:6)=x(4:6)+actual;
    else
        Case(k)="Loss of continuous control";
        % xPre describes nominal state at the chosen cutoff instant.
    end
    [pc,breach,tb]=driftCheck(x,0,cfg.screen.passiveHorizon_s, ...
        cfg.screen.dt_s,cfg.n,hazards,cfg.screen.collisionRadius_m);
    PassiveMinClearance_m(k)=pc; PassiveBreach(k)=breach;
    FirstCollision_min(k)=tb/60;
    % Detection latency must itself remain collision-free.
    [~,badDelay,~]=driftCheck(x,0,cfg.detectionDelay_s, ...
        cfg.screen.dt_s,cfg.n,hazards,cfg.screen.collisionRadius_m);
    if badDelay,Status(k)="Collision before response";continue;end
    xd=hcwPhi(cfg.n,cfg.detectionDelay_s)*x;
    bestDV=inf;
    for mag=cfg.camDVgrid_mps(:)'
        if mag>cfg.camMaxDV_mps,continue;end
        for d=1:size(cfg.camDirections,2)
            dir=cfg.camDirections(:,d);
            % Nominal delta-V equivalent = thrust/mass * burn duration.
            Tburn=mag*mass/cfg.camForce_N;
            if Tburn>cfg.camMaxBurn_s,continue;end
            [xb,badBurn]=finiteBurn(xd,cfg.n,dir*cfg.camForce_N/mass, ...
                Tburn,cfg.screen.burnStep_s,hazards, ...
                cfg.screen.collisionRadius_m);
            if badBurn,continue;end
            [~,badDrift,~]=driftCheck(xb,0,cfg.screen.passiveHorizon_s, ...
                cfg.screen.dt_s,cfg.n,hazards,cfg.screen.collisionRadius_m);
            if badDrift,continue;end
            [safe,tExit,cMin]=checkShell(xb,cfg.n,hazards,cfg.screen);
            if safe && mag<bestDV
                bestDV=mag;CAMBurn_min(k)=Tburn/60;
                TimeToSafeShell_min(k)=(cfg.detectionDelay_s+Tburn+tExit)/60;
                CAMTotalTime_min(k)=(cfg.detectionDelay_s+Tburn)/60;
                PostCAMMinClearance_m(k)=cMin;
            end
        end
    end
    if isfinite(bestDV)
        CAMFeasible(k)=true;CAM_DV_mps(k)=bestDV;
        Status(k)="Candidate escape passed screening";
    else
        Status(k)="No feasible candidate escape";
    end
end
R.trials=table(Trial,EventIndex,Module,Vehicle,Phase,Case, ...
    PassiveMinClearance_m,PassiveBreach,FirstCollision_min, ...
    CAMFeasible,CAM_DV_mps,CAMBurn_min,CAMTotalTime_min, ...
    TimeToSafeShell_min,PostCAMMinClearance_m,Status);
R.summary=table(N,mean(PassiveBreach),mean(CAMFeasible), ...
    pct(CAM_DV_mps,95),pct(CAMTotalTime_min,95), ...
    pct(TimeToSafeShell_min,95),nnz(~CAMFeasible), ...
    'VariableNames',{'Trials','PassiveCollisionFraction', ...
    'CandidateEscapeFraction','P95FeasibleCAMDV_mps', ...
    'P95FeasibleCAMResponse_min','P95FeasibleTimeToShell_min', ...
    'UnresolvedCases'});
R.assumptions=cfg;R.events=events;
R.warning=['Screening: conditional artificial failure mixture, ', ...
    'spherical hazard proxy, fixed-direction finite burns, no closed-loop ', ...
    'recovery or independently validated holding orbit.'];
end

function [ok,tExit,minClear]=checkShell(x,n,hazards,screen)
% Must exit the PROTECTED shell within maxExitTime_s and never re-enter
% during the subsequent passive horizon; collision checked separately.
t=unique([0:screen.dt_s:screen.passiveHorizon_s, ...
    screen.passiveHorizon_s]);
clear=nan(size(t));minClear=inf;
for j=1:numel(t)
    xt=hcwPhi(n,t(j))*x;
    clear(j)=FABER_clearance_m(xt(1:3),hazards, ...
        screen.protectedRadius_m);
    minClear=min(minClear,clear(j));
end
ix=find(clear>=0,1,'first');
ok=~isempty(ix);tExit=nan;
if ~ok,return;end
tExit=t(ix);
ok=tExit<=screen.maxExitTime_s && all(clear(ix:end)>=0);
end

function [minClear,bad,first]=driftCheck(x,t0,H,dt,n,hazards,radius)
t=unique([t0:dt:H,H]);minClear=inf;bad=false;first=nan;
for k=1:numel(t)
    xx=hcwPhi(n,t(k))*x;
    c=FABER_clearance_m(xx(1:3),hazards,radius);
    minClear=min(c,minClear);
    if c<0 && ~bad,bad=true;first=t(k);break;end
end
end

function [x,bad]=finiteBurn(x,n,a,T,dt,hazards,radius)
% Constant LVLH acceleration for full firing interval; 4th-order RK.
bad=false;
if T==0,return;end
N=max(1,ceil(T/dt));h=T/N;
for j=1:N
    k1=rhs(x,n,a);k2=rhs(x+.5*h*k1,n,a);
    k3=rhs(x+.5*h*k2,n,a);k4=rhs(x+h*k3,n,a);
    x=x+(h/6)*(k1+2*k2+2*k3+k4);
    if FABER_clearance_m(x(1:3),hazards,radius)<0
        bad=true;return;
    end
end
end
function dx=rhs(x,n,a)
dx=[x(4:6);3*n^2*x(1)+2*n*x(5); ...
    -2*n*x(4);-n^2*x(3)];
dx(4:6)=dx(4:6)+a;
end
function Phi=hcwPhi(n,t)
nt=n*t;c=cos(nt);s=sin(nt);
A=[4-3*c,0,0;6*(s-nt),1,0;0,0,c];
B=[s,2*(1-c),0;2*(c-1),4*s-3*nt,0;0,0,s]/n;
C=[3*n*s,0,0;6*n*(c-1),0,0;0,0,-n*s];
D=[c,2*s,0;-2*s,4*c-3,0;0,0,c];
Phi=[A B;C D];
end
function v=tilt(dv,theta)
if norm(dv)==0,v=dv;return;end
u=dv/norm(dv);r=randn(3,1);r=r-u*(u'*r);
if norm(r)<1e-12,v=dv;return;end
r=r/norm(r);v=norm(dv)*(u*cos(theta)+r*sin(theta));
end
function v=pct(x,p)
x=sort(x(isfinite(x)));
if isempty(x),v=nan;elseif numel(x)==1,v=x;else
    v=interp1(1:numel(x),x,1+(numel(x)-1)*p/100,'linear');
end
end
