function R = FABER_safety_MC_core(events,cfg)
% CONDITIONAL failure-injection, finite CAM and optional recovery screening.
% Not a collision-risk probability, flight controller, or qualified budget.
% All vectors are expressed in the same LVLH frame as the HCW states.
% CAM force vectors must be supplied as PHYSICALLY VERIFIED firing groups.
assert(~isempty(events),'No failure events supplied.');
assert(isfinite(cfg.vehicleMass_kg)&&cfg.vehicleMass_kg>0);
assert(cfg.screen.collisionRadius_m>0 && ...
    cfg.screen.protectedRadius_m>cfg.screen.collisionRadius_m);
assert(cfg.screen.dt_s>0 && cfg.screen.passiveHorizon_s>0);
assert(cfg.detectionDelay_s>=0);
assert(isfield(cfg,'camForceVectorsByVehicle') || ...
    isfield(cfg,'camForceVectors_LVLH_N'), ...
    ['Provide physical CAM firing-group force vectors in LVLH. ', ...
     'Total installed thrust is not sufficient.']);
rng(cfg.seed,'twister');
N=cfg.N;
Trial=(1:N)'; EventIndex=zeros(N,1); Module=zeros(N,1);
Vehicle=strings(N,1); Phase=strings(N,1); Case=strings(N,1);
Mass_kg=nan(N,1); PhysicalClearanceStart_m=nan(N,1);
ContactUnresolved=false(N,1); PassiveMinClearance_m=nan(N,1);
PassiveBreach=false(N,1); FirstCollision_min=nan(N,1);
CAMRequired=false(N,1); CAMFeasible=false(N,1);

CAM_DV_mps=nan(N,1); CAMBurn_min=nan(N,1);
TimeToSafeShell_min=nan(N,1); CAMTotalTime_min=nan(N,1);
PostCAMMinClearance_m=nan(N,1);
RecoveryFeasible=false(N,1); RecoveryDV_mps=nan(N,1);
RecoveryTOF_min=nan(N,1); TotalContingencyDV_mps=nan(N,1);
Status=strings(N,1);
isBurn=arrayfun(@(e)all(isfinite(e.dvNom(:))),events);

% Inject a conditional fault at one recorded burn or coast event per trial.
for k=1:N
    pool=find(isBurn==mod(k,2));
    if isempty(pool),pool=1:numel(events);end
    ei=pool(randi(numel(pool))); e=events(ei);
    EventIndex(k)=ei; Vehicle(k)=string(e.vehicle);
    Phase(k)=string(e.phase); Module(k)=e.index;
    hazards=eventHazards(cfg,e);
    mass=cfg.vehicleMass_kg;
    if isfield(e,'mass_kg') && isfinite(e.mass_kg) && e.mass_kg>0
        mass=e.mass_kg;
    end

    Mass_kg(k)=mass;
    x=e.xPre(:); dv=e.dvNom(:);
    % Check baseline geometry BEFORE injecting navigation uncertainty.
    c0=FABER_clearance_m(x(1:3),hazards,cfg.screen.collisionRadius_m);
    PhysicalClearanceStart_m(k)=c0;
    if c0<=0 || (isfield(e,'contactUnresolved') && e.contactUnresolved)
        ContactUnresolved(k)=true;
        Status(k)="Unresolved contact geometry (not screened)";
        continue;
    end

    x(1:3)=x(1:3)+cfg.sigmaPos_m*randn(3,1);
    x(4:6)=x(4:6)+cfg.sigmaVel_mps*randn(3,1);
    if all(isfinite(dv))
        z=rand();
        if z<cfg.caseWeights(1)
            Case(k)="Missed burn"; actual=zeros(3,1);
        elseif z<sum(cfg.caseWeights(1:2))
            Case(k)="Thrust magnitude error";
            actual=dv*(1+cfg.sigmaThrustFrac*randn());
        else
            Case(k)="Pointing error";
            actual=tilt(dv,cfg.sigmaPointing_rad*randn());
        end

        x(4:6)=x(4:6)+actual;
    else
        Case(k)="Control cutoff / coast drift";
    end

    [pc,breach,tb]=driftCheck(x,cfg.screen.passiveHorizon_s, ...
        cfg.screen.dt_s,cfg.n,hazards,cfg.screen.collisionRadius_m);
    PassiveMinClearance_m(k)=pc;PassiveBreach(k)=breach;
    FirstCollision_min(k)=tb/60;
    % A drift that remains outside the protected shell, collision-free
    % throughout the horizon, requires NO CAM under this study policy.
    [passiveSafe,~,~]=checkShell(x,cfg.n,hazards,cfg.screen,0);
    CAMRequired(k)=breach || ~passiveSafe;
    if ~CAMRequired(k)
        CAMFeasible(k)=true;CAM_DV_mps(k)=0;CAMBurn_min(k)=0;
        CAMTotalTime_min(k)=0;TimeToSafeShell_min(k)=0;
        Status(k)="No CAM required: passive safe";
        % Already on a passive safe trajectory; recovery still optional.
        if cfg.enableActiveCAM && cfg.recovery.enabled
            forceVectors=eventForces(cfg,e);
            [RecoveryFeasible(k),RecoveryDV_mps(k),RecoveryTOF_min(k)] = ...
                recoveryCheck(x,mass,cfg,hazards,forceVectors);
        end

        if RecoveryFeasible(k),TotalContingencyDV_mps(k)=RecoveryDV_mps(k);end
        continue;
    end

    [~,badDelay,~]=driftCheck(x,cfg.detectionDelay_s, ...
        cfg.screen.burnStep_s,cfg.n,hazards,cfg.screen.collisionRadius_m);
    if badDelay
        Status(k)="Collision before detection/response";continue;
    end

    if ~cfg.enableActiveCAM
        Status(k)="CAM unassessed: no verified firing-group vectors";
        continue;
    end

    xd=hcwPhi(cfg.n,cfg.detectionDelay_s)*x;
    forceVectors=eventForces(cfg,e);
    bestDV=inf;bestState=[];bestTime=nan;
    for j=1:size(forceVectors,2)
        Fvec=forceVectors(:,j);F=norm(Fvec);
        if F<=0,continue;end
        dir=Fvec/F;
        for mag=cfg.camDVgrid_mps(:)'
            if mag>cfg.camMaxDV_mps,continue;end
            % Constant-mass Δv approximation; exhaust depletion omitted.
            Tburn=mag*mass/F;
            if Tburn>cfg.camMaxBurn_s,continue;end
            [xb,badBurn]=finiteBurn(xd,cfg.n,Fvec/mass,Tburn, ...
                cfg.screen.burnStep_s,hazards,cfg.screen.collisionRadius_m);
            if badBurn,continue;end
            [~,badDrift,~]=driftCheck(xb,cfg.screen.passiveHorizon_s, ...
                cfg.screen.dt_s,cfg.n,hazards,cfg.screen.collisionRadius_m);
            if badDrift,continue;end
            [safe,tExit,cMin]=checkShell(xb,cfg.n,hazards,cfg.screen, ...
                cfg.screen.maxExitTime_s);
            if safe && mag<bestDV
                bestDV=mag;bestState=xb;bestTime=Tburn;
                CAMBurn_min(k)=Tburn/60;
                TimeToSafeShell_min(k)=(cfg.detectionDelay_s+Tburn+tExit)/60;
                CAMTotalTime_min(k)=(cfg.detectionDelay_s+Tburn)/60;
                PostCAMMinClearance_m(k)=cMin;
            end
        end
    end

    if isfinite(bestDV)
        CAMFeasible(k)=true;CAM_DV_mps(k)=bestDV;
        Status(k)="Escape found; recovery not established";
        % Recovery starts after the passive safety horizon; departure state
        % is the freely drifted post-CAM state, not the nominal endpoint.
        xRecover=hcwPhi(cfg.n,cfg.screen.passiveHorizon_s)*bestState;
        [RecoveryFeasible(k),RecoveryDV_mps(k),RecoveryTOF_min(k)] = ...
            recoveryCheck(xRecover,mass,cfg,hazards,forceVectors);
        if RecoveryFeasible(k)
            TotalContingencyDV_mps(k)=bestDV+RecoveryDV_mps(k);
            Status(k)="Escape and preliminary recovery passed";
        end
    else
        Status(k)="No feasible candidate escape";
    end
end

R.trials=table(Trial,EventIndex,Module,Vehicle,Phase,Case,Mass_kg, ...
    PhysicalClearanceStart_m,ContactUnresolved,PassiveMinClearance_m, ...
    PassiveBreach,FirstCollision_min,CAMRequired,CAMFeasible,CAM_DV_mps, ...
    CAMBurn_min,CAMTotalTime_min,TimeToSafeShell_min, ...
    PostCAMMinClearance_m,RecoveryFeasible,RecoveryDV_mps, ...
    RecoveryTOF_min,TotalContingencyDV_mps,Status);

% Report unresolved contact geometry separately from screened fault outcomes.
valid=~ContactUnresolved;
required=valid & CAMRequired;
R.summary=table(N,nnz(valid),nnz(ContactUnresolved),nnz(required), ...
    nnz(required & CAMFeasible),nnz(required & ~CAMFeasible), ...
    pct(CAM_DV_mps(required & CAMFeasible),95), ...
    nnz(valid & RecoveryFeasible), ...
    pct(TotalContingencyDV_mps(valid & RecoveryFeasible),95), ...
    'VariableNames',{'Trials','ScreenedTrials','UnresolvedContactTrials', ...
    'CAMRequiredTrials','CAMFoundTrials','CAMUnresolvedTrials', ...
    'P95CAMDV_Conditional_mps','RecoveryFoundTrials', ...
    'P95FullContingencyDV_Conditional_mps'});
R.assumptions=cfg;R.events=events;
R.warning=['Conditional fault scenarios; no actual fault probabilities. ', ...
    'CAM force groups require verified orientation and torque allocation. ', ...
    'Near-contact geometry unresolved; recovery is an ideal HCW two-impulse ', ...
    'screen, not a flight recovery plan. P95 values exclude unresolved cases.'];
end

function hazards=eventHazards(cfg,e)
hazards=cfg.hazards;
if isfield(e,'hazardsCenters_m') && ~isempty(e.hazardsCenters_m)
    hazards.centers_m=e.hazardsCenters_m;
    hazards.radii_m=zeros(size(e.hazardsCenters_m,1),1);
end

if isfield(e,'hazardsSegmentsA_m')
    hazards.segmentsA_m=e.hazardsSegmentsA_m;
    hazards.segmentsB_m=e.hazardsSegmentsB_m;
    hazards.segmentRadii_m=e.hazardsSegmentRadii_m;
end
end

function F=eventForces(cfg,e)
% Explicit validated net force VECTORS of allowed simultaneous firing groups.
% Zero lateral / torque feasibility must be certified upstream.
if isfield(cfg,'camForceVectorsByVehicle') && ...
        isfield(cfg.camForceVectorsByVehicle,char(e.vehicle))
    F=cfg.camForceVectorsByVehicle.(char(e.vehicle));
else
    assert(isfield(cfg,'camForceVectors_LVLH_N'), ...
        'No CAM force-group vectors for vehicle %s',char(e.vehicle));
    F=cfg.camForceVectors_LVLH_N;
end

assert(size(F,1)==3 && ~isempty(F) && all(isfinite(F(:))) && ...
    all(vecnorm(F,2,1)>0), ...
    'Supply nonzero [3 x Ngroup] verified force vectors in LVLH.');
end

function [ok,dv,tofMin]=recoveryCheck(x,m,cfg,hazards,forceVectors)
% Separate first-pass recovery after escape/24h drift.
% Computes two HCW burns to an EXPLICIT safe 6-state; requires each burn
% direction to match a listed, hardware-verified constant-force group.
% Integrates finite burn/coast/burn, checking geometry and endpoint errors.
% NOT a guidance law, re-docking plan, or flight-qualified recovery.
ok=false;dv=nan;tofMin=nan;
if ~isfield(cfg,'recovery') || ~cfg.recovery.enabled,return;end
rf=cfg.recovery.targetState(:);T=cfg.recovery.TOF_s;n=cfg.n;
assert(numel(rf)==6 && T>0);
P=hcwPhi(n,T);B=P(1:3,4:6);
if rcond(B)<1e-9,return;end
vPlus=B\(rf(1:3)-P(1:3,1:3)*x(1:3));
dv1=vPlus-x(4:6);
vEnd=P(4:6,1:3)*x(1:3)+P(4:6,4:6)*vPlus;
dv2=rf(4:6)-vEnd;
dv=norm(dv1)+norm(dv2);tofMin=T/60;
if dv>cfg.recovery.maxDV_mps,return;end
[F1,t1,ok1]=selectRecoveryForce(dv1,m,forceVectors,cfg.recovery);
[F2,t2,ok2]=selectRecoveryForce(dv2,m,forceVectors,cfg.recovery);
if ~ok1 || ~ok2 || t1+t2>T,return;end
% Finite propagation from actual state. Candidate burn/arrival timing fixed.
[y,bad]=finiteBurn(x,n,F1/m,t1,cfg.screen.burnStep_s, ...
    hazards,cfg.screen.collisionRadius_m);
if bad,return;end
[y,bad]=coastCheck(y,n,T-t1-t2, ...
    cfg.screen.dt_s,hazards,cfg.screen.collisionRadius_m);
if bad,return;end
[y,bad]=finiteBurn(y,n,F2/m,t2,cfg.screen.burnStep_s, ...
    hazards,cfg.screen.collisionRadius_m);
if bad,return;end
if norm(y(1:3)-rf(1:3))>cfg.recovery.posTolerance_m || ...
        norm(y(4:6)-rf(4:6))>cfg.recovery.velTolerance_mps
    return;
end

if cfg.recovery.stayOutsideProtected
    % Check whole ideal impulsive transfer and terminal state for protected
    % geometry. Finite departure/arrival samples also need separate checks.
    xIdeal=[x(1:3);vPlus];
    for t=linspace(0,T,max(301,ceil(T/cfg.screen.dt_s)+1))
        xt=hcwPhi(n,t)*xIdeal;
        if FABER_clearance_m(xt(1:3),hazards, ...
                cfg.screen.protectedRadius_m)<0,return;end
    end
end

ok=true;
end

function [F,t,ok]=selectRecoveryForce(dv,m,forceVectors,rec)
F=zeros(3,1);t=0;ok=true;
mag=norm(dv);
if mag<1e-12,return;end
unit=dv/mag;
Fmag=vecnorm(forceVectors,2,1);
cosines=(unit'*forceVectors)./Fmag;
[best,j]=max(cosines);
if best<cosd(rec.directionTolerance_deg),ok=false;return;end
F=forceVectors(:,j);t=mag*m/Fmag(j);
if t>rec.maxBurn_s,ok=false;end
end

function [ok,tExit,minClear]=checkShell(x,n,hazards,screen,maxExit)
t=unique([0:screen.dt_s:screen.passiveHorizon_s,screen.passiveHorizon_s]);
clear=zeros(size(t));minClear=inf;
for j=1:numel(t)
    xt=hcwPhi(n,t(j))*x;
    clear(j)=FABER_clearance_m(xt(1:3),hazards,screen.protectedRadius_m);
    minClear=min(minClear,clear(j));
end

ix=find(clear>=0,1,'first');ok=~isempty(ix);tExit=nan;
if ~ok,return;end
tExit=t(ix);
ok=tExit<=maxExit && all(clear(ix:end)>=0);
end

function [minClear,bad,first]=driftCheck(x,H,dt,n,hazards,radius)
t=unique([0:dt:H,H]);minClear=inf;bad=false;first=nan;
for k=1:numel(t)
    xx=hcwPhi(n,t(k))*x;
    c=FABER_clearance_m(xx(1:3),hazards,radius);
    minClear=min(c,minClear);
    if c<0 && ~bad,bad=true;first=t(k);break;end
end
end

function [x,bad]=coastCheck(x,n,T,dt,hazards,radius)
bad=false;x0=x;
for tau=unique([0:dt:T,T])
    xt=hcwPhi(n,tau)*x0;
    if FABER_clearance_m(xt(1:3),hazards,radius)<0
        bad=true;return;
    end
end

x=hcwPhi(n,T)*x0;
end

function [x,bad]=finiteBurn(x,n,a,T,dt,hazards,radius)
bad=false;
if T==0,return;end
N=max(1,ceil(T/dt));h=T/N;
if FABER_clearance_m(x(1:3),hazards,radius)<0,bad=true;return;end
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
dx=[x(4:6);3*n^2*x(1)+2*n*x(5);-2*n*x(4);-n^2*x(3)];
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
