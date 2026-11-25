function V = FABER_RPO_finite_validation(R,mode)
% FABER_RPO_FINITE_VALIDATION
% Post-process analytically feasible nine-phase MC candidates.
% Piecewise ode45 finite thrust/coast for impulsive legs, or ideal HCW
% quintic feedforward for continuous legs. Uses nominal waypoint initial
% conditions separately for each leg. NO closed-loop correction or chained
% off-nominal drift and NO attitude/torque dynamics.
%
% USAGE:
% V = FABER_RPO_finite_validation(R);          % cfg mode
% V = FABER_RPO_finite_validation(R,'pareto'); % override
%
% V.attempted / V.passed are N x M masks. Unattempted is UNKNOWN,
% not infeasible! V.table provides per-candidate results. V.details{i,m}
% provides per-phase time histories, commanded thrust and mass.

cfg = R.cfg;
if nargin < 2 || isempty(mode)
    mode = cfg.finiteThrustValidation;
end

mode = lower(string(mode));
[N,M] = size(R.valid);
V = struct();
V.mode = mode;
V.attempted = false(N,M);
V.passed = false(N,M);
V.details = cell(N,M);
V.phaseReason = cell(N,M);
V.failureReason = strings(N,M);
V.finitePareto = false(N,M);
if mode == "none"
    V.table = table();
    return;
end

assert(ismember(mode,["selected","pareto","all"]), ...
    'Mode must be none, selected, pareto, or all.');

switch mode
    case "all"
        mask = R.valid;
    case "pareto"
        mask = R.valid & R.scenarioPareto;
    case "selected"
        mask = false(N,M);
        % Representative nominal choice for every accepted trial.
        for i=1:N
            m=R.selectedArchitecture(i);
            if isfinite(m) && R.valid(i,m)
                mask(i,m)=true;
            end
        end

        % Ensure minimum DV/TOF, each architecture's P95 DV sample,
        % and maximum phase burn-fraction sample also represented.
        for m=1:M
            ids=find(R.valid(:,m));
            if isempty(ids), continue; end
            [~,a]=min(R.CycleDV_mps(ids,m)); mask(ids(a),m)=true;
            [~,a]=min(R.CycleTOF_min(ids,m)); mask(ids(a),m)=true;
            [~,order]=sort(R.CycleDV_mps(ids,m));
            mask(ids(order(max(1,ceil(.95*numel(order))))),m)=true;
            [~,a]=max(R.MaxPhaseBurnFraction(ids,m)); mask(ids(a),m)=true;
        end
end

%% STRATIFIED FINITE-THRUST CANDIDATE SELECTION
% Spread validation across maneuver architectures.
% Reproducible random ordering within each architecture.

[iiAll,mmAll] = find(mask);

if isfield(cfg,'finiteMaxCandidates') && ...
        isfinite(cfg.finiteMaxCandidates)

    maxCandidates = max(0,floor(cfg.finiteMaxCandidates));

    % Preserve caller's random-number generator state
    oldRng = rng;
    rngCleanup = onCleanup(@() rng(oldRng));
    rng(cfg.seed + 1000,'twister');

    % Randomized trial lists for each architecture
    candidates = cell(M,1);

    for m = 1:M
        trials = find(mask(:,m));

        if ~isempty(trials)
            trials = trials(randperm(numel(trials)));
        end

        candidates{m} = trials;
    end

    nextCandidate = ones(M,1);

    ii = zeros(0,1);
    mm = zeros(0,1);

    % Round-robin selection across architectures
    while numel(ii) < maxCandidates

        added = false;

        for m = 1:M

            k = nextCandidate(m);

            if k <= numel(candidates{m})

                ii(end+1,1) = candidates{m}(k);
                mm(end+1,1) = m;

                nextCandidate(m) = k+1;
                added = true;

                if numel(ii) >= maxCandidates
                    break;
                end
            end
        end

        if ~added
            break;
        end
    end

    fprintf(['\nFinite validation: %d candidates ' ...
        'selected across %d architectures.\n'], ...
        numel(ii),numel(unique(mm)));

else

    % No cap: validate every eligible candidate
    ii = iiAll;
    mm = mmAll;

end

rows = cell(0,11);
for q=1:numel(ii)
    i=ii(q); m=mm(q);
    rec=R.records{i,m};
    if isempty(rec), continue; end
    V.attempted(i,m)=true;
    candidateDetails=cell(2,5);
    phaseReason=repmat("Not run",2,5);
    coupled=isfield(rec,'ReturnTransfer') && ~isempty(rec.ReturnTransfer);
    ok=true;
    maxPos=0; maxVel=0; maxSpeed=0; maxClosing=0;
    totalProp=0; totalBurn=0;
    why="";
    mass=cfg.massTaxi_kg; % continuous propellant mass across BOTH approaches
    for leg=1:2
        legCfg = cfg;
        if leg==1
            approach=rec.Client;
            if isfield(cfg,'clientH_km'), legCfg.h_km=cfg.clientH_km; end
            if isfield(cfg,'clientGeometry'), legCfg.pathGeometry=cfg.clientGeometry; end
        else
            approach=rec.ATLANTE;
            if isfield(cfg,'hostGeometry'), legCfg.pathGeometry=cfg.hostGeometry; end
            if coupled
                % The checked return mass includes payload and the far-field burns.
                mass=rec.massReturnStart_kg;
            elseif isfield(cfg,'enablePayloadMassChange') && cfg.enablePayloadMassChange
                mass=mass+cfg.payloadMass_kg;
            end
        end

        for k=1:5

            z = approach.transferDetails{k};

            % Map five transfer legs to nine-phase indices
            phaseNumber = 2*k-1;

            % Select main or RCS propulsion for this phase
            phaseCfg = FABER_RPO_phase_propulsion( ...
                legCfg,phaseNumber);

            % Propagate with the correct propulsion system
            [P,fail] = validateLeg(z,phaseCfg,mass,k);

            candidateDetails{leg,k}=P;
            if ~isempty(P)
                mass=P.massFinal_kg;
                totalProp=totalProp+P.propellant_kg;
                totalBurn=totalBurn+P.burnTime_s;
                maxPos=max(maxPos,P.positionError_m);
                maxVel=max(maxVel,P.velocityError_mps);
                maxSpeed=max(maxSpeed,P.maxSpeed_mps);
                maxClosing=max(maxClosing,P.maxClosing_mps);
            end

            phaseReason(leg,k)=fail;
            if strlength(fail)>0
                if ok, why=sprintf('leg%d_phase%d_%s',leg,2*k-1,fail); end
                ok=false;
                if ~coupled, break; end
                if isempty(P), mass=mass-approach.fuelPhase_kg(phaseNumber); end
            end
            if coupled && k<5
                % Holds consume their nominal allowance. Only burns are replayed.
                holdPropellant=approach.fuelPhase_kg(phaseNumber+1);
                mass=mass-holdPropellant;
                totalProp=totalProp+holdPropellant;
            end
        end

        if ~ok && ~coupled, break; end
    end

    if totalProp > cfg.maxPropellantCycle_kg
        ok=false; why="finite_cycle_propellant_limit";
    end

    V.passed(i,m)=ok;
    V.failureReason(i,m)=why;
    V.details{i,m}=candidateDetails;
    V.phaseReason{i,m}=phaseReason;
    rows(end+1,:)={i,m,ok,why,R.CycleDV_mps(i,m), ...
        R.CycleTOF_min(i,m),totalProp,totalBurn,maxPos,maxVel, ...
        maxClosing}; %#ok<AGROW>
end

V.table=cell2table(rows,'VariableNames', ...
    {'Trial','Architecture','FiniteFeasible','Reason','IdealDV_mps', ...
    'IdealTOF_min','FinitePropellant_kg','FiniteBurnTime_s', ...
    'MaxPositionError_m','MaxVelocityError_mps','MaxClosingSpeed_mps'});
V.nominalValidCount=nnz(R.valid);
V.attemptedCount=nnz(V.attempted);
V.passedCount=nnz(V.passed);
V.attemptedFraction=V.attemptedCount/max(1,V.nominalValidCount);

% Recompute nondominated set ONLY among those successfully validated.
for i=1:N
    ids=find(V.passed(i,:));
    for m=ids
        dv=R.CycleDV_mps(i,m); tof=R.CycleTOF_min(i,m);
        dominates=R.CycleDV_mps(i,ids)<=dv+1e-10 & ...
            R.CycleTOF_min(i,ids)<=tof+1e-10 & ...
            (R.CycleDV_mps(i,ids)<dv-1e-10 | ...
             R.CycleTOF_min(i,ids)<tof-1e-10);
        V.finitePareto(i,m)=~any(dominates);
    end
end

fprintf('\nFINITE-THRUST VALIDATION: %s\n',mode);
fprintf('Analytically feasible: %d | validated: %d | passed: %d\n', ...
    V.nominalValidCount,V.attemptedCount,V.passedCount);
if ~isempty(V.table), disp(V.table); end
if cfg.saveCSV && ~isempty(V.table)
    if ~exist(cfg.outputFolder,'dir'), mkdir(cfg.outputFolder); end
    writetable(V.table,fullfile(cfg.outputFolder,'finite_thrust_validation.csv'));
end
end

function [P,fail] = validateLeg(z,cfg,mass0,phase)
P=[]; fail="";
T=z.TOF_s;
Isp=cfg.Isp_s; ve=9.80665*Isp;
n=sqrt(cfg.mu_km3s2/(cfg.Re_km+cfg.h_km)^3);
Nout=max(21,round(cfg.finiteOutputSamples));
opts=odeset('RelTol',cfg.finiteRelTol,'AbsTol',cfg.finiteAbsTol);
xinit=[z.startState(:);mass0];

timeline=[]; states=[]; thrust=[];
if z.isImpulsive
    d1=z.burn1(:); d2=z.burn2(:);
    [F1,t1,m1,why]=burnCommand(d1,mass0,cfg);
    if strlength(why)>0, fail=why; return; end
    [F2,t2,~,why]=burnCommand(d2,m1,cfg);
    if strlength(why)>0, fail=why; return; end
    if (t1+t2)>T+1e-7
        fail="burns_exceed_transfer"; return;
    end
    if (t1+t2)/T>cfg.maxBurnFraction+1e-12
        fail="finite_impulsive_fraction"; return;
    end

    segments={ [0,t1],F1;[t1,T-t2],zeros(3,1);[T-t2,T],F2 };
    x=xinit;
    for k=1:3
        t0=segments{k,1}(1); tf=segments{k,1}(2);
        if tf-t0 <= 1e-10, continue; end
        tt=linspace(t0,tf,Nout);
        force=segments{k,2};
        [t,Y]=ode45(@(~,y) dynamics(y,n,force,ve,[],0),tt,x,opts);
        if ~isempty(timeline)
            t=t(2:end); Y=Y(2:end,:);
        end

        timeline=[timeline;t]; %#ok<AGROW>
        states=[states;Y]; %#ok<AGROW>
        thrust=[thrust;repmat(force(:)',numel(t),1)]; %#ok<AGROW>
        x=Y(end,:)';
    end

    totalBurn=t1+t2;
else
    % Exact commanded continuous feedforward u(t) in LVLH.
    tt=linspace(0,T,Nout);
    [timeline,states]=ode45(@(t,y) dynamics(y,n,[],ve, ...
        @(tau) quinticAccel(z.startState,z.targetState,T,tau,n),t), ...
        tt,xinit,opts);
    thrust=zeros(numel(timeline),3);
    for j=1:numel(timeline)
        accel=quinticAccel(z.startState,z.targetState,T,timeline(j),n);
        thrust(j,:)=states(j,7)*accel(:)';
        amag=norm(accel);
        if amag>1e-12
            Favail=directionalThrust(accel/amag,cfg);
            if states(j,7)*amag>Favail+1e-8
                fail="continuous_thrust_exceeded";
            end
        end
    end

    totalBurn=trapz(timeline,double(vecnorm(thrust,2,2)>1e-9));
end

if isempty(states), fail="empty_propagation"; return; end
massFinal=states(end,7);
r=states(:,1:3); v=states(:,4:6);
rho=vecnorm(r,2,2);
rhoDot=sum(r.*v,2)./max(rho,eps);
maxClosing=max([0;-rhoDot]);
maxSpeed=max(vecnorm(v,2,2));
maxAllowed=[cfg.maxSpeed_homing1_mps,cfg.maxSpeed_homing2_mps, ...
    cfg.maxSpeed_closing1_mps,cfg.maxSpeed_closing2_mps, ...
    cfg.maxSpeed_final_mps];
if phase==3 || phase==4
    speedMode=lower(string(cfg.closingSpeedMode));
    if speedMode=="radial", observed=maxClosing;
    else, observed=maxSpeed; end
else
    observed=maxSpeed;
end

if observed>maxAllowed(phase)+1e-7, fail="finite_speed"; end

% Zone safety on FULL finite-thrust path including burn arcs.
tol=max(1e-5,cfg.zoneTolerance_m);
switch phase
    case {1,2}
        if any(rho<cfg.AZ_m-tol), fail="premature_AZ"; end
    case 3
        if any(rho<cfg.KOZ_m-tol), fail="premature_KOZ"; end
        if rho(end)>=cfg.AZ_m, fail="AZ_endpoint"; end
    case 4
        if any(rho>cfg.AZ_m+tol), fail="AZ_reexit"; end
    case 5
        if any(rho>cfg.KOZ_m+tol), fail="KOZ_reexit"; end
end

if phase>=4 && any(rho<cfg.final_m(1)-tol)
    fail="target_passage";
end

if isfield(cfg,'pathGeometry') && ~isempty(fieldnames(cfg.pathGeometry)) && ...
        FABER_path_clearance(r',cfg.pathGeometry)<=0
    fail="finite_obstacle_clearance";
end
if isfield(cfg,'terminalCorridorHalfWidth_m') && isfinite(cfg.terminalCorridorHalfWidth_m)
    inside = rho<=cfg.KOZ_m;
    if any(r(inside,2)>0 | hypot(r(inside,1),r(inside,3))>cfg.terminalCorridorHalfWidth_m)
        fail="finite_terminal_corridor";
    end
end

posErr=norm(states(end,1:3)'-z.targetState(1:3));
velErr=norm(states(end,4:6)'-z.targetState(4:6));
if phase==5
    posTol=cfg.finiteFinalPosTolerance_m;
    velTol=cfg.finiteFinalVelTolerance_mps;
else
    posTol=cfg.finitePosTolerance_m;
    velTol=cfg.finiteVelTolerance_mps;
end

if posErr>posTol, fail="finite_position_error"; end
if velErr>velTol, fail="finite_velocity_error"; end

P=struct('time_s',timeline,'state',states(:,1:6), ...
    'thrust_LVLH_N',thrust,'mass_kg',states(:,7), ...
    'massFinal_kg',massFinal,'propellant_kg',mass0-massFinal, ...
    'burnTime_s',totalBurn,'positionError_m',posErr, ...
    'velocityError_mps',velErr,'maxSpeed_mps',maxSpeed, ...
    'maxClosing_mps',maxClosing,'rho_m',rho,'rhoDot_mps',rhoDot, ...
    'phase',phase,'method',z.method);
end

function dy=dynamics(y,n,F,ve,desiredAccel,t)
r=y(1:3); v=y(4:6); m=y(7);
if isempty(desiredAccel)
    accel=F/max(m,eps);
    Fmag=norm(F);
else
    accel=desiredAccel(t);
    Fmag=m*norm(accel);
end

free=[3*n^2*r(1)+2*n*v(2);-2*n*v(1);-n^2*r(3)];
dy=[v;free+accel;-Fmag/ve];
end

function u=quinticAccel(x0,xf,T,t,n)
% Same quintic path + HCW feedforward as 9-phase generator.
s=t/T; r0=x0(1:3);v0=x0(4:6);rf=xf(1:3);vf=xf(4:6);
D=rf-r0-v0*T; W=(vf-v0)*T;
A=10*D-4*W; B=-15*D+7*W; C=6*D-3*W;
r=r0+v0*t+A*s^3+B*s^4+C*s^5;
v=v0+(3*A*s^2+4*B*s^3+5*C*s^4)/T;
a=(6*A*s+12*B*s^2+20*C*s^3)/T^2;
u=[a(1)-3*n^2*r(1)-2*n*v(2);a(2)+2*n*v(1);a(3)+n^2*r(3)];
end

function [F,tb,mFinal,fail]=burnCommand(dv,m,cfg)
fail=""; F=zeros(3,1); tb=0;mFinal=m;
v=norm(dv);
if v<1e-12, return; end
Fmag=directionalThrust(dv/v,cfg);
if Fmag<=1e-10, fail="no_thrust_in_burn_direction"; return; end
mFinal=m*exp(-v/(cfg.Isp_s*9.80665));
tb=(m-mFinal)*cfg.Isp_s*9.80665/Fmag;
F=Fmag*dv/v;
if tb>cfg.maxBurnDuration_s+1e-8, fail="burn_duration_limit"; end
end

function F=directionalThrust(direction,cfg)
if isempty(cfg.thrusterDirections_LVLH)
    F=cfg.thrust_N*cfg.nActiveThrusters*cfg.thrustEfficiency;
    return;
end

d=direction(:)/norm(direction);
U=cfg.thrusterDirections_LVLH;
U=U./vecnorm(U,2,1);
A=U.*reshape(cfg.thrusterForces_N,1,[]);
Q=eye(3)-d*d';
options=optimoptions('linprog','Display','none');
[w,~,flag]=linprog(-(d'*A)',[],[],Q*A,zeros(3,1), ...
    zeros(size(A,2),1),ones(size(A,2),1),options);
if flag<=0 || isempty(w), F=0;
else, F=max(0,d'*A*w)*cfg.thrustEfficiency; end
end
