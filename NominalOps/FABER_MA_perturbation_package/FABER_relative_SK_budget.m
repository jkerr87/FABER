function S=FABER_relative_SK_budget(opts)
% FABER_RELATIVE_SK_BUDGET  Open-loop, discrete relative-orbit maintenance.
% Phase-A Mission Analysis: no GNC controller or navigation estimator.
%
% Chief (host) and deputy (Taxi/ATLANTE) independently propagated in ECI
% under central gravity + optional J2, drag and SRP. Compare to a nominal
% bounded HCW holding ellipse in the chief's instantaneous LVLH frame.
% At scheduled inspections, IF the error exceeds a deadband, solve a
% TWO-IMPULSE recovery to the reference ellipse at t+correctionTOF_s.
% Record DV, thrust/fuel and recovery feasibility. During correction the
% chief/deputy are propagated with the same environmental dynamics.
%
% Example:
% opts=struct('h_km',680.12,'inclination_deg',98.12,...
%     'duration_days',7,'centerAlongTrack_m',-1000,...
%     'radialAmplitude_m',20,'checkInterval_s',5400);
% S=FABER_relative_SK_budget(opts);
%
% No feedback: scheduled geometric correction, ideal burns, no estimation.
% 'correctionTOF_s' is executed DURING the idle interval; no additional
% elapsed mission time is charged unless an operational conflict occurs.
if nargin<1,opts=struct();end
opts=dflt(opts,'env',FABER_MA_dynamics('defaults'));
opts=dflt(opts,'h_km',680.12);
opts=dflt(opts,'inclination_deg',98.12);
opts=dflt(opts,'raan_deg',0);
opts=dflt(opts,'argumentLatitude_deg',0);
opts=dflt(opts,'duration_days',7);

opts=dflt(opts,'architecture','ellipse'); % 'ellipse' or 'vbar'
opts=dflt(opts,'centerAlongTrack_m',-1000);
opts=dflt(opts,'radialAmplitude_m',20);
opts=dflt(opts,'checkInterval_s',5400);
opts=dflt(opts,'correctionTOF_s',1800);
opts=dflt(opts,'positionTolerance_m',10);

opts=dflt(opts,'velocityTolerance_mps',0.01);
opts=dflt(opts,'minClearance_m',300);
opts=dflt(opts,'thrust_N',25);
opts=dflt(opts,'Isp_s',220); % placeholder; set actual Isp
opts=dflt(opts,'maxBurnFraction',0.05);
opts=dflt(opts,'maxStep_s',300);
assert(opts.env.deputy.mass_kg>0 && opts.env.chief.mass_kg>0);
assert(opts.duration_days>0 && opts.checkInterval_s>0 && opts.correctionTOF_s>0);
assert(opts.thrust_N>0 && opts.Isp_s>0);
e=opts.env; n=sqrt(e.mu_m3s2/(e.Re_m+opts.h_km*1000)^3);
chief=FABER_MA_dynamics('chiefInitial',opts.h_km,...
    opts.inclination_deg,opts.raan_deg,opts.argumentLatitude_deg,e);
X0=reference(0,n,opts);
dep=FABER_MA_dynamics('deputyECI',chief,X0);
Y=[chief;dep];
t=0; tf=opts.duration_days*86400;
maxErr=0; maxVelErr=0; minClear=inf; maxRange=0;
dvtot=0; propTot=0; burnTime=0; eventRows=cell(0,10);
corrected=0; missed=0; initialMass=e.deputy.mass_kg;
logTime=[]; logRange=[]; logError=[];
while t<tf-1e-6
    dt=min(opts.checkInterval_s,tf-t);
    [tt,ZZ]=propagate(Y,dt,e,opts.maxStep_s);
    for a=1:numel(tt)
        X=FABER_MA_dynamics('relative',ZZ(a,1:6)',ZZ(a,7:12)');
        target=reference(t+tt(a),n,opts);
        err=X-target;
        maxErr=max(maxErr,norm(err(1:3)));
        maxVelErr=max(maxVelErr,norm(err(4:6)));
        r=norm(X(1:3));
        minClear=min(minClear,r); maxRange=max(maxRange,r);
        logTime(end+1,1)=t+tt(a); %#ok<AGROW>
        logRange(end+1,1)=r; %#ok<AGROW>
        logError(end+1,1)=norm(err(1:3)); %#ok<AGROW>
    end

    t=t+dt; Y=ZZ(end,:)';
    X=FABER_MA_dynamics('relative',Y(1:6),Y(7:12));
    targetNow=reference(t,n,opts);
    errNow=X-targetNow;
    trigger=norm(errNow(1:3))>opts.positionTolerance_m || ...
            norm(errNow(4:6))>opts.velocityTolerance_mps;
    if ~trigger,continue;end
    if tf-t<opts.correctionTOF_s
        missed=missed+1;
        eventRows(end+1,:)={t/86400,norm(errNow(1:3)), ...
            norm(errNow(4:6)),NaN,NaN,NaN,NaN,NaN,NaN,"no_time"}; %#ok<AGROW>
        continue;
    end

    Tc=opts.correctionTOF_s;
    P=FABER_MA_dynamics('phi',n,Tc);
    if rcond(P(1:3,4:6))<1e-10
        error('Correction TOF singular or ill-conditioned; change correctionTOF_s');
    end

    desired=reference(t+Tc,n,opts);
    vPost=P(1:3,4:6)\(desired(1:3)-P(1:3,1:3)*X(1:3));
    dv1=vPost-X(4:6);
    [Q,~]=FABER_MA_dynamics('frame',Y(1:6));
    Y(10:12)=Y(10:12)+Q*dv1;
    [tt,ZZ]=propagate(Y,Tc,e,opts.maxStep_s);
    for a=1:numel(tt)
        Xi=FABER_MA_dynamics('relative',ZZ(a,1:6)',ZZ(a,7:12)');
        tar=reference(t+tt(a),n,opts);
        maxErr=max(maxErr,norm(Xi(1:3)-tar(1:3)));
        r=norm(Xi(1:3)); minClear=min(minClear,r);maxRange=max(maxRange,r);
        logTime(end+1,1)=t+tt(a); %#ok<AGROW>
        logRange(end+1,1)=r; %#ok<AGROW>
        logError(end+1,1)=norm(Xi(1:3)-tar(1:3)); %#ok<AGROW>
    end

    t=t+Tc; Y=ZZ(end,:)';
    Xend=FABER_MA_dynamics('relative',Y(1:6),Y(7:12));
    desire=reference(t,n,opts);
    dv2=desire(4:6)-Xend(4:6);
    [Q,~]=FABER_MA_dynamics('frame',Y(1:6));
    Y(10:12)=Y(10:12)+Q*dv2;
    dv=norm(dv1)+norm(dv2);
    massBefore=e.deputy.mass_kg;
    massAfter=massBefore*exp(-dv/(opts.Isp_s*9.80665));
    dm=massBefore-massAfter;
    bt=(dm*opts.Isp_s*9.80665)/opts.thrust_N;
    e.deputy.mass_kg=massAfter;
    posResidual=norm(Xend(1:3)-desire(1:3));
    status="ok";
    if bt/Tc>opts.maxBurnFraction,status="burn_fraction";end
    if posResidual>opts.positionTolerance_m,status="residual_position";end
    if minClear<opts.minClearance_m,status="clearance";end
    dvtot=dvtot+dv;propTot=propTot+dm;burnTime=burnTime+bt;
    corrected=corrected+1;
    eventRows(end+1,:)={t/86400,norm(errNow(1:3)), ...
        norm(errNow(4:6)),norm(dv1),norm(dv2),dv,dm,bt,posResidual,status}; %#ok<AGROW>
end

S=struct();S.options=opts;
S.corrections=cell2table(eventRows,'VariableNames',{'Epoch_days','PreCorrectionPosError_m',...
    'PreCorrectionVelError_mps','Burn1_mps','Burn2_mps','DV_mps',...
    'Propellant_kg','BurnTime_s','TerminalPositionError_m','Status'});
S.Ncorrections=corrected;S.NinsufficientTime=missed;
S.TotalDV_mps=dvtot;S.Propellant_kg=propTot;
S.BurnTime_s=burnTime;S.DVperYear_mps=dvtot/opts.duration_days*365.25;
S.MaxPositionError_m=maxErr;S.MaxVelocityError_mps=maxVelErr;
S.MinSeparation_m=minClear;S.MaxSeparation_m=maxRange;

S.clearanceFeasible=minClear>=opts.minClearance_m;
S.correctionFeasible=all(string(S.corrections.Status)=="ok") && S.clearanceFeasible && missed==0;
S.time_s=logTime;S.separation_m=logRange;S.positionError_m=logError;
S.initialMass_kg=initialMass;S.finalMass_kg=e.deputy.mass_kg;
S.MissionTimePenalty_s=0; % correction executed within allocated idle interval
fprintf('Relative SK %s: %.5f m/s in %.2f d, %d corrections, min rho %.2f m\n',...
    opts.architecture,dvtot,opts.duration_days,corrected,minClear);
if ~S.correctionFeasible
    warning('Some corrections failed preliminary timing/clearance/tolerance checks.');
end
end

function X=reference(t,n,o)
switch lower(char(o.architecture))
    case 'ellipse'
        A=o.radialAmplitude_m; yc=o.centerAlongTrack_m;
        X=[A*cos(n*t);yc-2*A*sin(n*t);0; ...
            -A*n*sin(n*t);-2*A*n*cos(n*t);0];
    case 'vbar'
        X=[0;o.centerAlongTrack_m;0;0;0;0];
    otherwise
        error('Unsupported architecture %s',o.architecture);
end
end

function [tt,Z]=propagate(Y,dt,e,maxStep)
settings=odeset('RelTol',1e-9,'AbsTol',1e-7,'MaxStep',maxStep);
[tt,Z]=ode45(@(t,z)rhs(z,e),[0 dt],Y,settings);
end

function d=rhs(z,e)
a1=FABER_MA_dynamics('accelerations',z(1:6),e,'chief');
a2=FABER_MA_dynamics('accelerations',z(7:12),e,'deputy');
d=[z(4:6);a1;z(10:12);a2];
end

function s=dflt(s,k,val)
if ~isfield(s,k)||isempty(s.(k)),s.(k)=val;end
end
