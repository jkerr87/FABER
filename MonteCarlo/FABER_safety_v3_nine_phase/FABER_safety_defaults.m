function cfg = FABER_safety_defaults(cfg,mu,Re,h,vehicle)
% Preliminary defaults. Actual CAM firing groups, collision geometry,
% masses, detection times and recovery targets MUST be mission-verified.
if nargin<1 || isempty(cfg),cfg=struct();end
D=struct();
D.seed=4242;D.N=100;D.n=sqrt(mu/(Re+h)^3);
D.maxSampleStep_s=600;
D.sigmaPos_m=0.5;D.sigmaVel_mps=0.002;
D.sigmaThrustFrac=0.03;D.sigmaPointing_rad=deg2rad(1);
D.caseWeights=[1/3 1/3 1/3];
D.detectionDelay_s=60;
D.screen=struct('passiveHorizon_s',24*3600,'dt_s',300, ...
    'protectedRadius_m',30,'collisionRadius_m',5, ...
    'burnStep_s',2,'maxExitTime_s',3600);
D.camDVgrid_mps=[0.002 0.005 0.01 0.02 0.03 0.05 ...
    0.075 0.1 0.15 0.2 0.5 1];
D.camMaxDV_mps=1;D.camMaxBurn_s=300;
D.vehicleMass_kg=NaN;
D.hazards=struct('centers_m',[0 0 0]);
% Physical net force of independently executable simultaneous firing groups,
% in the relevant LVLH frame: cfg.camForceVectors_LVLH_N=[fx;fy;fz] columns.
% No directions are assumed or fabricated by default.
D.camForceVectors_LVLH_N=[];
% Recovery is explicitly optional and separate from CAM. No target inferred.
D.recovery=struct('enabled',false,'targetState',[], ...
    'TOF_s',3600,'maxDV_mps',2,'maxBurn_s',300, ...
    'posTolerance_m',5,'velTolerance_mps',0.02, ...
    'directionTolerance_deg',5,'stayOutsideProtected',true);
if strcmpi(vehicle,'assembly')
    D.screen.protectedRadius_m=40;
end
f=fieldnames(D);
for k=1:numel(f)
    if ~isfield(cfg,f{k}),cfg.(f{k})=D.(f{k});end
end
sf=fieldnames(D.screen);
for k=1:numel(sf)
    if ~isfield(cfg.screen,sf{k})
        cfg.screen.(sf{k})=D.screen.(sf{k});
    end
end
rf=fieldnames(D.recovery);
for k=1:numel(rf)
    if ~isfield(cfg.recovery,rf{k})
        cfg.recovery.(rf{k})=D.recovery.(rf{k});
    end
end
assert(isfinite(cfg.vehicleMass_kg) && cfg.vehicleMass_kg>0, ...
    'Set actual reference safety vehicle mass: cfg.vehicleMass_kg');
assert(all(cfg.caseWeights>=0) && ...
    abs(sum(cfg.caseWeights)-1)<1e-12, ...
    'Failure case weights must sum to one.');
if cfg.recovery.enabled
    assert(numel(cfg.recovery.targetState)==6 && ...
        all(isfinite(cfg.recovery.targetState)), ...
        'Set a verified 6-state safe recovery target.');
end
end
