function cfg = FABER_safety_defaults(cfg,mu,Re,h,vehicle)
% Placeholder screening assumptions MUST be replaced by FABER requirements.
D=struct();
D.seed=4242; D.N=100;D.n=sqrt(mu/(Re+h)^3);
D.maxSampleStep_s=600; % <=10 min for assembly, selectable
D.sigmaPos_m=0.5; D.sigmaVel_mps=0.002;
D.sigmaThrustFrac=0.03;D.sigmaPointing_rad=deg2rad(1);
D.caseWeights=[1/3 1/3 1/3]; % CONDITIONAL stress-test mixture, not event probabilities
D.detectionDelay_s=60;
D.screen=struct('passiveHorizon_s',24*3600,'dt_s',300, ...
    'protectedRadius_m',30,'collisionRadius_m',5);
D.camDVgrid_mps = [ ...
    0.002 0.005 0.01 0.02 0.03 ...
    0.05 0.075 0.1 0.15 0.2 0.5 1];
D.camMaxDV_mps=1;
D.camMaxBurn_s=300;
D.screen.burnStep_s=2;
D.screen.maxExitTime_s=3600;
D.camForce_N=NaN;
D.vehicleMass_kg=NaN;
D.hazards=struct('centers_m',[0 0 0]);
if strcmpi(vehicle,'assembly')
    D.screen.protectedRadius_m=40; % test value, NOT a verified hold requirement
end
f=fieldnames(D);
for k=1:numel(f)
    if ~isfield(cfg,f{k}),cfg.(f{k})=D.(f{k});end
end
sf=fieldnames(D.screen);
for kk=1:numel(sf)
    if ~isfield(cfg.screen,sf{kk})
        cfg.screen.(sf{kk})=D.screen.(sf{kk});
    end
end
assert(isfinite(cfg.camForce_N) && cfg.camForce_N>0,...
    'Set actual available emergency CAM thrust: cfg.camForce_N');
assert(isfinite(cfg.vehicleMass_kg) && cfg.vehicleMass_kg>0,...
    'Set actual loaded/unloaded mass: cfg.vehicleMass_kg');
assert(abs(sum(cfg.caseWeights)-1)<1e-12,'caseWeights must sum to 1');
vec=[eye(3),-eye(3)];
for i=1:3
  for j=i+1:3
      v=zeros(3,1);v([i j])=1;vec=[vec,v/norm(v),-v/norm(v)]; %#ok<AGROW>
      v([i j])=[1;-1];vec=[vec,v/norm(v),-v/norm(v)]; %#ok<AGROW>
  end
end
cfg.camDirections=vec;
end
