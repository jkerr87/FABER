%% FABER_RAAN_Phase_Joint30d_smoketest
% Run once after copying both .m files into the same MATLAB folder.
cfg = struct();
cfg.makeFigures = true;
cfg.exportFiles = false;
R = FABER_RAAN_Phase_Joint30d(cfg);

idxR = find(R.cfg.raanDiff_deg==0,1);
idxP0 = find(R.cfg.phaseDiff_deg==0,1);
idxP20 = find(R.cfg.phaseDiff_deg==20,1);

assert(R.feasible(idxP0,idxR));
assert(R.minDV_mps(idxP0,idxR)==0);
assert(R.feasible(idxP20,idxR));
assert(R.minDV_mps(idxP20,idxR)<R.cfg.DVcap_mps);

f=R.feasible;
assert(all(R.minDV_mps(f)<=R.cfg.DVcap_mps+1e-5));
assert(all(R.bestTOF_days(f)<=R.cfg.Tmax_days+1e-8));
assert(all(abs(R.bestResidualRAAN_deg(f)) <= ...
           R.cfg.maxPoweredRAAN_deg+1e-7));
assert(all(abs(R.bestPhaseResidual_deg(f))<1e-7));
assert(all(abs(R.minDV_mps(f)- ...
    R.bestOrbitChangeDV_mps(f)-R.bestPoweredPlaneDV_mps(f))<1e-6));

fprintf('\nAll FABER joint-screening smoke tests passed.\n');
