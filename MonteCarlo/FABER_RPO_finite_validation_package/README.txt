FABER nine-phase RPO -- finite-thrust validation upgrade

MAIN FILES
  ASV_RPO_MC_config.m
  FABER_RPO_9phase_MC.m
  FABER_RPO_finite_validation.m  (NEW)
  plot_FABER_RPO_finite.m        (NEW)
  FABER_RPO_Pareto.m
  plot_FABER_RPO_P95.m
  main_mc.m
  TECHNICAL_SKELETON.md        (equations, logic, limitations, tests)

QUICK START
  cfg = ASV_RPO_MC_config();
  cfg.N = 30; cfg.makePlots = false; cfg.saveCSV = false;
  cfg.finiteThrustValidation = 'selected';
  R = FABER_RPO_9phase_MC(cfg);
  V = FABER_RPO_finite_validation(R);
  disp(V.table);

SELECT MODES
  'none' / 'selected' / 'pareto' / 'all'
  cfg.finiteMaxCandidates = 100;  % cap for faster debugging
  cfg.finiteMaxCandidates = inf;  % full selected mode / pareto / all

The analytic 5% burn fraction remains a hard pre-screen. A longer burn
could still be physically feasible: relax the 5% gate in a separate
low-thrust design sensitivity case before concluding infeasibility.

FINITE VALIDATION STATUS
  V.attempted : candidates that were integrated
  V.passed    : candidates meeting propagated geometric and terminal checks
  V.failureReason, V.table, V.details : diagnostics and time histories
  V.finitePareto : nondominated only among tested and passed candidates

REQUIRES ORIGINAL analytical_manos.m on MATLAB path.
MATLAB was not available in the build environment, so run checkcode
and a small debug experiment before using results in your report.

Read TECHNICAL_SKELETON.md for the full equation and algorithm review.
