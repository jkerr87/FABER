# Sunlight, eclipse and RPO hold trade

The study runs the same Phase A geometry in winter (1 February 2038) and
summer (21 June 2038). Each season compares zero, one and two added
holding opportunities during Taxi Homing 1. Up to five retrievals can
be specified.

```matlab
addpath('Environment/SunlightEclipseRPO');
trade = run_FABER_sunlight_eclipse_study([], ...
    fullfile(pwd,'SunlightEclipseResults'), 5, true);
```

Pass an existing Phase A budget as the first argument to keep its
orbital geometry and selected trajectories. Pass `false` as the fourth
argument for CSVs without the individual MATLAB figures.

- Active RPO is sunlight gated. This includes added H1 sublegs and
  reacquisition after a hold. The terminal closing chain, HP4 fly-through
  and the applicable sunlit service still reserve a single sunlight window.
- Far-field Lambert coast can span eclipse. The far-field table logs
  modeled coast eclipse and the instantaneous lighting at each burn
  epoch. The full-sun column audits the actual timed leg and does not
  represent a reoptimized sunlight-only far-field schedule. Optical
  navigation or full finite-duration burn requirements must be assessed
  before authorizing eclipse burns.
- An added gate is stationary on the negative V-bar outside the approach
  zone. Homing 1 is rebuilt as shorter impulses and independently
  screened for speed, clearance, braking distance and burn duration.
  The existing HCW drift and reacquisition screen applies to eclipse waits.
  An unsuccessful preview is discarded. The scheduler can try fewer
  gates or fall back to the original route. Added gates are experimental,
  and this trade does not replace finite-thrust or closed-loop validation.
- The trade does not insert a passive hold into final approach or a
  fly-through. If the terminal sequence cannot fit, the existing shorter
  RPO candidate search and park-departure retiming remain applicable.

`Sunlight_Eclipse_Comparison.csv` summarizes scheduled waits, delivery,
assembly and the strict sunlight audit. `Far_Field_Lighting.csv` lists
coast and burn-epoch lighting for each executed transfer. Each case
folder contains the original event, hold, lighting and EPS exports.
