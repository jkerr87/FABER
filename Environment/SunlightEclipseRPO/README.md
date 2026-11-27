# Sunlight, eclipse and RPO hold trade

## Environmental schedule and delta-v sweep

```matlab
addpath('Environment/SunlightEclipseRPO');
study = MA_main('MA_Duty_Cycle',budget);
```

The optional third argument overrides the study configuration. The default
generates 120 Monte Carlo trials, schedules up to 12 feasible alternate paths
plus the Phase A baseline, and tests departure delays from 0 to 180 minutes
in 15-minute steps. Rows in `cfg.sweepHoldStrategies` give the extra Taxi
homing holds and extra assembly holds, respectively. Use `cfg.cluster="SSO"`
with the corresponding Phase A budget for an SSO run.

The trial CSV records passes and first blockers. The convergence CSV shows
acceptance and best and P95 modeled delta-v on nested subsets of the grid.
The full EPS and sampled sunlight replay applies to the selected schedule.
Vehicle totals reconcile known raw Phase A values by budget unit when the
selected replay passes. Phase A margins and open items remain visible.

## Winter Kiki duty cycle

Run the two-delivery winter case and open the tabbed mission, assembly,
and per-installation duty cycles in one MATLAB window.

```matlab
addpath('Environment/SunlightEclipseRPO');
[study,budget] = main_FABER_kiki_winter('Kiki_Winter_30m');
study45 = main_FABER_kiki_winter('Kiki_Winter_45m',1.5,budget);
```

Pass the folder containing the existing `Winter_On_Demand_Five` and
`Kiki_Winter_Duty_Cycles` exports to show all five deliveries and replay four
alternating assembly round trips per module in the light-mode duty-cycle window.

```matlab
study = main_FABER_kiki_winter('Winter_Delivery_Preview',[],[],resultsFolder);
```

Pass `50` as the fifth argument for the full planning case. The mission tab
covers the entire replay. Each module tab shows all its trips and their
recharge windows, with a separate close-up of its first host recharge when available.
The replay uses the exported winter activity durations, including
holds, and does not recompute sunlight or Atlante battery reserve.
Each Kiki returns to Atlante after every round trip. The replay schedules a
Kiki top-up at berth when its next trip would reduce its battery below reserve.

Use the optional sixth argument to vary the assembly duty cycle.

```matlab
options = struct('respectSunlight',false,'leaveRangoOnHost',true, ...
    'rangoAssembly_min',32.17,'rangoChargeConnect_min',5);
study = main_FABER_kiki_winter('Assembly_Host_Replay',[],[],resultsFolder,3,options);
```

The sunlight case preserves the template waits. Set `respectSunlight=false`
to omit those waits in a timing comparison. In the host case, Rango stays
at the module location between visits by its assigned Kiki and returns on
that pair's last trip for the module. Rango's recharge is a stationary service
at the host with its Kiki supplying power. Kiki itself can recharge at Atlante.
The replay exports the timeline, per-trip installation times, and a charge
ledger with the sender, receiver, and location. The cases are timing and
energy screens using saved trajectories and do not qualify new RPO epochs.
The `rangoWorkBetweenCharges_min` option defaults to 90 minutes per Rango.
Installations stay intact, so charging occurs before any installation that
would exceed this cumulative work limit. The EPS reserve can require earlier
charging, and a resident-host plan can still fail the standby reserve check.
Set `options.makePlots=false` for a trip-count sweep. The resulting
`study.deliveryCapacity` reports the number of trips completed before each
following Taxi handover and the remaining margin in minutes.
By default, `options.independentKikiPorts=true` allows both pairs to load,
recharge or fly concurrently. Host approaches, departure and installation
still reserve the shared assembly site. Set the option to false to compare
against the previous single-transit schedule. The timing replay does not
screen simultaneous flight paths for separation or replay Atlante's total
power demand for simultaneous charging.

Each completed assembly gets a detail timeline for Taxi and the assigned
Kiki/Rango pair. Dotted markers show arrival at added free-flight hold
points. All timelines share the mission time axis and are exported as PNGs.

## GMAT sunlight reference

Set `cfg.useGMAT = true` and `cfg.cluster = "ISS"` or `"SSO"` before
calling `run_faber_environment(budget,cfg)` with the corresponding Phase A
orbit. The study epoch and horizon must fit inside the 2038 SunOrbit report.
The planner interpolates the GMAT Sun vectors and computes illumination along
each simulated trajectory. The EclipseLocator intervals for the selected
reference orbit are available in `study.gmatEclipseIntervals` and its CSV
export. The trajectory-specific intervals remain in `study.eclipseIntervals`.

## Seasonal hold comparison

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
- Kiki is assumed mechanically attached after reaching Atlante or the
  assembly site. Attached waits can span eclipse. The schedule retains
  Kiki and Rango at the correct location while they wait for sunlit
  departure, deployment, handover or recovery. Intermediate free-flight
  gate holds still require drift and clearance checks.
- If the normal Kiki relocation cannot be scheduled, the scheduler can
  try up to two stationary V-bar holds in its gate-to-gate flight. It
  keeps the original departure and terminal paths and checks the new
  sublegs for Phase A speed and gate limits, obstacle clearance, and an
  impulsive burn-duration estimate. Each scheduled eclipse hold must
  also pass the drift and reacquisition checks. Set `assemblyExtraHolds`
  to 0 to disable the alternatives, and `assemblyHoldSublegMax_min` to
  bound a candidate subleg. `Environment_assemblyHoldPoints.csv` gives
  the Atlante-relative coordinates of the committed added gates.
- Transient free-flight waits, including waits at original gates, keep
  their full modeled uncertainty envelope at least 30 m from Atlante and
  the modeled attachment interfaces. The required distance is the larger
  of 30 m, `assemblyHoldKeepOut_m`, and `assemblyHoldKOZFactor` times the
  Phase A KOZ radius. Specify additional structure geometry through
  `assemblyHoldGeometry`. Attached waits at an interface use the attached
  state rather than the free-flight hold buffer.

`Sunlight_Eclipse_Comparison.csv` summarizes scheduled waits, delivery,
assembly and the strict sunlight audit. `Far_Field_Lighting.csv` lists
coast and burn-epoch lighting for each executed transfer. Each case
folder contains the original event, hold, lighting and EPS exports.
