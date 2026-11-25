# FABER

## Integrated Taxi RPO and mission budget

### Delta-V margin policy

Each deterministic powered event receives an added allowance of the larger
of 5% of its raw delta-V or 10 m/s. The floor is applied once to a complete
event. For injection, inter-client transfer, and each far-field rendezvous,
the detailed component rows receive proportional shares of that event's
allowance. The assembly allowance is applied per powered leg per trip,
before summing the campaign. Events with zero modeled delta-V remain zero.
The event allowance is applied before multiplying by the event count.

Stochastic RPO and contingency lines receive a 100% allowance. The integrated
Monte Carlo budget doubles each accepted RPO cycle's cost and combines it
with the already-margined non-RPO ledger. This produces raw and margined
conditional P95 results without applying a second allowance to the total.
Unresolved costs remain open. The rule applies to delta-V and leaves time
and trial feasibility unchanged.

### Per-vehicle phase tables

```matlab
cfg = struct('vehicleBudgetFolder', 'FABER_Per_Vehicle_Tables');
out = MA_vehicle_budget_phaseA(cfg);
```

The function returns and optionally exports one table each for Taxi, Atlante,
Kiki and Rango. Maneuver costs are quoted once per event, including each
far-field leg, RPO approach, inter-client transfer, and parking entry or
departure. Vehicle subtotals group like reporting units and sum one of
each listed event, regardless of campaign size. The default
`cfg.interClientRAANReserve_mps = 100` is added once to the margined
inter-client transfer cost per vehicle. Module retrieval uses the
modeled plane difference and the usual maneuver margin. Taxi close-field rows
come from the paired trial nearest the accepted-cycle delta-V P95. Kiki
relocation costs are the mean per sampled leg, with momentum unloading
quoted over one relocation round trip. Rango is carried by Kiki
during orbital transfers, so the Rango table flags unresolved surface
operations and the Kiki table has an open host-access transport line.

Phase A evaluates end-of-life disposal using `PassiveDeorbit/PassiveDisposal.m`
at the configured circular orbit altitude. The default budgets the low solar
activity profile. The model assigns zero delta-V if natural decay reaches
the cutoff within five years, or an insertion burn followed by drag decay
otherwise. Set the disposal masses and effective projected areas in
`FABER_phaseA_defaults` or override them for a specific release configuration.
Atlante's reference mass includes the two attached Kiki/Rango pairs. Both
density profiles and their elapsed decay times are reported in `out.disposal`,
`Passive_Disposal.csv` and the workbook's **Passive Disposal** tab.

The baseline budgets Taxi's initial 60-minute continuous-thrust RCS
release from 5 m to 350 m behind Atlante. Taxi's parking departure
alignment remains in the parking scenarios and phase ledger. The
release budget row is a standalone estimate. The outbound Lambert
still starts at the configured parking handoff state.

The Phase A export includes nominal RPO checks, first-rejection counts,
per-phase ranges and speeds, burn use, and a finite-thrust replay of the
budget representative and the trajectories at each nominal leg P95. Use
`rpoFiniteValidation="all"` to replay every
nominally passing pair or `"none"` to skip the replay. The replay tests
powered phases from their nominal starts, using the coupled arrival masses
and local chief altitudes. It reports finite errors separately from the
nominal delta-V allocation. Set mission-specific geometry, terminal corridor,
or cycle and burn limits through the `rpo...` inputs in `FABER_phaseA_defaults`.
The checks report absent geometry and optional limits as unconfigured.
The Excel workbook lists the nominal Client and Atlante RPO leg P95s in
the maneuver budget and the leg P95s by architecture and acceptance population in
**RPO P95**. Detailed cases, phase rows, and checks remain in the MATLAB
result and the optional CSV export. The workbook's **Relative SK** tab shows
the raw annual rate and the raw and margined totals over the modeled exposure.
**RPO Propulsion** reports the powered phases of each leg's P95 trajectory.
It records the selected trial, thruster group, configured active and installed
counts, rated thrust, Isp, and the modeled firing seconds from finite replay.
The full replay adds the P95 trajectories among the finite-thrust passes.

To also refresh the `MA_DV_Vehicle.xlsx` report layout, pass the original
template path and a separate destination to Phase A. The run refreshes
the ISS-like tab for a prograde orbit or the SSO-like tab for a retrograde
orbit. Run the second orbit with the same destination to fill its tab.
The heading, row layout and formatting come from the template. Managed
cells contain modeled values, including separate nominal RPO leg P95s.
Unresolved costs and times remain blank. In the ToF column, modeled
maneuver times are in minutes and disposal decay is explicitly in years.
Each orbit also gets its own RPO propulsion tab.
The **Transfer Firing** tabs cover injection and client changes, parking
maneuvers, the initial Taxi release, both burns on each Lambert leg,
disposal insertion and the mean loaded/unloaded Kiki relocation. The Taxi
release uses its simulated firing time. Powered impulsive Taxi transfers
use equivalent firing time at the configured thrust and Isp, excluding
coast time and delta-V margin. The Lambert rows carry the mass between
burns and use the captured-stack mass at return departure. Atlante uses
LEROS 4-ET main propulsion for injection, client changes and disposal,
and HBT-1 RCS for parking. Kiki uses DST-12 thrusters for translation.
The defaults assume one active Atlante main engine and four active thrusters
in each RCS group, with installed hardware reported separately. Override
`atlanteTransferPropulsion`, `atlanteParkingPropulsion`, or
`kikiTransferPropulsion` with the confirmed firing allocation. Supply
`type`, `nActiveThrusters`, `thrust_N` per thruster, `Isp_s`, and
`thrustEfficiency`, with optional `nInstalled`.

```matlab
cfg = struct('vehicleReportTemplate','C:/path/to/MA_DV_Vehicle.xlsx', ...
    'vehicleReportWorkbook',fullfile(pwd,'MA_DV_Vehicle_Simulation.xlsx'));
iss = MA_vehicle_budget_phaseA(cfg);
cfg.h_target_km = 680.12;
cfg.i_target_deg = 98.12;
cfg.hClientSpread_km = 72.84;
cfg.iClientSpread_deg = 0.31;
sso = MA_vehicle_budget_phaseA(cfg);
```

The original Phase A workbook still exports with `print_FABER_budget(iss)`.
The CSV option `vehicleBudgetFolder` can also be used in the same run.

Phase A uses `ORBMECH/lambertMR.m` to connect Taxi's parking departure
handoff to the module-relative RPO gate `cfg.x0_client`. The module starts
`phaseAngle_deg` ahead in argument of latitude. Set the module orbit
offsets relative to Atlante at round trip departure. The transfer uses a
two-body circular target ephemeris and budgets both the departure burn
and the gate velocity match. The propagated post-burn state seeds each
outbound RPO trial. Those trials still sample approach design variables.

See `out.outboundLambert` for the transfer and its inertial/LVLH handoff.
The search defaults to `phasingRevs` revolutions and compares candidate
durations and branches above `outboundMinimumPerigee_km`. Override the
search with `outboundLambertTOF_s` and `outboundLambertRevolutions`.
`outboundDepartureEpoch_s` sets the scenario clock at the pre-burn handoff.
The initial orbit orientations are set at that epoch.

The coupled return uses a five-minute capture dwell by default. Override
`captureWait_min` for a different scenario. Each successful nominal
client approach supplies its terminal LVLH state and elapsed time. Phase A
models instantaneous capture by combining Taxi and module position and
momentum, then propagates the resulting stack and Atlante through the
specified dwell. A second Lambert transfer targets `cfg.x0_atlante` at
the achieved arrival epoch. The return RPO then uses the propagated
post-burn state and the mass after the return burns. Both far-field legs
estimate fuel using `farFieldIsp_s`.

The selected pair's capture, departure and arrival are in
`out.returnLambert`. Use `returnLambertTOF_s` and
`returnLambertRevolutions` to override the default search. With the
default unknown dwell, the return rows stay open and Phase A can still
report a nominal outbound-only RPO budget. The model assumes ideal
capture at the RPO endpoint and an impulsive departure from the stack
state. Local separation and departure clearance, perturbations, and
finite-burn validation remain separate from this calculation.

### Coupled Taxi round trip safety and sensitivity

Run the study independently of the maneuver budget. The example below saves
CSV files and plots in the requested folder.

```matlab
cfg = struct('N',30,'outputFolder','FABER_Round_Trip_Study');
cfg.phaseA = struct('captureWait_min',45);
study = main_FABER_round_trip_safety(cfg);
```

Each case runs `N` trials. The default includes a nominal control,
dispersions, missed burns and capture failure. Directed faults default
to zero dispersions so the study can reach and isolate each injection
point. Set `directedUsesDispersions = true` to combine them with the same
draws used in the dispersed case. Set `cfg.cases` to a string array to
choose a subset or include `"Outbound departure partial"` and
`"Return departure partial"`. Matched trials also share an RPO seed. Each accepted client
approach supplies its achieved endpoint and time to capture and dwell.
The return is retargeted from that trial's stack and Atlante state, and
its executed post-burn state and mass seed the return RPO. Faults and
rejected gates leave downstream states open instead of substituting a
nominal arrival.

Input dispersions are illustrative 1-sigma values. `sigmaInjectionLVLH`
and `sigmaModuleLVLH` are six-vectors in m and m/s, applied at the round trip
departure to Taxi's residual injection state and the module state.
`sigmaDepartureTime_s` slips departure, `sigmaBurnScale` is a fractional
magnitude error, and `sigmaBurnVector_mps` adds independent inertial vector
components at each far-field burn. A missed departure cancels the
scheduled arrival burn. `truthJ2` can be set to `1.08262668e-3` to stress
the two-body plan with J2 in inertial execution. Close-field propagation
uses the existing nominal HCW model. Use `cfg.rpo` for existing
feasibility inputs such as `clientGeometry`, `hostGeometry`,
`terminalCorridorHalfWidth_m`, and speed or propulsion limits. Set
`gatePositionTolerance_m`, `gateVelocityTolerance_mps`,
`farProtectedRadius_m`, `contactRadius_m`, and `passiveHorizon_s` to
mission-specific screening values.

`study.trials` contains stage, reason, epochs, gate residuals, closest
approach, and executed delta-V. P95 cost and time use complete round trips.
`study.sensitivity` reports descriptive rank correlations and sample
counts for dispersed trials. Later-stage metrics are conditional on
reaching that stage. The plots show LVLH ensembles, terminal paths and
configured boundaries, Lambert burn delta-V versus time, gate error
scatter, range and closing speed, passive drift, and distributions.

`cfg.rpoFaultTrials` controls an additional conditional fault screen on
recorded burn and coast events. The screen uses dispersed trajectories
when available and falls back to the nominal control. The default is
30 samples per approach when events are available. CAM searches require
verified net firing-group vectors in
`cfg.rpoSafety.camForceVectors_LVLH_N`. Without those vectors, the screen
reports passive drift and leaves CAM unassessed. CAM and optional recovery
are separate branches from the round trip simulation. Directed case counts
and conditional fault outcomes do not represent operational fault
probabilities.

The model uses nominal RPO retargeting, ideal attachment, and illustrative
clearance and uncertainty inputs. It supports preliminary trajectory
screening. A navigation estimator, closed-loop controller, capture
mechanics, and vehicle-specific qualification are separate work.

### Environment, sunlight scheduling and EPS duty cycle

From the ASMAD root, pass a Phase A result with its selected retrieval and
assembly trajectories into the environment study.

```matlab
addpath('Environment');

% Illustrative handover, loading and dwell durations in minutes.
budget = MA_vehicle_budget_phaseA(struct('N_RPO_MC',30, ...
    'captureWait_min',45,'handoverWait_min',10,'veicLoading_min',5));

cfg = struct('epochUTC',datetime(2038,1,1,'TimeZone','UTC'), ...
    'assemblyModuleIndices',1,'kikiHandover_min',5, ...
    'outputFolder','FABER_Environment_Study');
study = run_faber_environment(budget,cfg);
```

The epoch is the origin of the Phase A retrieval timeline. The schedule
maintains one explicit cargo inventory and operates one Kiki/Rango pair
at a time. Select real assembly interface indices, use `preloadedModules`
for independent assembly work, and use `modulesPerRetrieval` to specify
the cargo delivered by the default single retrieval. Unknown handover and
loading durations stop the dependent schedule and appear in
`study.openItems`.

For a complete two-delivery example with one active assembly location,
run from the ASMAD root:

```matlab
addpath('Environment');
study = run_faber_environment_campaign('FABER_Environment_Campaign');
% Use 3 for the 20-, 15-, and 10-degree case.
three = run_faber_environment_campaign('FABER_Environment_Three',3);
```

For five modules, run `run_faber_environment_campaign('FABER_Environment_Five',5)`.
The two Kiki/Rango pairs alternate one assembly site at a time. Kiki
flies each loaded and empty relocation with Rango retained. Rango
installs for 32.17 minutes from the PM3/PM4 duty cycle while Kiki waits
at the site. The schedule uses the PM4 5-minute cargo load, handover and
Atlante electrical berthing, plus the 16.67-minute deployment and
recovery proxies. Relocation times come from the Phase A trajectories.

The default charging planner checks modeled energy and reserve through
each sortie and adds charging at Atlante only when required. The
optional `chargeAsNeeded=false` setting restores the fixed charge
windows. Docked standby power is accounted for separately.
See `study.assemblySorties`, `study.chargeSummary`, and
`study.epsSummary` for the schedule and power checks.

The two-delivery example uses 20-degree and 10-degree initial phases.
Both examples use a June 2038 epoch with Atlante RAAN at 180 degrees and
assume 5 minutes of active capture followed by 40 minutes attached. The mission duty cycle and its
detail tabs all share the same hours-since-epoch time axis. Reevaluate
sunlight availability when changing the epoch or orbit geometry.

To schedule more than one capture, supply the retrievals in flight order.

```matlab
cfg.retrievals = struct('moduleId',{'Module A','Module B'}, ...
    'phaseAngle_deg',{budget.cfg.phaseAngle_deg,10},'cargoCount',{1,1});
cfg.assemblyModuleIndices = [1 2];
cfg.captureActive_min = 5;
study = run_faber_environment(budget,cfg);
study.retrievalSummary
```

Each angle is the module's argument of latitude offset from Atlante at
`cfg.epochUTC`. All modules share the client orbit geometry and Phase A
payload mass and dwell assumptions. Set the angles to the actual module
states for the chosen epoch. After each delivery, Taxi aligns from its
achieved Atlante berth to the parking departure gate in sunlight. The
alignment uses an impulsive HCW model and sampled clearance checks.
`taxiDeparture_min` sets its modeled duration. This preliminary departure
model needs finite-burn and docking-interface validation.

The next outbound Lambert solve uses the alignment end state and epoch.
Taxi carries its remaining wet mass into each retrieval and sheds the
Phase A payload mass at handover. No intermediate propellant top-up is
modeled. Supply `taxiDryMass_kg` to apply a dry-mass limit. Assembly uses
the cargo available at each handover, and one EPS replay spans all
retrievals and assembly events. `study.retrievalSummary` reports delivery
status, times, delta-v, and mass by retrieval. A blocked or pending trip
has no reported final mass. `study.retrievalArcs` stores outbound and
return transfers in columns one and two, and `study.returnLambert`
exposes the most recent attempted return solve when available. The full
duty-cycle timeline marks each completed delivery and opens a detail
figure for each started retrieval.

The reference replay uses the selected Phase A RPO durations without lighting delays. The gated
replay requires full sunlight on the spacecraft and reference target
during active RPO, Kiki relocation, and handovers. It preserves
fly-through gates and reports a blocked operation when no suitable
continuous window fits or its preceding hold cannot be used safely.
It retargets the return Lambert arc after outbound lighting delays.

If the original Taxi path is blocked, the environment study tries
stationary V-bar entry holds and shorter nine-phase candidates under the
mission's existing feasibility limits. It schedules the terminal closing
and fly-through as one continuous sunlit chain on both approaches.
Set `replanSunlightRPO=false` to replay only the original Phase A design.
The fallback plan is selected for the environment study and does not
replace the saved Phase A budget. See `study.rpoPlanning` for the outcome.

When Taxi is blocked by lighting, the study also tests delays to its
first Lambert departure after trying the normal and shorter RPO paths.
`taxiDepartureSearch_min` sets the candidate delays, with defaults of
15, 30, 45 and 60 minutes. Use `[]` to disable the search or
`taxiInitialWait_s` to specify an initial wait. The parking drift and
clearance must pass before Taxi retargets from the propagated parking
state. See `study.rpoPlanning` for the delay selected. The closing chain
still reserves one sunlight window through the following sunlit capture
or handover.

Before Kiki departs the loading berth or installation site, the
scheduler checks the complete three-leg relocation and each intervening
hold, including any required reacquisition. If a preview fails, it
checks later departures from the attached location. Set
`assemblyDepartureSearchStep_s` to refine the search and
`maxWait_s` to set the departure wait limit. Only a passing preview
enters the event ledger and EPS simulation. If none passes, Kiki stays
attached and `study.openItems` reports the blocking leg.

By default the entire Phase A capture duration needs sunlight. If the
capture duration includes a confirmed attached wait, set
`captureActive_min` to the part that requires sunlight. For example, a
five-minute active capture and forty-minute attached wait use
`captureWait_min=45` and `captureActive_min=5`. The example is an input
assumption to validate against robotics operations.

A relative hold may span eclipse if sampled HCW drift, including the
specified initial-state error box, stays inside the allowed displacement
and outside the keep-outs. The sampled screen checks approach-zone
membership for RPO holds. Where a delay changes the nominal departure
state, a sunlit two-impulse reacquisition must pass the clearance, speed
and delta-V limits before resuming. These checks assume the supplied
geometry and sample resolution. Use `clientHoldGeometry`,
`atlanteHoldGeometry` and `assemblyHoldGeometry` to add actual
spheres or capsule geometry.

The study reports sun aspect relative to `aspectAxisLVLH` (default
along-track +V), orbit beta, sunlight, penumbra and umbra. One standalone
window contains all plots as tabs. The trajectory tabs compare the
reference and sunlight-gated timelines. The Kiki and Rango assembly
timeline shows the scheduled assembly work and reports any blocker, such
as a missing Taxi handover or a failed eclipse-hold screen. If available,
the Aerospace Toolbox eclipse audit uses the sampled inertial
ephemerides and writes its intervals separately. A Toolbox or ephemeris
failure is reported, and the geometric schedule uses the approximate
Sun ephemeris when the planetary data are unavailable.

The EPS replay uses the supplied PM4 hardware and power-mode tables,
constrained bus dispatch, battery efficiencies and rainflow cycling.
It uses measured lighting from the analysis, EOL reserve limits, and
Atlante-to-Kiki and Kiki-to-Rango charging. By default it assumes the
proposed retained-pair architecture and available docked standby power.
The Rango deploy/recover entries are attached-interface time proxies,
not independent flight paths. Set `retainRangoForRelocations=false` for
a prepositioned Rango, or `dockedPairSupport=false` to assess
battery-only standby. Enable `includeRefueling` for the optional PM4
fluid-service windows and top-ups. `chargeSummary` records port energy,
battery energy gain, and any residual time to the requested SoC.

`study.ready` requires both a complete lighting schedule and a passing
EPS energy/reserve replay. Check `study.turnaroundReady` separately for
whether all spacecraft reach the configured end-of-cycle SoC target.
Charge windows are configurable and the replay flags shortfalls so
they can be revised. Battery lifetime outputs
give observed depth of discharge and equivalent cycles. Set
`journeyFraction` only when the fraction of mission exposure represented
by the cycle is known. An exposure estimate is not a calibrated battery
life prediction. The LVLH aspect axis is a reference direction. PV
defaults to PM4 pointing factors unless `pvUsesAspect` is selected.

### Parking departure and relative stationkeeping

Phase A estimates a two-impulse HCW transfer from each parking option to
a pre-burn LVLH handoff. The default handoff is 100 m behind the reference
at zero relative velocity. Set `cfg.taxiParkingDepartureState` or
`cfg.atlanteParkingDepartureState` to supply a different six-component
state in m and m/s. `cfg.parkingDepartureTime_s` defaults to 3600 s.
The vehicle-specific `taxiParkingWaitBeforeDeparture_s` and
`atlanteParkingWaitBeforeDeparture_s` advance the ellipse phase from
acquisition to departure. Each scenario reports an acquisition and a
departure per event in `out.annualSK`.
The far-field or inter-client budget carries the subsequent orbit-transfer
burn. The HCW alignment is a budget estimate pending trajectory clearance
and thrust verification.

`MA_deltav.m` runs the same differential-drag function for Taxi relative
to Atlante and Atlante relative to the host, with separate lifetime budget
rows. Both use the full mission duration by default. Set
`missionParkingConfig.TaxiYears` and `missionParkingConfig.AtlanteYears`
for scenario-specific occupancy. The script uses its original idle Taxi
properties and phase-A assumptions for the Atlante stack and host.

Phase A runs the differential-drag model from
`NominalOps/run_FABER_Taxi_idle_SK.m` for Taxi relative to Atlante and for
Atlante relative to the host using the configured masses, drag areas, and
drag coefficients. Taxi uses the independently supplied
`cfg.taxiParkingYears`, which defaults to continuous parking over
`cfg.missionYears`. Atlante uses the full mission operating period. Each vehicle has one relative stationkeeping row
in the baseline. Both parking alternatives share that allowance because
the constant-drag HCW model produces the same corrections for either
geometry. The minimum sampled separation and approach-zone violations are
reported in `out.annualSK` and `Parking_Stationkeeping.csv`.

Set `cfg.taxiParkingEllipse` and `cfg.atlanteParkingEllipse` independently.
Each accepts `ellipseCenter_m` (positive distance behind the chief),
`ellipseRadialAmplitude_m`, and `ellipsePhase_rad`. For example:

```matlab
cfg.taxiParkingEllipse = struct('ellipseCenter_m',450, ...
    'ellipseRadialAmplitude_m',25,'ellipsePhase_rad',0);
cfg.atlanteParkingEllipse = struct('ellipseCenter_m',250, ...
    'ellipseRadialAmplitude_m',20);
out = MA_vehicle_budget_phaseA(cfg);
```

The standalone function takes two vehicle structs and an options struct:
`SK = run_FABER_Taxi_idle_SK(follower,reference,opts)`. Both vehicles
require `mass_kg`, `area_m2` and `Cd`, and can carry separate relative
ellipse geometry in a common circular-chief frame. The follower pays
for the corrections. Use `opts.altitude_m` and `opts.rho_kgm3` to specify
the reference orbit and constant density. The model uses a 2:1 in-plane
HCW ellipse and periodic two-impulse corrections, with a six-hour interval
by default. The Atlante–host approach-zone radius is unresolved by default.

The absolute drag allocation is still listed separately. Summing the two
allocations is conservative if the same burns can serve both purposes.
The model has no jointly targeted burn schedule or eccentric-chief
dynamics. The legacy aggregate `MA_deltav` budget still represents its own
idle configuration and retains its original stationkeeping model.

### Integrated mission totals

Run `MA_deltav` from the ASMAD root. The script writes the deterministic
budget and two additional CSVs: `FABER_MC_RPO_Trials.csv` (one row per
paired close-field cycle, including rejection reason) and
`FABER_MC_Mission_Trials.csv` (one row per resampled delivery campaign).
It uses the nine-phase safety package explicitly to avoid MATLAB path
ambiguity. The budget removes both deterministic close-field approaches
before adding the sampled paired cycle once per delivery. Remaining line
items, including legacy abort and CAM allowances, are deterministic.
Reported P95 delta-V and known time are conditional on successful cycles.
Mission acceptance assumes independent, identically distributed deliveries;
uncertain costs after a rejected attempt, shared biases, retries, and
statistical confidence bounds require separate models.
If any assembly round trip is unresolved, the affected assembly rows and
full mission delta-V remain NaN. The output separately reports a P95 subtotal
for known rows, excluding the unresolved assembly campaign.

The default Taxi architecture follows the locally available FABER PM3
report (ASV-3.pdf, Figure 1 and Table 5). The default 2,000 kg payload
and 10% impulse-time screen follow the PM3 sizing case. The output includes
the radial/along-track projection with AZ/KOZ and HP1–HP4, three-dimensional
trajectories, range and radial closing speed, acceptance and delta-V versus
TOF. The representative cycle is one actual accepted trial nearest the
conditional cycle delta-V P95. PM3's 36-architecture study, 1,000 scenarios,
and wet masses and orbits need to be reproduced separately before comparing
published percentiles.

Optional MATLAB workspace inputs:

- `runMissionMC = false` runs only the deterministic budget.
- `missionMCConfig` overrides fields from the nine-phase `ASV_RPO_MC_config`.
  Use `missionArchitecture` to select a column of `MC_RPO.methods` when
  supplying multiple architectures. This column is fixed across deliveries.
- `missionMCConfig.arrivalStatesClient` and `.arrivalStatesHost` accept
  paired 6-by-N target-centered LVLH states in metres and m/s. Preserve
  trial indexing and generate each host state from the matching client
  realization where possible. `FABER_relative_state` converts synchronous
  inertial and rotating RTN/LVLH states in both directions. Current
  analytical far-field budget legs still provide costs and durations only.
- `missionMCConfig.clientGeometry`, `.hostGeometry`, and
  `assemblyOptions.pathGeometry` take the same sphere/capsule fields as
  `FABER_clearance_m`, with optional `extraRadius_m` for footprint and
  margin. Only include obstacles whose clearance is mandatory on the whole
  checked path. Authorized surfaces need a separate approach/contact model.
- `missionMCConfig.terminalCorridorHalfWidth_m` selects a negative V-bar
  lateral clearance check inside the KOZ. Supply a mission-derived width.
- `assemblyOptions.maxTransferDV_mps` caps the sum of the gate-to-gate
  impulses. Existing speed and controlled-acceleration caps remain
  configurable. The budget now uses v4 and accounts for departure and
  arrival times as well as the coast. Inspect
  `assemblyResult.rejectionByDirection` for the first failed segment by
  interface, method, and direction (outbound, return).

Model limitations for the PM3 output set:

1. To close the handoff, carry a common epoch, target and chaser inertial
   states, vehicle mass, and maneuver execution errors through both
   legs. Perturb inertial injection/target orbit and burn errors, propagate
   both spacecraft, and transform the resulting correlated states at
   the entry event. Propagate the capture/release and departure states
   before solving the next far-field leg. Handoffs built from independent
   offsets miss these correlations and cannot quantify interface error.
2. The nominal RPO model assumes GO decisions and ideal holds. Zones and
   centerline clearance alone cannot assess timed gates, hold uncertainty,
   navigation outages, loss of control, or retreat. Keep PM3 Table 6's
   conditional fault and CAM analysis separate from the reserve until fault
   likelihoods, recovery paths, and effective firing groups are known.
3. Assembly feasibility still uses sampled translation and ideal impulses.
   Add as-built geometry per installation step, loaded/unloaded swept
   volumes, docking/berthing contact rules, attitude-dependent force and
   torque, and missed-burn drift. Refine near-boundary trajectories with
   event detection. The current default geometry and acceleration limits
   require mission inputs before making clearance or propulsion claims.
4. PM3 Table 15 and Appendix G require relative navigation and pointing.
   Add attitude and target geometry for line of sight, field of view,
   occlusion and lighting. An estimator/control simulation should log
   truth versus estimate, covariance consistency, tracking error and
   actuator saturation before assessing capture tolerances.
5. The lifetime delta-V sum combines Taxi and assembly line items for
   bookkeeping. Use vehicle-specific mass histories and Isp for propellant
   sizing, and a resource-constrained schedule for mission duration.
   Neither the raw P95 nor summed known TOF is a vehicle propellant or
   calendar-duration estimate.

The rest of this README describes the earlier deterministic study and
its historical Monte Carlo code.

# FABER — Facility for Autonomous Building and Expansion Robotics 

Preliminary MATLAB mission analysis for an Assembly and Servicing Vehicle (ASV), including a mission-level **delta-V budget**, transfer times, rendezvous/proximity operations (RPO), assembly relocations, and a Monte Carlo study of RPO strategies.

**Main entry point:** [`MA_deltav.m`](MA_deltav.m). The file is a MATLAB **script**, not a function; edit its input section and run it from the repository root.

> **Status:** preliminary engineering model. Results depend on scenario assumptions, chosen structure configuration, and helper code. They are not a flight-validated propellant budget or collision-safety demonstration.

## 1. Requirements and dependencies

- MATLAB (release/toolboxes not yet verified). The script uses standard MATLAB tables, plots, `writetable`, and local functions.
- The `astroConstants` function (used for Earth's gravitational parameter and radius), normally provided through the orbital-mechanics course/toolbox files.
- Helper functions and configuration files called by `MA_deltav_v2.m`. The script attempts to add `ORBMECH`, `PM_V2`, and `PM_V2/configs` to the MATLAB path.
- The repository also contains `HCW_solver.m`, `MonteCarlo/`, and `ProximityOps/`. Ensure all required helper functions are on the path. These directories **do not necessarily replace** the `ORBMECH` and `PM_V2` dependencies referenced in the script.

Examples of invoked helpers include `phasing_manoeuvre`, `proximity_operations_v2`, `rpo_holdingpt`, `assembly_relocations_v2`, `station_keeping`, `momentum_unloading`, `ASV_RPO_MC_config`, `ASV_RPO_MonteCarlo`, `FABER_RPO_validate_baseline`, and `FABER_RPO_analyze_results`. `hohmannTransfer` is defined locally at the end of the main script.

**Before running**, verify dependencies in MATLAB:

```matlab
which astroConstants
which phasing_manoeuvre
which proximity_operations_v2
which ASV_RPO_MC_config
which ASV_RPO_MonteCarlo
which solar_array
```

Each command should return a valid `.m` file path. If it reports a function is not found, add its directory or obtain the missing dependency. The GitHub repository's publicly visible top-level listing does not include `ORBMECH` or `PM_V2`; therefore a fresh clone may not run unmodified.

## 2. How to run

1. Clone/download the repository:

   ```bash
   git clone https://github.com/jkerr87/FABER.git
   cd FABER
   ```

2. Open MATLAB and set **Current Folder** to the cloned `FABER` root.
3. Place the required external helper folders where the script expects them, or adjust its `addpath` statements. To make the checked-in folders available recursively, you can run:

   ```matlab
   addpath(genpath(pwd));
   % Add external dependency directories as needed:
   % addpath(genpath('C:/path/to/ORBMECH'));
   % addpath(genpath('C:/path/to/PM_V2'));
   ```

4. Edit **Section 1 — Constants and Mission Inputs** of `MA_deltav.m` and the RPO settings further down the script. In particular, select the intended `structureType` and make sure its configuration script exists.
5. Run in the MATLAB Command Window:

   ```matlab
   MA_deltav
   ```

   Or open the script in the MATLAB Editor and click **Run**. `clear; clc; close all;` at the start resets the workspace, command window, and open figures.

**To run just the deterministic mission budget:** the final block beginning `%% RPO MONTE CARLO - OPTIMIZATION STUDY` can be commented out temporarily if Monte Carlo helpers are unavailable. Earlier parts of the script still call validation and other helper functions and therefore require those dependencies.

## 3. Scenario inputs (current values in the main script)

The following are **defaults from the checked-in script**, not universal requirements:

| Parameter | MATLAB variable | Current value | Units / meaning |
|---|---|---:|---|
| Target altitude | `h_target` | 457.25 | km |
| Target inclination | `i_target` | 51.02 | deg |
| Target RAAN | `RAAN_target` | 0 | deg |
| Disposal perigee altitude | `h_disposal` | 250 | km |
| Selected structure | `structureType` | `solar_array` | config script name |
| Client altitude spread | `h_SD` | 67.47 | km, ± from target |
| Client inclination spread | `i_SD` | 2.81 | deg, ± from target |
| Client RAAN spread | `RAAN_SD` | 0 | deg |
| Mission duration | `mission_years` | 7 | years |
| Planning delivery rate | `modules_per_year` | 5 | modules/year |
| Number of clients | `N_clients` | 3 | clients |
| Launcher altitude error | `alt_error_launch` | 15 | km |
| Launcher inclination error | `inc_error_launch` | 0.1 | deg |
| Launcher RAAN error | `RAAN_error_launch` | 0.1 | deg |
| Abort fraction | `abort_fraction` | 0.5 | fraction of deliveries |
| Avoidance events | `CAM_per_year` | 4 | per year |
| Abort delta-V | `DV_abort` | 1 | m/s per event |
| Avoidance delta-V | `DV_CAM` | 1 | m/s per event |
| General delta-V margin | `margin_fraction` | 0.05 | 5% |
| Minimum deterministic powered-event margin | `margin_minimum` | 10 | m/s |

**Important:** `N_modules_planned = mission_years * modules_per_year` is a *planning* value. The actual number of deliveries used in the budget is `N_deliveries = N_containers`, obtained from the selected structure configuration; `N_modules` is set to the same count. Edit the configuration script if you want to change the assembly architecture, container capacity, or installation positions.

The code contains commented alternative orbital values. Select a coherent scenario rather than mixing defaults from different mission cases.

### Proximity operations and assembly

The current client-side RPO profile is:

| Phase | Target or assumption | Duration |
|---|---|---:|
| Initial relative state | `[50; -2000; 20; 0.01; 0.04; -0.01]` (m, m/s) | — |
| Homing | stationary hold at LVLH `y = -1000 m` | 3600 s |
| Holding point 1 | `y = -1000 m` | 30 min |
| Closing | `Vbar` to `y = -100 m` | 1800 s |
| Holding point 2 | `y = -100 m` | 10 min |
| Final approach | `Vbar` to `y = -5 m` | 600 s |

These are **specified profile values**, not optimized maneuver times. The script also models return-to-host RPO using separately defined targets, and a seven-day illustrative Taxi assembly hold. Assembly transfers are computed from `modulePoints` in the selected configuration, using `assembly_relocations_v2`. The configured Veicoletto host point is `[0; -50; 0]` m and `timesMin = [5 10 15 20]` min.

## 4. Model methodology and assumptions

### Orbital transfers

- Earth-centered motion; the transfer visualization and Hohmann estimates use **ideal two-body impulsive mechanics**.
- Initial and target orbits for transfer calculations are treated as **circular**.
- Launcher injection correction is estimated from a target-altitude-plus-error injection orbit. Altitude, inclination, and RAAN corrections are computed **separately and summed**, rather than solving a combined minimum-delta-V plane-change problem.
- Client-orbit transfer covers extreme altitude/inclination/RAAN cases derived from the specified spreads. The full inter-plane angle accounts for inclination and RAAN together; it is not added again to the separate inclination-only estimate.
- Far-field phasing is evaluated using `phasing_manoeuvre` with `phase_deg = 180` and `N_rev = 100` in the checked-in scenario.
- A deorbit maneuver is **one retrograde burn** from the target circular orbit to an ellipse with perigee at `h_disposal`. **No circularization burn** at disposal perigee is included, and atmospheric re-entry is not propagated.

### Rendezvous, proximity operations, and return

- Close-field phases are represented through the project's proximity-operations routines, including homing, holding, closing and final approach. The main script specifies V-bar methods for closing/final approach; feasible trajectories and collision constraints must be assessed separately.
- For the far-field return leg, **J2 secular RAAN precession** is estimated for the host and captured module, assuming circular orbits, fixed other elements, and an identical initial epoch. J2 drift is applied to estimate the return-plane mismatch; this is not a full numerical J2 trajectory propagation.
- Captured-module inclination and RAAN are currently example offsets derived from launcher errors. Capture/handling time is set to **zero** (`T_capture = 0`).
- Some burns are impulsive and assigned zero duration; these are **not** estimates of real finite-thrust burn times.

### Stationkeeping, momentum unloading, and assembly

- Stationkeeping is sized at `SK_altitude = h_target - (alt_error_launch + h_SD)` as a conservative lower-altitude case.
- Taxi properties are estimates: `mass_Taxi_dry = 8256.93 kg`, `mass_Veicoletto = 875.17 kg`, `Cd_Taxi = 2.2`, and equivalent density `100 kg/m^3`. A cube-equivalent size determines reference frontal area and inertia.
- Atmospheric density is approximated using an **exponential model**, with reference density `9.158e-12 kg/m^3` at `350 km` and scale height `53.298 km` in the current script. This is an approximate, scenario-dependent model; it is not a solar-cycle-dependent atmosphere simulation.
- Momentum unloading uses approximate geometry/disturbance assumptions, including 10% inertia asymmetry, center-of-pressure offset of 10% of the equivalent body size, and `Cr = 1.5`.
- Assembly relocations are aggregated across all `modulePoints` from the chosen structure configuration. The assembly out/return totals already represent the configured assembly campaign; they are counted **once** in the lifetime budget, not once per delivery.

### Event counts and budget margins

- A single injection correction and a single disposal maneuver are budgeted.
- Client orbit transfers are multiplied by `N_clients`; outbound and inbound Taxi RPO/far-field phases by `N_deliveries`.
- Expected abort count is `abort_fraction * N_deliveries`; collision-avoidance count is `CAM_per_year * mission_years`. An expected count can be fractional—it is a budgeting assumption, not a literal event schedule.
- The current delta-V margin policy is specified above. Momentum unloading uses the deterministic rule, and aborts and collision avoidance use the stochastic rule.
- Some TOFs (e.g., aborts, collision avoidance, stationkeeping and momentum unloading) are `NaN`/unknown and **excluded** from the reported total modeled time. Therefore the displayed mission TOF is **not** a full seven-year mission timeline.
- The script does **not** convert the final delta-V budget into propellant mass or account for changing vehicle mass using a full rocket-equation mass budget.

## 5. Outputs

Running the script prints phase-specific delta-V and time-of-flight results to the MATLAB Command Window and builds `BudgetTable`, including:

- Maneuver name and event count
- Raw and margined delta-V **per event** [m/s]
- Raw and margined **lifetime** delta-V [m/s]
- TOF per event and aggregate known TOF [min]

It also calculates mission-level `DV_total_raw` and `DV_total_margin`, and produces horizontal bar charts for lifetime delta-V allocation and raw-versus-margined totals.

The main CSV is written to the **current MATLAB working directory**:

```text
ASV_SevenYear_DeltaV_Budget.csv
```

Despite the fixed filename, the contents reflect the currently configured mission duration and event counts. Change the filename in `writetable` if comparing multiple cases to avoid overwriting prior results.

The final section runs an RPO Monte Carlo study via `ASV_RPO_MonteCarlo`, analyzes results with `FABER_RPO_analyze_results`, and extracts 95th-percentile cycle delta-V/time (`DV_P95`, `TOF_P95`). It also plots transfer-count relationships. **These statistics are separate from the main deterministic `BudgetTable`**; the script does not automatically replace the baseline budget with its Monte Carlo P95 result.

## 6. Known limitations / future work

- Verify and package the external helper folders and document tested MATLAB version/toolbox requirements so a fresh clone is reproducible.
- Validate selected transfer geometries, line-of-sight, keep-out zones, approach-speed bounds and abort trajectories against mission safety requirements.
- Replace illustrative RPO durations and capture time with mission-specific operations scheduling.
- Evaluate optimized combined plane changes and more complete perturbation/propagation models as needed.
- Refine atmospheric density, geometry, mass evolution, stationkeeping and momentum-unloading allocations.
- Revisit contingency frequencies and margins with evidence or uncertainty analysis rather than fixed assumptions.
- Separate reusable mission settings from the main script into a single configuration file, and save each run under a scenario-specific filename.

## 7. Repository

Source: https://github.com/jkerr87/FABER

For traceability, record the selected structure, orbit settings, MATLAB version and Git commit hash with every published result.
