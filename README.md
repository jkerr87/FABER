# FABER

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
   git clone https://github.com/jkerr87/ASMAD.git
   cd ASMAD
   ```

2. Open MATLAB and set **Current Folder** to the cloned `ASMAD` root.
3. Place the required external helper folders where the script expects them, or adjust its `addpath` statements. To make the checked-in folders available recursively, you can run:

   ```matlab
   addpath(genpath(pwd));
   % Add external dependency directories as needed:
   % addpath(genpath('C:/path/to/ORBMECH'));
   % addpath(genpath('C:/path/to/PM_V2'));
   ```

4. Edit **Section 1 — Constants and Mission Inputs** of `MA_deltav_v2.m` and the RPO settings further down the script. In particular, select the intended `structureType` and make sure its configuration script exists.
5. Run in the MATLAB Command Window:

   ```matlab
   MA_deltav_v2
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
| Minimum major-event margin | `margin_minimum` | 10 | m/s |

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
- The default delta-V margin is 5%. For **injection, client orbit transfer and disposal**, the added allowance is the larger of 5% or 10 m/s. Abort and momentum-unloading line items are **doubled**. This policy is preliminary, not derived from Monte Carlo confidence intervals.
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

Source: https://github.com/jkerr87/ASMAD

For traceability, record the selected structure, orbit settings, MATLAB version and Git commit hash with every published result.
