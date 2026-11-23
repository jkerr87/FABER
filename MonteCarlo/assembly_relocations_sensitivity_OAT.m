
function S = assembly_relocations_sensitivity_OAT( ...
    mu,Re,h,Npts,timesMin,hostPoint,modulePoints,baseOpts)
% ASSEMBLY_RELOCATIONS_SENSITIVITY_OAT
%
% One-at-a-time sensitivity analysis for FABER assembly relocations.
%
% Evaluates the influence of:
%   1. Exit/entrance gate distance
%   2. Maximum controlled-exit speed
%   3. Maximum HCW transfer speed
%   4. Maximum controlled-arrival speed
%
% In each experiment, ONLY ONE parameter changes.
% All remaining parameters retain their baseline values.
%
% Uses existing assembly_relocations_v4.m.
%
% VEHICLE ALLOCATION:
%   Kiki  : modules 1,3,5,...
%   Rango : modules 2,4,6,...
%
% OUTPUTS:
%   S.baseline     Baseline assembly relocation result
%   S.table        Complete sensitivity comparison
%   S.phaseTable   DV breakdown by maneuver phase
%   S.runs         Individual relocation results
%   S.baseOpts     Baseline settings
%
% No graphs.
% No CSV exports.
% No Monte Carlo.
%
% All DV values in m/s.
% All TOF values in minutes.

%% ================================================================
% 1. INPUT VALIDATION
% ================================================================

if nargin < 8 || isempty(baseOpts)
    baseOpts = struct();
end

assert(exist('assembly_relocations_v4','file')==2, ...
    'assembly_relocations_v4.m not found.');

assert(size(modulePoints,2)==3, ...
    'modulePoints must be Nmodules x 3.');

assert(numel(hostPoint)==3, ...
    'hostPoint must contain three coordinates.');

assert(Npts>=3 && Npts==floor(Npts), ...
    'Npts must be an integer >= 3.');

assert(all(timesMin>0) && all(isfinite(timesMin)), ...
    'timesMin must contain positive finite durations.');

%% ================================================================
% 2. BASELINE SETTINGS
% ================================================================

% Defaults correspond to the current Phase A configuration.
% Any values provided through baseOpts override these defaults.

defaults = struct();

defaults.methods = {'GeneralHCW','Vbar','Rbar'};

defaults.gateDistance_m = 20;

defaults.exitTimesMin = 5:5:30;
defaults.finalTimesMin = 5:5:30;

defaults.maxExitSpeed_mps = 0.10;
defaults.maxTransferSpeed_mps = 0.30;
defaults.maxFinalSpeed_mps = 0.10;

defaults.maxExitAccel_mps2 = Inf;
defaults.maxFinalAccel_mps2 = Inf;

defaults.minTargetRange_m = 0;
defaults.printSummary = false;

opts = defaults;

fields = fieldnames(baseOpts);

for k = 1:numel(fields)
    opts.(fields{k}) = baseOpts.(fields{k});
end

% Prevent individual relocation printouts.
opts.printSummary = false;

%% ================================================================
% 3. DEFINE ONE-AT-A-TIME SWEEPS
% ================================================================

% Edit these arrays to change the sensitivity ranges.

gateDistances_m = [10 20 30 40];

exitSpeeds_mps = [0.05 0.10 0.15 Inf];

transferSpeeds_mps = [0.10 0.12 0.14 0.16 0.18 0.20 0.30];

arrivalSpeeds_mps = [0.05 0.10 0.15 Inf];

% Parameter descriptions and corresponding opts fields.

sweeps = struct( ...
    'name',{}, ...
    'field',{}, ...
    'unit',{}, ...
    'values',{});

sweeps(1).name = 'Gate distance';
sweeps(1).field = 'gateDistance_m';
sweeps(1).unit = 'm';
sweeps(1).values = gateDistances_m;

sweeps(2).name = 'Exit speed';
sweeps(2).field = 'maxExitSpeed_mps';
sweeps(2).unit = 'm/s';
sweeps(2).values = exitSpeeds_mps;

sweeps(3).name = 'HCW speed';
sweeps(3).field = 'maxTransferSpeed_mps';
sweeps(3).unit = 'm/s';
sweeps(3).values = transferSpeeds_mps;

sweeps(4).name = 'Arrival speed';
sweeps(4).field = 'maxFinalSpeed_mps';
sweeps(4).unit = 'm/s';
sweeps(4).values = arrivalSpeeds_mps;

%% ================================================================
% 4. RUN BASELINE
% ================================================================

fprintf('\n========================================\n');
fprintf('FABER ASSEMBLY SENSITIVITY ANALYSIS\n');
fprintf('========================================\n');

fprintf('Modules: %d\n',size(modulePoints,1));
fprintf('Altitude: %.2f km\n',h);

fprintf('\nBASELINE CONSTRAINTS:\n');
fprintf('Gate: %.2f m\n',opts.gateDistance_m);

fprintf('Exit speed: %.3f m/s\n', ...
    opts.maxExitSpeed_mps);

fprintf('HCW speed: %.3f m/s\n', ...
    opts.maxTransferSpeed_mps);

fprintf('Arrival speed: %.3f m/s\n', ...
    opts.maxFinalSpeed_mps);

B0 = assembly_relocations_v4( ...
    mu,Re,h,Npts,timesMin, ...
    hostPoint,modulePoints,opts);

Nmodules = size(modulePoints,1);

valid0 = isfinite(B0.DV_perModule) & ...
         isfinite(B0.TOF_perModule_min);

if ~all(valid0)

    error(['Baseline is infeasible: only %d/%d ', ...
        'modules have complete round trips.'], ...
        nnz(valid0),Nmodules);

end

baselineDV = sum(B0.DV_perModule);
baselineTOF = sum(B0.TOF_perModule_min);

fprintf('\nBASELINE RESULTS:\n');
fprintf('Total DV: %.6f m/s\n',baselineDV);
fprintf('Total TOF: %.2f min\n',baselineTOF);

%% ================================================================
% 5. PREALLOCATE SENSITIVITY RESULTS
% ================================================================

% Include one baseline row.
nRuns = 1;

for k = 1:numel(sweeps)

    field = sweeps(k).field;

    values = sweeps(k).values;

    baselineValue = opts.(field);

    % Exclude repeated baseline values.
    nRuns = nRuns + nnz(values~=baselineValue);

end

Scenario = strings(nRuns,1);
Parameter = strings(nRuns,1);
Value = NaN(nRuns,1);
Unit = strings(nRuns,1);

FeasibleModules = zeros(nRuns,1);
Feasible = false(nRuns,1);

KikiDV_mps = NaN(nRuns,1);
RangoDV_mps = NaN(nRuns,1);
TotalDV_mps = NaN(nRuns,1);

KikiTOF_min = NaN(nRuns,1);
RangoTOF_min = NaN(nRuns,1);
TotalTOF_min = NaN(nRuns,1);

DVchange_pct = NaN(nRuns,1);
TOFchange_pct = NaN(nRuns,1);

ExitDV_mps = NaN(nRuns,1);
HCW_DV_mps = NaN(nRuns,1);
ArrivalDV_mps = NaN(nRuns,1);

% Retain full results for every simulation.
runs = cell(nRuns,1);

%% ================================================================
% 6. EVALUATE BASELINE
% ================================================================

row = 1;

Scenario(row) = "Baseline";
Parameter(row) = "Baseline";
Value(row) = NaN;
Unit(row) = "-";

runs{row} = B0;

[ ...
    FeasibleModules(row), ...
    Feasible(row), ...
    KikiDV_mps(row), ...
    RangoDV_mps(row), ...
    TotalDV_mps(row), ...
    KikiTOF_min(row), ...
    RangoTOF_min(row), ...
    TotalTOF_min(row), ...
    ExitDV_mps(row), ...
    HCW_DV_mps(row), ...
    ArrivalDV_mps(row) ...
] = summarizeAssembly(B0);

DVchange_pct(row) = 0;
TOFchange_pct(row) = 0;

%% ================================================================
% 7. RUN ONE-AT-A-TIME SENSITIVITY
% ================================================================

for p = 1:numel(sweeps)

    name = sweeps(p).name;
    field = sweeps(p).field;
    unit = sweeps(p).unit;
    values = sweeps(p).values;

    fprintf('\n----------------------------------------\n');
    fprintf('SWEEP: %s\n',name);
    fprintf('----------------------------------------\n');

    for k = 1:numel(values)

        value = values(k);

        % Baseline already calculated.
        if value == opts.(field)
            continue;
        end

        row = row+1;

        % Start from original baseline every time.
        optsTest = opts;

        % Modify ONLY the parameter being examined.
        optsTest.(field) = value;

        Scenario(row) = name + " = " + num2str(value);
        Parameter(row) = name;
        Value(row) = value;
        Unit(row) = unit;

        fprintf('Testing %s = %.4g %s ...\n', ...
            name,value,unit);

        try

            B = assembly_relocations_v4( ...
                mu,Re,h,Npts,timesMin, ...
                hostPoint,modulePoints,optsTest);

            runs{row} = B;

            [ ...
                FeasibleModules(row), ...
                Feasible(row), ...
                KikiDV_mps(row), ...
                RangoDV_mps(row), ...
                TotalDV_mps(row), ...
                KikiTOF_min(row), ...
                RangoTOF_min(row), ...
                TotalTOF_min(row), ...
                ExitDV_mps(row), ...
                HCW_DV_mps(row), ...
                ArrivalDV_mps(row) ...
            ] = summarizeAssembly(B);

            if Feasible(row)

                DVchange_pct(row) = 100 * ...
                    (TotalDV_mps(row)-baselineDV) ...
                    /baselineDV;

                TOFchange_pct(row) = 100 * ...
                    (TotalTOF_min(row)-baselineTOF) ...
                    /baselineTOF;

                fprintf('  DV: %.6f m/s\n', ...
                    TotalDV_mps(row));

                fprintf('  TOF: %.2f min\n', ...
                    TotalTOF_min(row));

                fprintf('  DV change: %+.2f %%\n', ...
                    DVchange_pct(row));

                fprintf('  TOF change: %+.2f %%\n', ...
                    TOFchange_pct(row));

            else

                fprintf('  Infeasible: %d/%d modules\n', ...
                    FeasibleModules(row),Nmodules);

            end

        catch ME

            % Do not silently treat numerical errors
            % as physical infeasibility.

            fprintf('  ERROR: %s\n',ME.message);

            Scenario(row) = ...
                Scenario(row) + " [ERROR]";

        end

    end
end

%% ================================================================
% 8. MAIN SENSITIVITY TABLE
% ================================================================

SensitivityTable = table( ...
    Scenario, ...
    Parameter, ...
    Value, ...
    Unit, ...
    FeasibleModules, ...
    Feasible, ...
    KikiDV_mps, ...
    RangoDV_mps, ...
    TotalDV_mps, ...
    KikiTOF_min, ...
    RangoTOF_min, ...
    TotalTOF_min, ...
    DVchange_pct, ...
    TOFchange_pct);

fprintf('\n========================================\n');
fprintf('FINAL OAT SENSITIVITY RESULTS\n');
fprintf('========================================\n');

disp(SensitivityTable);

%% ================================================================
% 9. PHASE-BY-PHASE DELTA-V BREAKDOWN
% ================================================================

PhaseTable = table( ...
    Scenario, ...
    ExitDV_mps, ...
    HCW_DV_mps, ...
    ArrivalDV_mps, ...
    TotalDV_mps);

fprintf('\n========================================\n');
fprintf('ASSEMBLY DELTA-V BY PHASE\n');
fprintf('========================================\n');

disp(PhaseTable);

%% ================================================================
% 10. PRINT INDIVIDUAL PARAMETER COMPARISONS
% ================================================================

for p = 1:numel(sweeps)

    idx = Parameter==sweeps(p).name;

    fprintf('\n========================================\n');
    fprintf('%s SENSITIVITY\n',upper(sweeps(p).name));
    fprintf('========================================\n');

    disp(SensitivityTable(idx, ...
        {'Value','Feasible', ...
         'TotalDV_mps','TotalTOF_min', ...
         'DVchange_pct','TOFchange_pct'}));

end

%% ================================================================
% 11. OUTPUT STRUCTURE
% ================================================================

S = struct();

S.baseline = B0;

S.table = SensitivityTable;
S.phaseTable = PhaseTable;

S.runs = runs;

S.baseOpts = opts;
S.sweeps = sweeps;

S.baselineDV_mps = baselineDV;
S.baselineTOF_min = baselineTOF;

fprintf('\n========================================\n');
fprintf('SENSITIVITY ANALYSIS COMPLETE\n');
fprintf('========================================\n');

end

%% ================================================================
% LOCAL FUNCTION: SUMMARIZE ASSEMBLY RESULTS
% ================================================================

function [ ...
    nFeasible,allFeasible, ...
    kikiDV,rangoDV,totalDV, ...
    kikiTOF,rangoTOF,totalTOF, ...
    exitDV,hcwDV,arrivalDV ...
] = summarizeAssembly(B)

N = B.Nmodules;

kikiIdx = 1:2:N;
rangoIdx = 2:2:N;

valid = isfinite(B.DV_perModule) & ...
        isfinite(B.TOF_perModule_min);

nFeasible = nnz(valid);

allFeasible = all(valid);

% Never report partial-campaign sums as complete budgets.

kikiDV = NaN;
rangoDV = NaN;
totalDV = NaN;

kikiTOF = NaN;
rangoTOF = NaN;
totalTOF = NaN;

exitDV = NaN;
hcwDV = NaN;
arrivalDV = NaN;

if ~allFeasible
    return;
end

% Vehicle-level Delta V
kikiDV = sum(B.DV_perModule(kikiIdx));

rangoDV = sum(B.DV_perModule(rangoIdx));

totalDV = kikiDV+rangoDV;

% Vehicle-level time of flight
kikiTOF = sum(B.TOF_perModule_min(kikiIdx));

rangoTOF = sum(B.TOF_perModule_min(rangoIdx));

totalTOF = kikiTOF+rangoTOF;

% Phase Delta V, outbound and return combined
exitDV = sum(B.DV_exit_mps(:));

hcwDV = sum(B.DV_transfer_mps(:));

arrivalDV = sum(B.DV_final_mps(:));

% Internal consistency check
tol = 1e-6*max(1,totalDV);

assert(abs(totalDV-(exitDV+hcwDV+arrivalDV))<tol, ...
    'Phase DV components do not match total assembly DV.');

end
