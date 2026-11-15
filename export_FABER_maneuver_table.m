
function M = export_FABER_maneuver_table(out, filename, options)
% EXPORT_FABER_MANEUVER_TABLE
%
% Creates an Excel-ready FABER mission Delta-V budget.
%
% USAGE:
%   out = MA_vehicle_budget_phaseA();
%   M = export_FABER_maneuver_table(out);
%
% OPTIONAL:
%   options = struct();
%   options.architecture = 'Vbar';
%   options.externalRelativeSK_mpsYr = 40;
%   options.assemblyRelocationsPerVehicle = 25;
%
%   M = export_FABER_maneuver_table( ...
%       out,'FABER_Maneuver_Budget.xlsx',options);
%
% MAIN FEATURES:
% - Uses current mission analysis results
% - Consolidates injection corrections
% - Consolidates inter-client orbit transfers
% - Preserves all Taxi RPO phases
% - Separates Kiki and Rango relocations
% - Includes momentum unloading per vehicle
% - Includes parking drag and acquisition
% - Lists external relative SK separately
% - Avoids counting both parking architectures
% - Exports results directly to Excel

%% DEFAULT SETTINGS

if nargin < 2 || isempty(filename)
    filename = 'FABER_Maneuver_Budget.xlsx';
end

if nargin < 3 || isempty(options)
    options = struct();
end

if ~isfield(options,'architecture')
    options.architecture = 'Vbar';
end

if ~isfield(options,'externalRelativeSK_mpsYr')
    options.externalRelativeSK_mpsYr = 40;
end

if ~isfield(options,'assemblyRelocationsPerVehicle')
    options.assemblyRelocationsPerVehicle = 25;
end

%% LOAD SOURCE RESULTS

T = out.phaseTable;

assert(istable(T), ...
    'out.phaseTable must be a MATLAB table.');

requiredFields = { ...
    'Vehicle', ...
    'Segment', ...
    'Phase', ...
    'Architecture', ...
    'Count', ...
    'DV_perEvent_mps', ...
    'TOF_perEvent_min', ...
    'IncludeInBaseline'};

assert(all(ismember( ...
    requiredFields,T.Properties.VariableNames)), ...
    'Expected mission ledger fields missing.');

mfrac = out.cfg.marginFraction;

%% INITIALIZE OUTPUT ARRAYS

S = strings(0,1);
V = strings(0,1);
Seg = strings(0,1);

C = zeros(0,1);
D = zeros(0,1);
ToF = zeros(0,1);

Status = strings(0,1);
Basis = strings(0,1);

Included = false(0,1);
Notes = strings(0,1);

used = false(height(T),1);

%% ================================================================
% CONSOLIDATE INJECTION AND ORBIT TRANSFER CORRECTIONS
% ================================================================

vehicles = ["Taxi","ATLANTE"];

for v = vehicles

    % Initial orbit correction:
    % Altitude + inclination + RAAN

    mask = string(T.Vehicle)==v & ...
        strcmpi(string(T.Segment),'Deployment');

    if any(mask)

        appendMerged( ...
            "Initial orbit correction", ...
            v, ...
            "Deployment", ...
            mask, ...
            "Altitude + inclination + RAAN corrections");

    end

    % Inter-client orbit transfer:
    % Altitude + plane change

    mask = string(T.Vehicle)==v & ...
        strcmpi(string(T.Segment),'Inter-client');

    if any(mask)

        appendMerged( ...
            "Client orbit transfer", ...
            v, ...
            "Inter-client", ...
            mask, ...
            "Altitude + plane-change corrections");

    end

end

%% ================================================================
% IMPORT ALL OTHER MISSION MANEUVERS
% ================================================================

for i = 1:height(T)

    if used(i)
        continue;
    end

    v = string(T.Vehicle(i));
    seg = string(T.Segment(i));
    p = string(T.Phase(i));
    arch = string(T.Architecture(i));

    %% PARKING SCENARIO SELECTION

    if strcmpi(seg,'Parking')

        % Relative station keeping is externally budgeted.
        % Do not include it again here.

        if contains(lower(p),'relative maintenance')
            continue;
        end

        % Keep only one parking architecture.

        desired = v + " " + ...
            string(options.architecture);

        if arch ~= desired
            continue;
        end

    end

    %% EXTRACT DATA

    count = T.Count(i);
    dv = T.DV_perEvent_mps(i);
    tof = T.TOF_perEvent_min(i);

    label = p;
    note = "";

    %% ASSEMBLY RELOCATIONS

    if strcmpi(seg,'Assembly') && ...
            (v=="Kiki" || v=="Rango") && ...
            (contains(lower(p),'relocation') || ...
             contains(lower(p),'return'))

        % Original mission ledger contains total
        % campaign Delta-V with Count = 1.
        %
        % Convert to average Delta-V per relocation.
        % This preserves the total campaign budget.

        trips = options.assemblyRelocationsPerVehicle;

        assert( ...
            isscalar(trips) && ...
            trips>0 && ...
            trips==floor(trips), ...
            'Assembly relocation count must be positive integer.');

        if count == 1

            count = trips;
            dv = dv/trips;
            tof = tof/trips;

            note = ...
                "Average per relocation; lifetime DV unchanged";

        else

            warning( ...
                'Assembly row %d already has Count ~= 1.',i);

        end

    end

    %% PARKING LABELS

    if strcmpi(seg,'Parking')

        label = p + " (" + ...
            string(options.architecture) + ")";

        note = ...
            "Parking scenario; Vbar and orbit are alternatives";

    end

    %% ADD ROW

    append( ...
        label, ...
        v, ...
        seg, ...
        count, ...
        dv, ...
        tof, ...
        string(T.Status(i)), ...
        string(T.Basis(i)), ...
        logical(T.IncludeInBaseline(i)), ...
        note);

end

%% ================================================================
% EXTERNAL RELATIVE STATION KEEPING
% ================================================================

% Relative SK is already accounted for externally.
% It is shown for traceability but excluded
% from the internal vehicle budget.

if options.externalRelativeSK_mpsYr > 0

    for v = vehicles

        append( ...
            "Relative station keeping (external)", ...
            v, ...
            "External SK", ...
            out.cfg.missionYears, ...
            options.externalRelativeSK_mpsYr, ...
            NaN, ...
            "Externally budgeted", ...
            "External literature assumption", ...
            false, ...
            "Not included in MATLAB baseline totals");

    end

end

%% ================================================================
% CALCULATE MARGINS AND LIFETIME TOTALS
% ================================================================

MarginRate = repmat( ...
    mfrac,numel(C),1);

MarginedPerEvent = D.*(1+MarginRate);

RawLifetime = C.*D;

MarginedLifetime = C.*MarginedPerEvent;

LifetimeTOF = C.*ToF;

%% ================================================================
% CREATE MAIN EXPORT TABLE
% ================================================================

M = table( ...
    S, ...
    V, ...
    Seg, ...
    C, ...
    D, ...
    MarginRate, ...
    MarginedPerEvent, ...
    RawLifetime, ...
    MarginedLifetime, ...
    ToF, ...
    LifetimeTOF, ...
    Included, ...
    Status, ...
    Basis, ...
    Notes, ...
    'VariableNames',{ ...
        'Maneuver', ...
        'Vehicle', ...
        'Segment', ...
        'Count', ...
        'RawPerEvent_mps', ...
        'MarginRate', ...
        'MarginedPerEvent_mps', ...
        'RawLifetime_mps', ...
        'MarginedLifetime_mps', ...
        'TOF_perEvent_min', ...
        'LifetimeTOF_min', ...
        'IncludeInBaseline', ...
        'Status', ...
        'Basis', ...
        'Notes'});

%% CONSISTENCY CHECK

valid = isfinite(M.RawLifetime_mps);

assert(all( ...
    abs( ...
        M.RawLifetime_mps(valid) - ...
        M.Count(valid).* ...
        M.RawPerEvent_mps(valid)) < 1e-7), ...
    'Lifetime Delta-V consistency check failed.');

%% ================================================================
% EXCEL EXPORT
% ================================================================

if isfile(filename)
    delete(filename);
end

% Main maneuver budget
writetable( ...
    M,filename, ...
    'Sheet','Maneuver Budget');

% Original detailed mission ledger
writetable( ...
    T,filename, ...
    'Sheet','Source Phase Ledger');

% Parking results
if isfield(out,'annualSK') && ...
        istable(out.annualSK)

    writetable( ...
        out.annualSK,filename, ...
        'Sheet','Parking Scenarios');

end

% Momentum unloading
if isfield(out,'momentumUnloading') && ...
        istable(out.momentumUnloading)

    if ~isempty(out.momentumUnloading)

        writetable( ...
            out.momentumUnloading,filename, ...
            'Sheet','Momentum Unloading');

    end

end

% RPO phase results
if isfield(out,'RPO_MC_PhaseTable') && ...
        istable(out.RPO_MC_PhaseTable)

    writetable( ...
        out.RPO_MC_PhaseTable,filename, ...
        'Sheet','RPO Phases');

end

% Original vehicle totals
if isfield(out,'vehicleTotals') && ...
        istable(out.vehicleTotals)

    writetable( ...
        out.vehicleTotals,filename, ...
        'Sheet','Original Vehicle Totals');

end

%% ================================================================
% BASELINE TOTALS BY VEHICLE
% ================================================================

% Only include rows explicitly marked as baseline.
% External relative SK and alternative parking
% scenarios are excluded.

select = M.IncludeInBaseline & ...
    isfinite(M.RawLifetime_mps);

summary = table( ...
    ["Taxi";"ATLANTE";"Kiki";"Rango"], ...
    zeros(4,1), ...
    zeros(4,1), ...
    'VariableNames',{ ...
        'Vehicle', ...
        'BaselineRawDV_mps', ...
        'BaselineMarginedDV_mps'});

for j = 1:height(summary)

    keep = select & ...
        M.Vehicle==summary.Vehicle(j);

    summary.BaselineRawDV_mps(j) = ...
        sum(M.RawLifetime_mps(keep));

    summary.BaselineMarginedDV_mps(j) = ...
        sum(M.MarginedLifetime_mps(keep));

end

writetable( ...
    summary,filename, ...
    'Sheet','Baseline By Vehicle');

%% PRINT RESULTS

fprintf('\n========================================\n');
fprintf('FABER MANEUVER BUDGET EXPORTED\n');
fprintf('========================================\n');

fprintf('Excel file: %s\n',filename);
fprintf('Total maneuver rows: %d\n',height(M));

fprintf('Parking architecture: %s\n', ...
    string(options.architecture));

fprintf('External relative SK: %.2f m/s/year\n', ...
    options.externalRelativeSK_mpsYr);

fprintf('Relative SK excluded from baseline totals.\n');

disp(summary);

%% ================================================================
% LOCAL FUNCTION - MERGE MANEUVERS
% ================================================================

    function appendMerged(label,v,segment,mask,notesText)

        idx = find(mask);

        counts = T.Count(idx);

        if any(abs(counts-counts(1))>1e-10)
            error( ...
                'Cannot merge %s: inconsistent event counts.', ...
                label);
        end

        addDV = sum(T.DV_perEvent_mps(idx));

        % Components of a single orbit correction
        % do not require separate transfer times.

        ts = T.TOF_perEvent_min(idx);

        if all(isnan(ts))
            mergedTOF = NaN;
        else
            mergedTOF = max(ts,[],'omitnan');
        end

        append( ...
            label, ...
            v, ...
            segment, ...
            counts(1), ...
            addDV, ...
            mergedTOF, ...
            "Estimated", ...
            "Summed from mission component ledger", ...
            all(logical(T.IncludeInBaseline(idx))), ...
            notesText);

        used(idx) = true;

    end

%% ================================================================
% LOCAL FUNCTION - ADD ONE MANEUVER
% ================================================================

    function append(label,v,segment,count,dv,tof, ...
            status,basis,included,note)

        S(end+1,1) = string(label);
        V(end+1,1) = string(v);
        Seg(end+1,1) = string(segment);

        C(end+1,1) = count;
        D(end+1,1) = dv;
        ToF(end+1,1) = tof;

        Status(end+1,1) = string(status);
        Basis(end+1,1) = string(basis);

        Included(end+1,1) = logical(included);
        Notes(end+1,1) = string(note);

    end

end
