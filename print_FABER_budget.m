
function B = print_FABER_budget(out)
% PRINT_FABER_BUDGET
%
% Prints and exports the FABER Phase A maneuver budget.
%
% Margin policy:
%   General: DV + max(5% DV, 10 m/s)
%   RPO:     100% margin (2 x raw DV)
%
% Includes:
%   - Initial orbit correction
%   - Inter-client orbit transfer
%   - Disposal
%   - External relative station keeping
%   - General SK + parking acquisition
%   - RPO to Client
%   - RPO to ATLANTE
%   - Assembly relocations (Kiki and Rango)
%   - Momentum unloading (all vehicles)
%
% Usage:
%   out = MA_vehicle_budget_phaseA();
%   B = print_FABER_budget(out);

%% SETTINGS

T = out.phaseTable;
SK = out.annualSK;

Nyears = out.cfg.missionYears;

% External relative SK allocation
relativeSK_mpsYr = 40;

% Choose one parking architecture
parkingArch = "Vbar";

% Number of complete assembly campaigns.
% Current assembly result already covers 50 relocations.
NassemblyCampaigns = 1;

% Excel output
filename = 'FABER_PhaseA_Maneuver_Budget.xlsx';

%% INITIALIZE TABLE

data = cell(0,9);

%% ================================================================
% INITIAL ORBIT CORRECTION
% ================================================================

for veh = ["ATLANTE","Taxi"]

    idx = string(T.Vehicle)==veh & ...
          string(T.Segment)=="Deployment";

    dv = sum(T.DV_perEvent_mps(idx),'omitnan');
    tof = max(T.TOF_perEvent_min(idx),[],'omitnan');

    addRow("Initial orbit correction",veh, ...
        1,dv,tof,"General", ...
        "Calculated from altitude, inclination and RAAN");

end

%% ================================================================
% CLIENT ORBIT TRANSFER
% ================================================================

for veh = ["ATLANTE","Taxi"]

    idx = string(T.Vehicle)==veh & ...
          string(T.Segment)=="Inter-client";

    dv = sum(T.DV_perEvent_mps(idx),'omitnan');
    tof = max(T.TOF_perEvent_min(idx),[],'omitnan');

    addRow("Client orbit transfer",veh, ...
        out.cfg.Nclients,dv,tof,"General", ...
        "Extreme altitude and plane transfer");

end


%% ================================================================
% FAR-FIELD RENDEZVOUS - OUTBOUND AND INBOUND
% ================================================================
%
% Combines:
%   1. Hohmann altitude correction
%   2. Inclination correction
%   3. RAAN correction
%   4. Phasing maneuver
%
% All components are already calculated in
% MA_vehicle_budget_phaseA.m.
%
% Margin: max(5%, 10 m/s), applied to the
% combined Delta-V per rendezvous event.

farSegments = ["Outbound far","Inbound far"];

farNames = [ ...
    "Far-field rendezvous to Client", ...
    "Far-field rendezvous to ATLANTE"];

for j = 1:2

    idx = string(T.Vehicle)=="Taxi" & ...
        string(T.Segment)==farSegments(j);

    assert(any(idx), ...
        'Missing far-field phase: %s',farSegments(j));

    % Verify all components have the same event count
    counts = T.Count(idx);

    assert(all(counts==counts(1)), ...
        'Inconsistent far-field event counts.');

    % Total Delta-V per rendezvous
    dvFar = sum(T.DV_perEvent_mps(idx),'omitnan');

    % Total TOF:
    % Sum Hohmann transfer and phasing durations.
    % Inclination and RAAN corrections have zero
    % additional duration in the current model.
    tofFar = sum(T.TOF_perEvent_min(idx),'omitnan');

    % Add one consolidated budget row
    addRow( ...
        farNames(j), ...
        "Taxi", ...
        counts(1), ...
        dvFar, ...
        tofFar, ...
        "General", ...
        "Hohmann + inclination + RAAN + phasing");

end


%% ================================================================
% DISPOSAL
% ================================================================

for veh = ["ATLANTE","Taxi"]

    idx = string(T.Vehicle)==veh & ...
          string(T.Segment)=="End of life" & ...
          string(T.Phase)=="Deorbit";

    if ~any(idx)
        continue;
    end

    dv = sum(T.DV_perEvent_mps(idx),'omitnan');
    tof = max(T.TOF_perEvent_min(idx),[],'omitnan');

    addRow("Disposal",veh, ...
        1,dv,tof,"General", ...
        "Nominal impulsive deorbit");

end

%% ================================================================
% RELATIVE STATION KEEPING - EXTERNAL ALLOCATION
% ================================================================

for veh = ["ATLANTE","Taxi"]

    addRow("Relative SK (external)",veh, ...
        Nyears,relativeSK_mpsYr,NaN,"General", ...
        "40 m/s/year; externally allocated");

end

%% ================================================================
% GENERAL SK + PARKING ACQUISITION
% ================================================================

for veh = ["Taxi","ATLANTE"]

    idx = string(SK.Vehicle)==veh & ...
          string(SK.Architecture)==parkingArch;

    assert(nnz(idx)==1, ...
        'Expected exactly one parking scenario for %s.',veh);

    % Absolute orbital drag compensation over operating life
    dvDrag = SK.DragLifetime_mps(idx);

    % Total acquisition DV already accounts for all entries
    dvAcquisition = SK.EntryLifetime_mps(idx);

    % Relative SK intentionally excluded
    dv = dvDrag + dvAcquisition;

    addRow("SK general + parking acquisition",veh, ...
        1,dv,NaN,"General", ...
        "Lifetime drag + parking entries; relative SK excluded");

end

%% ================================================================
% RPO - SUM ALL SEVEN PHASES INTO TWO ROWS
% ================================================================

rpoSegments = ["Outbound RPO","Inbound RPO"];

rpoNames = ["RPO to Client","RPO to ATLANTE"];

for j = 1:2

    idx = string(T.Vehicle)=="Taxi" & ...
          string(T.Segment)==rpoSegments(j);

    if ~any(idx)
        continue;
    end

    counts = T.Count(idx);

    assert(all(counts==counts(1)), ...
        'Inconsistent RPO phase counts.');

    % Sum phases of the selected P95-DV trajectory
    dv = sum(T.DV_perEvent_mps(idx),'omitnan');

    tof = sum(T.TOF_perEvent_min(idx),'omitnan');

    addRow(rpoNames(j),"Taxi", ...
        counts(1),dv,tof,"RPO", ...
        "100% margin; summed selected RPO phases");

end

%% ================================================================
% ASSEMBLY RELOCATIONS - PER VEICOLETTO
% ================================================================

for veh = ["Kiki","Rango"]

    idx = string(T.Vehicle)==veh & ...
          string(T.Segment)=="Assembly";

    if ~any(idx)
        continue;
    end

    % Existing source rows represent aggregate campaign costs.
    % Count is one in the original phase ledger.

    campaignDV = sum(T.DV_lifetime_mps(idx),'omitnan');
    campaignTOF = sum(T.ActivityTime_min(idx),'omitnan');

    addRow("Assembly relocations",veh, ...
        NassemblyCampaigns,campaignDV, ...
        campaignTOF,"General", ...
        "One campaign contains 25 outbound-return relocations");

end

%% ================================================================
% MOMENTUM UNLOADING - PER VEHICLE
% ================================================================

for veh = ["ATLANTE","Taxi","Kiki","Rango"]

    idx = string(T.Vehicle)==veh & ...
          string(T.Segment)=="ADCS" & ...
          string(T.Phase)=="Momentum unloading";

    if ~any(idx)
        continue;
    end

    dv = sum(T.DV_lifetime_mps(idx),'omitnan');

    addRow("Momentum unloading",veh, ...
        1,dv,NaN,"General", ...
        "Calculated for allocated ADCS operating duration");

end

%% ================================================================
% CREATE OUTPUT TABLE
% ================================================================

B = cell2table(data, ...
    'VariableNames',{ ...
        'Maneuver', ...
        'Vehicle', ...
        'Count', ...
        'RawPerEvent_mps', ...
        'MarginedPerEvent_mps', ...
        'MarginedLifetime_mps', ...
        'TOFperEvent_min', ...
        'MarginType', ...
        'Notes'});

%% ================================================================
% DISPLAY RESULTS
% ================================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf('FABER PHASE A - MANEUVER DELTA V BUDGET\n');
fprintf('====================================================\n');
fprintf('General margin: max(5%%, 10 m/s)\n');
fprintf('RPO margin: 100%%\n');
fprintf('Relative SK: externally allocated\n');
fprintf('Parking architecture: %s\n',parkingArch);
fprintf('====================================================\n\n');

disp(B(:,1:7));

%% ================================================================
% TOTALS BY VEHICLE
% ================================================================

vehicleNames = ["Taxi";"ATLANTE";"Kiki";"Rango"];

rawTotals = zeros(numel(vehicleNames),1);
marginedTotals = zeros(numel(vehicleNames),1);

for j = 1:numel(vehicleNames)

    idx = B.Vehicle==vehicleNames(j);

    rawTotals(j) = sum( ...
        B.Count(idx).*B.RawPerEvent_mps(idx), ...
        'omitnan');

    marginedTotals(j) = sum( ...
        B.MarginedLifetime_mps(idx),'omitnan');

end

VehicleTotals = table( ...
    vehicleNames,rawTotals,marginedTotals, ...
    'VariableNames',{ ...
        'Vehicle', ...
        'RawLifetime_mps', ...
        'MarginedLifetime_mps'});

fprintf('\n========================================\n');
fprintf('TOTAL DELTA V BY VEHICLE\n');
fprintf('========================================\n');

disp(VehicleTotals);

%% ================================================================
% EXCEL EXPORT
% ================================================================

if isfile(filename)
    delete(filename);
end

writetable(B,filename, ...
    'Sheet','Maneuver Budget');

writetable(VehicleTotals,filename, ...
    'Sheet','Vehicle Totals');

writetable(T,filename, ...
    'Sheet','Original Phase Ledger');

writetable(SK,filename, ...
    'Sheet','Station Keeping');

if isfield(out,'momentumUnloading') && ...
        ~isempty(out.momentumUnloading)

    writetable(out.momentumUnloading,filename, ...
        'Sheet','Momentum Unloading');

end

if isfield(out,'RPO_MC_PhaseTable')

    writetable(out.RPO_MC_PhaseTable,filename, ...
        'Sheet','RPO Details');

end

fprintf('\nExcel exported: %s\n',filename);

%% ================================================================
% NESTED FUNCTION - ADD A BUDGET ROW
% ================================================================

    function addRow(name,vehicle,count,dv,tof, ...
            marginType,note)

        % Apply margins per event
        if strcmpi(string(marginType),"RPO")

            % 100% margin
            dvMargined = 2*dv;

        else

            % General margin:
            % larger of 5% raw DV or 10 m/s
            dvMargined = dv + max(0.05*dv,10);

        end

        lifetimeDV = count*dvMargined;

        data(end+1,:) = { ...
            string(name), ...
            string(vehicle), ...
            count, ...
            dv, ...
            dvMargined, ...
            lifetimeDV, ...
            tof, ...
            string(marginType), ...
            string(note)};

    end

end
