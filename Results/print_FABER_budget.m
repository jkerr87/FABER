function B = print_FABER_budget(out,filename)
% PRINT_FABER_BUDGET
%
% Prints and exports the FABER Phase A maneuver budget.
%
% Margin policy:
%   General: DV + max(5% D5V, 10 m/s)
%   RPO:     100% margin (2 x raw DV)
%
% Includes:
%   - Initial orbit correction
%   - Inter-client orbit transfer
%   - Disposal
%   - External relative station keeping
%   - Operating-period drag compensation and per-event parking
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
marginFloor = 10;
if isfield(out.cfg,'marginMinimum_mps')
    marginFloor = out.cfg.marginMinimum_mps;
end

% External relative SK allocation
relativeSK_mpsYr = 40;

% Choose one parking architecture
parkingArch = "Vbar";

% Number of complete assembly campaigns.
% Current assembly result already covers 50 relocations.

% Excel output
if nargin<2, filename = 'FABER_PhaseA_Maneuver_Budget.xlsx'; end

%% INITIALIZE TABLE

data = cell(0,8);

%% ================================================================
% INITIAL ORBIT CORRECTION
% ================================================================

for veh = ["ATLANTE","Taxi"]

    idx = string(T.Vehicle)==veh & ...
          string(T.Segment)=="Deployment";

    dv = sum(T.DV_perEvent_mps(idx));
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

    dv = sum(T.DV_perEvent_mps(idx));
    tof = max(T.TOF_perEvent_min(idx),[],'omitnan');

    if ismember('MarginedDV_perEvent_mps',T.Properties.VariableNames)
        margined = sum(T.MarginedDV_perEvent_mps(idx));
    else
        margined = FABER_apply_DV_margin(dv,false,out.cfg.marginFraction,marginFloor);
    end

    note = "Extreme altitude and plane transfer";
    reserve = idx & string(T.Phase)=="RAAN reserve";
    if any(reserve)
        note = note + sprintf(". Includes %.2f m/s RAAN reserve.", ...
            sum(T.MarginedDV_perEvent_mps(reserve)));
    end

    addRow("Client orbit transfer",veh,1,dv,tof,"General",note,margined);

end

%% INITIAL TAXI RELEASE FROM ATLANTE

initialDeparture = T.Vehicle=="Taxi" & string(T.Segment)=="Initial departure" & ...
    string(T.Phase)=="Atlante Departure";
if any(initialDeparture)
    addRow("Atlante Departure","Taxi",1,T.DV_perEvent_mps(initialDeparture), ...
        T.TOF_perEvent_min(initialDeparture),"General",string(T.Basis(initialDeparture)), ...
        T.MarginedDV_perEvent_mps(initialDeparture));
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
    dvFar = sum(T.DV_perEvent_mps(idx));

    basis = join(unique(string(T.Basis(idx)),'stable')," + ");
    tofFar = sum(T.TOF_perEvent_min(idx));

    % Add one consolidated budget row
    addRow( ...
        farNames(j), ...
        "Taxi", ...
        counts(1), ...
        dvFar, ...
        tofFar, ...
        "General", ...
        basis);

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

    dv = sum(T.DV_perEvent_mps(idx));
    tof = max(T.TOF_perEvent_min(idx),[],'omitnan');

    addRow("Disposal",veh, ...
        1,dv,tof,"General",join(string(T.Basis(idx)), " + "), ...
        sum(T.MarginedDV_perEvent_mps(idx)));

end

%% RELATIVE STATION KEEPING

relative = table();
for veh = ["ATLANTE","Taxi"]
    idx = string(T.Vehicle)==veh & string(T.Segment)=="Relative stationkeeping" ...
        & T.IncludeInBaseline;
    if any(idx)
        sk = SK(string(SK.Vehicle)==veh & string(SK.Architecture)==parkingArch,:);
        note = string(T.Basis(idx)) + sprintf(". %.6g m/s/year raw over %.3g years.", ...
            sk.RelativeAnnual_mps,sk.ParkingYears);
        addRow("Relative stationkeeping",veh,T.Count(idx), ...
            T.DV_perEvent_mps(idx),T.TOF_perEvent_min(idx),"General", ...
            note,T.MarginedDV_perEvent_mps(idx));
        annualEquivalent = NaN;
        if sk.ParkingYears>0, annualEquivalent = T.MarginedDV_perEvent_mps(idx)/sk.ParkingYears; end
        relative = [relative; table(veh,sk.ParkingYears,sk.RelativeAnnual_mps, ...
            annualEquivalent,T.DV_perEvent_mps(idx),T.MarginedDV_perEvent_mps(idx), ...
            'VariableNames',{'Vehicle','Exposure_years','RawDV_mps_per_year', ...
            'EquivalentMarginedDV_mps_per_year','ExposureRawDV_mps','ExposureMarginedDV_mps'})]; %#ok<AGROW>
        continue;
    end

    addRow("Relative SK (external)",veh, ...
        1,Nyears*relativeSK_mpsYr,NaN,"General", ...
        "40 m/s/year over the operating period; externally allocated");
end

%% ================================================================
% GENERAL SK + PARKING ACQUISITION
% ================================================================

for veh = ["Taxi","ATLANTE"]

    idx = string(SK.Vehicle)==veh & ...
          string(SK.Architecture)==parkingArch;

    assert(nnz(idx)==1, ...
        'Expected exactly one parking scenario for %s.',veh);

    addRow("Absolute stationkeeping",veh,1,SK.DragLifetime_mps(idx), ...
        NaN,"General","Drag compensation over the parking exposure");
    addRow("Parking acquisition",veh,1,SK.EntryDV_mps(idx), ...
        SK.EntryTOF_min(idx),"General","One parking acquisition");

    departure = T.Vehicle==veh & T.Architecture==veh+" "+parkingArch & ...
        string(T.Phase)=="Departure state alignment correction";
    % Taxi's parking alignment remains in the parking scenarios and ledger.
    if any(departure) && (veh~="Taxi" || ~any(initialDeparture))
        addRow("Parking departure alignment",veh,T.Count(departure), ...
            T.DV_perEvent_mps(departure),T.TOF_perEvent_min(departure), ...
            "General","HCW transfer to pre-burn handoff",T.MarginedDV_perEvent_mps(departure));
    end

end

%% RPO LEG P95

if isfield(out,'RPO_MC')
    mc = out.RPO_MC;
    paired = mc.valid(:,1);
    names = ["RPO to Client","RPO to ATLANTE"];
    note = "Leg P95 among nominal paired passes; independent time P95; 100% delta-V margin";
    for leg = 1:2
        dv = legP95(mc.DV_mps(paired,1,leg));
        tof = legP95(mc.TOF_min(paired,1,leg));
        addRow(names(leg),"Taxi",1,dv,tof,"RPO",note);
    end
end

%% ================================================================
% ASSEMBLY RELOCATIONS - PER VEICOLETTO
% ================================================================

idx = find(ismember(string(T.Vehicle),["Kiki","Rango"]) & ...
    string(T.Segment)=="Assembly" & T.IncludeInBaseline);
for k = idx'
    margined = NaN;
    if ismember('MarginedDV_perEvent_mps',T.Properties.VariableNames)
        margined = T.MarginedDV_perEvent_mps(k);
    else
        margined = FABER_apply_DV_margin(T.DV_perEvent_mps(k),false, ...
            out.cfg.marginFraction,marginFloor);
    end

    addRow(string(T.Phase(k)),T.Vehicle(k),T.Count(k),T.DV_perEvent_mps(k), ...
        T.TOF_perEvent_min(k),"General",string(T.Basis(k)),margined);
end

idx = find(string(T.Segment)=="Contingency" & T.IncludeInBaseline);
for k = idx'
    addRow(string(T.Phase(k)),T.Vehicle(k),T.Count(k),T.DV_perEvent_mps(k), ...
        T.TOF_perEvent_min(k),"Stochastic",string(T.Basis(k)));
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

    dv = sum(T.DV_lifetime_mps(idx));

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
        'BudgetUnit', ...
        'RawDV_mps', ...
        'MarginedDV_mps', ...
        'ModeledTime_min', ...
        'MarginType', ...
        'Notes'});

% The report-template exporter reuses the calculated rows without writing
% or printing the standard workbook.
if strlength(string(filename))==0, return; end

%% ================================================================
% DISPLAY RESULTS
% ================================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf('FABER PHASE A - MANEUVER DELTA V BUDGET\n');
fprintf('====================================================\n');
fprintf('General margin: max(5%%, 10 m/s)\n');
fprintf('RPO margin: 100%%\n');
fprintf('Parking architecture: %s\n',parkingArch);
fprintf('====================================================\n\n');

disp(B);

%% ================================================================
% TOTALS BY VEHICLE
% ================================================================

VehicleTotals = table();
for vehicle = unique(B.Vehicle,'stable')'
    for unit = unique(B.BudgetUnit(B.Vehicle==vehicle),'stable')'
        idx = B.Vehicle==vehicle & B.BudgetUnit==unit;
        VehicleTotals = [VehicleTotals; table(vehicle,unit, ...
            sum(B.RawDV_mps(idx),'omitnan'),sum(B.MarginedDV_mps(idx),'omitnan'), ...
            nnz(isnan(B.RawDV_mps(idx)) | isnan(B.MarginedDV_mps(idx))), ...
            'VariableNames',{'Vehicle','BudgetUnit','KnownRawDV_mps', ...
            'KnownMarginedDV_mps','UnresolvedDV'})];
    end
end

fprintf('\nBUDGET SUBTOTALS BY UNIT (one of each listed event)\n');
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
    'Sheet','Budget by Unit');

ledger = removevars(T,{'Count','DV_lifetime_mps','MarginedDV_mps','ActivityTime_min'});
ledger = renamevars(ledger,{'DV_perEvent_mps','MarginedDV_perEvent_mps','TOF_perEvent_min'}, ...
    {'RawDV_mps','MarginedDV_mps','ModeledTime_min'});
% The detailed trial and phase data remain available through the CSV export.
isRPO = ismember(string(ledger.Segment),["Outbound RPO","Inbound RPO"]);
writetable(ledger(~isRPO,:),filename,'Sheet','Phase Ledger');

writetable(SK,filename, ...
    'Sheet','Station Keeping');

if isfield(out,'disposal')
    writetable(out.disposal,filename,'Sheet','Passive Disposal');
end

if isfield(out,'momentumUnloading') && ...
        ~isempty(out.momentumUnloading)

    writetable(out.momentumUnloading,filename, ...
        'Sheet','Momentum Unloading');

end

if ~isempty(relative)
    writetable(relative,filename,'Sheet','Relative SK');
end

if isfield(out,'transferFiring')
    writetable(out.transferFiring,filename,'Sheet','Transfer Firing');
end

if isfield(out,'rpoChecks')
    writetable(rpoP95Table(out.rpoChecks),filename,'Sheet','RPO P95');
    if isfield(out.rpoChecks,'propulsion')
        writetable(out.rpoChecks.propulsion,filename,'Sheet','RPO Propulsion');
    end
end

fprintf('\nExcel exported: %s\n',filename);

%% ================================================================
% NESTED FUNCTION - ADD A BUDGET ROW
% ================================================================

    function addRow(name,vehicle,~,dv,tof, ...
            marginType,note,dvMargined)
        if nargin < 8
            stochastic = any(string(marginType)==["RPO","Stochastic"]);
            dvMargined = FABER_apply_DV_margin(dv,stochastic, ...
                out.cfg.marginFraction,marginFloor);
        end

        unit = "Per event";
        if contains(string(name),"stationkeeping","IgnoreCase",true) || ...
                string(name)=="Relative SK (external)" || ...
                (string(name)=="Momentum unloading" && ismember(string(vehicle),["Taxi","ATLANTE"]))
            unit = "Operating period";
        elseif string(name)=="Momentum unloading" && string(vehicle)=="Kiki"
            unit = "Per relocation round trip";
        end

        data(end+1,:) = { ...
            string(name), ...
            string(vehicle), ...
            unit, ...
            dv, ...
            dvMargined, ...
            tof, ...
            string(marginType), ...
            string(note)};

    end

end

function summary = rpoP95Table(checks)
cases = checks.cases;
architecture = checks.architectures;
rows = cell(0,9);
for j = 1:height(architecture)
    subset = cases(cases.Architecture==j,:);
    nominal = subset.NominalStatus=="Paired pass";
    passed = subset.FiniteStatus=="Pass";
    attempted = nnz(passed | subset.FiniteStatus=="Fail");
    label = join([architecture.methodH2(j),architecture.methodC1(j), ...
        architecture.methodC2(j),architecture.gateMode(j)]," / ");
    population = ["Nominal paired pass","Finite-thrust pass"];
    selection = {nominal,passed};
    coverage = [height(subset),attempted];
    for k = 1:2
        if k==2 && attempted==0, continue; end
        sample = subset(selection{k},:);
        client = legP95(sample.ClientRawDV_mps);
        atlante = legP95(sample.AtlanteRawDV_mps);
        rows(end+1,:) = {j,label,population(k),coverage(k),height(sample), ...
            client,2*client,atlante,2*atlante}; %#ok<AGROW>
    end
end
summary = cell2table(rows,'VariableNames',{'Architecture','Methods','Population', ...
    'TrialCount','CasesInP95','ClientRawP95_mps','ClientMarginedP95_mps', ...
    'AtlanteRawP95_mps','AtlanteMarginedP95_mps'});
end

function value = legP95(values)
values = sort(values(isfinite(values)));
value = NaN;
if ~isempty(values)
    value = values(ceil(0.95*numel(values)));
end
end
