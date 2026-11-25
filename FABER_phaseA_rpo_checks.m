function checks = FABER_phaseA_rpo_checks(R,selected,mode)
% Audit nominal feasibility and replay finite burns on budget candidates.
mode = string(validatestring(mode,{'representative','all','none'}));
cfg = R.cfg;
paired = R.valid(:,1);
candidate = R;
if mode=="representative"
    candidate.valid(:) = false;
    if isfinite(selected), candidate.valid(selected,1) = true; end
    for j = 1:size(R.valid,2)
        for leg = 1:2
            trial = p95LegTrial(R,j,leg,R.valid(:,j));
            if isfinite(trial), candidate.valid(trial,j) = true; end
        end
    end
end

finiteMode = "all";
if mode=="none" || ~any(candidate.valid,'all'), finiteMode = "none"; end
candidate.cfg.finiteMaxCandidates = Inf;
V = FABER_RPO_finite_validation(candidate,finiteMode);
V.nominalValidCount = nnz(R.valid);
V.attemptedCount = nnz(V.attempted);
V.passedCount = nnz(V.passed);
V.attemptedFraction = V.attemptedCount/max(1,V.nominalValidCount);

available = "Open";
finiteStatus = "Not evaluated";
finiteReason = "";
phases = table();
if isfinite(selected)
    available = "Pass";
    if V.attempted(selected,1)
        finiteStatus = "Pass";
        finiteReason = V.failureReason(selected,1);
        if ~V.passed(selected,1), finiteStatus = "Fail"; end
    end
    phases = phaseChecks(R.records{selected,1},R.phaseNames,cfg,V,selected);
end
if mode=="none", finiteStatus = "Disabled"; end

names = ["Monte Carlo trials";"Nominal client approaches passed";"Nominal paired trials passed"; ...
    "Budget representative trial";"Finite-thrust candidates replayed";"Finite-thrust candidates passed"];
values = [numel(paired);nnz(~cellfun(@isempty,R.clientCandidates(:,1)));nnz(paired); ...
    selected;nnz(V.attempted);nnz(V.passed)];
counts = table(names,values,'VariableNames',{'Metric','Value'});

test = ["Lambert perigee and endpoint";"Relative speed and transfer geometry"; ...
    "Approach zone, keep-out zone and target passage";"Directional thrust and burn duration"; ...
    "Impulsive burn fraction";"Braking distance";"Client obstacle clearance"; ...
    "Atlante obstacle clearance";"Terminal corridor";"Peak acceleration cap"; ...
    "Cycle propellant cap";"Cycle delta-v cap";"Finite-thrust replay"];
status = repmat(available,numel(test),1);
notes = strings(numel(test),1);
notes(1) = "Outbound and return arcs for the paired case";
notes(2) = "Nominal nine-phase screening";
notes(4) = "Phase-specific main/RCS propulsion";
notes(5) = sprintf('Maximum %g%% per powered phase',100*cfg.maxBurnFraction);
notes(6) = "Nominal radial stopping-distance screen";
optional = [hasGeometry(cfg.clientGeometry);hasGeometry(cfg.hostGeometry); ...
    isfinite(cfg.terminalCorridorHalfWidth_m);isfinite(cfg.maxPeakAccel_mps2); ...
    isfinite(cfg.maxPropellantCycle_kg);isfinite(cfg.maxDV_perCycle_mps)];
status(find(~optional)+6) = "Not configured";
notes(9) = "LVLH terminal lateral half-width";
status(end) = finiteStatus;
notes(end) = finiteReason;
if ~cfg.checkBrakingDistance, status(6) = "Disabled"; end

failures = string(R.rejectionReason(:,1));
failures = failures(failures~="");
[reasons,~,group] = unique(failures);
trials = zeros(numel(reasons),1);
if ~isempty(group), trials = accumarray(group,1); end
rejections = table(reasons,trials, ...
    'VariableNames',{'FirstRejection','Trials'});
rejections = sortrows(rejections,'Trials','descend');

[cases,casePhases] = caseTables(R,V,selected);
architectures = R.summary;
architectures.FiniteAttempted = sum(V.attempted,1)';
architectures.FinitePassed = sum(V.passed,1)';
architectures.P95CycleDVFinitePass_mps = nan(height(architectures),1);
for j = 1:height(architectures)
    dv = sort(R.CycleDV_mps(V.passed(:,j),j));
    if ~isempty(dv)
        architectures.P95CycleDVFinitePass_mps(j) = dv(ceil(0.95*numel(dv)));
    end
end
checks = struct('summary',table(test,status,notes,'VariableNames',{'Check','Status','Notes'}), ...
    'counts',counts,'phases',phases,'rejections',rejections,'finite',V, ...
    'finiteStatus',finiteStatus,'mode',mode,'cases',cases,'casePhases',casePhases, ...
    'architectures',architectures,'propulsion',propulsionRows(R,V,mode));
end

function trial = p95LegTrial(R,architecture,leg,eligible)
ids = find(eligible);
trial = NaN;
if ~isempty(ids)
    [~,order] = sort(R.DV_mps(ids,architecture,leg));
    trial = ids(order(ceil(0.95*numel(ids))));
end
end

function rows = propulsionRows(R,V,mode)
data = cell(0,15);
legs = ["Client","Atlante"];
cfg = R.cfg;
for j = 1:size(R.valid,2)
    populations = "Nominal leg P95";
    pools = {R.valid(:,j)};
    if mode=="all" && any(V.passed(:,j))
        populations(end+1) = "Finite-pass leg P95";
        pools{end+1} = V.passed(:,j);
    end
    for i = 1:numel(populations)
        for leg = 1:2
            trial = p95LegTrial(R,j,leg,pools{i});
            if ~isfinite(trial), continue; end
            rec = R.records{trial,j};
            approaches = {rec.Client,rec.ATLANTE};
            a = approaches{leg};
            for phase = 1:2:9
                k = (phase+1)/2;
                p = FABER_RPO_phase_propulsion(cfg,phase);
                type = "Legacy";
                installed = NaN;
                if isfield(cfg,'propulsionMode') && strcmpi(cfg.propulsionMode,'phase')
                    system = cfg.phasePropulsion{phase};
                    type = upper(string(system));
                    installed = cfg.propulsion.(system).nInstalled;
                end
                firing = NaN;
                status = "Not replayed";
                reason = "";
                if V.attempted(trial,j)
                    reason = V.phaseReason{trial,j}(leg,k);
                    status = "Pass";
                    if reason~="", status = "Fail"; end
                    finite = V.details{trial,j}{leg,k};
                    if ~isempty(finite), firing = finite.burnTime_s; end
                end
                nominal = a.burnTimePhase_s(phase);
                if ~a.transferDetails{k}.isImpulsive, nominal = NaN; end
                data(end+1,:) = {populations(i),legs(leg),j,trial,string(R.phaseNames{phase}), ...
                    a.dvPhase(phase),type,installed,p.nActiveThrusters,p.thrust_N,p.Isp_s, ...
                    nominal,firing,status,reason}; %#ok<AGROW>
            end
        end
    end
end
rows = cell2table(data,'VariableNames',{'Population','Leg','Architecture','Trial','Phase', ...
    'RawDV_mps','ThrusterType','InstalledThrusters','ConfiguredActiveThrusters', ...
    'RatedThrustEach_N','Isp_s','NominalImpulsiveFiring_s','FiniteFiring_s','FiniteStatus','FiniteReason'});
end

function t = phaseChecks(record,phaseNames,cfg,V,selected)
rows = cell(0,17);
approaches = {record.Client,record.ATLANTE};
legs = ["Client","Atlante"];
geometry = {cfg.clientGeometry,cfg.hostGeometry};
limits = [cfg.maxSpeed_homing1_mps,cfg.maxSpeed_homing2_mps,cfg.maxSpeed_closing1_mps, ...
    cfg.maxSpeed_closing2_mps,cfg.maxSpeed_final_mps];
for leg = 1:2
    a = approaches{leg};
    for phase = 1:9
        x = a.trace{phase};
        v = x(4:6,:);
        velocity = max(vecnorm(v));
        radius = vecnorm(x(1:3,:));
        closing = max([0,-sum(x(1:3,:).*v,1)./max(radius,eps)]);
        speedLimit = NaN;
        clearance = NaN;
        finiteStatus = "Not applicable";
        detail = "";
        err = [NaN NaN];
        if mod(phase,2)==1
            k = (phase+1)/2;
            speedLimit = limits(k);
            if hasGeometry(geometry{leg})
                clearance = FABER_path_clearance(x(1:3,:),geometry{leg});
            end
            finiteStatus = "Not evaluated";
            if V.attempted(selected,1)
                detail = V.phaseReason{selected,1}(leg,k);
                finiteStatus = "Pass";
                if detail~="", finiteStatus = "Fail"; end
                p = V.details{selected,1}{leg,k};
                if ~isempty(p), err = [p.positionError_m,p.velocityError_mps]; end
            end
        end
        rows(end+1,:) = {legs(leg),string(phaseNames{phase}),a.dvPhase(phase), ...
            a.tofPhase(phase),velocity,closing,speedLimit,min(radius),max(radius),clearance, ...
            a.burnTimePhase_s(phase),a.burnFractionPhase(phase),a.fuelPhase_kg(phase), ...
            finiteStatus,detail,err(1),err(2)}; %#ok<AGROW>
    end
end
t = cell2table(rows,'VariableNames',{'Leg','Phase','RawDV_mps','TOF_min', ...
    'MaxRelativeSpeed_mps','MaxClosingSpeed_mps','SpeedLimit_mps','MinRange_m','MaxRange_m', ...
    'MinObstacleClearance_m','NominalBurnTime_s','NominalBurnFraction','NominalPropellant_kg', ...
    'FiniteStatus','FiniteReason','FinitePositionError_m','FiniteVelocityError_mps'});
end

function [cases,phases] = caseTables(R,V,selected)
% Preserve individual leg estimates even when the paired trial is rejected.
[n,m] = size(R.valid);
[trial,architecture] = ndgrid(1:n,1:m);
methods = R.methods(architecture(:),:);
dv = R.DV_mps;
tof = R.TOF_min;
phaseDV = R.DV_phase_mps;
phaseTOF = R.TOF_phase_min;
candidates = {R.clientCandidates,cell(n,m)};
if isfield(R,'hostCandidates'), candidates{2} = R.hostCandidates; end
for leg = 1:2
    entries = candidates{leg};
    [ii,mm] = find(~R.valid & ~cellfun(@isempty,entries));
    for j = 1:numel(ii)
        a = entries{ii(j),mm(j)};
        dv(ii(j),mm(j),leg) = a.dv;
        tof(ii(j),mm(j),leg) = a.tof;
        phaseDV(ii(j),mm(j),leg,:) = reshape(a.dvPhase,1,1,1,[]);
        phaseTOF(ii(j),mm(j),leg,:) = reshape(a.tofPhase,1,1,1,[]);
    end
end

finite = repmat("Not replayed",n,m);
finite(~R.valid) = "Not eligible";
finite(V.attempted) = "Fail";
finite(V.passed) = "Pass";
accepted = repmat("Paired pass",n,m);
accepted(~R.valid) = "Rejected";
accepted(~R.valid & isfinite(dv(:,:,1)) & ~isfinite(dv(:,:,2))) = "Client only";
cycle = sum(dv,3);
totalTOF = sum(tof,3);
client = dv(:,:,1);
host = dv(:,:,2);
clientTime = tof(:,:,1);
hostTime = tof(:,:,2);
cases = table(trial(:),architecture(:),methods(:,1),methods(:,2),methods(:,3),methods(:,4), ...
    trial(:)==selected & architecture(:)==1,accepted(:),client(:),2*client(:),host(:),2*host(:),cycle(:),2*cycle(:), ...
    clientTime(:),hostTime(:),totalTOF(:),string(R.rejectionReason(:)),finite(:),V.failureReason(:), ...
    'VariableNames',{'Trial','Architecture','Homing2','Closing1','Closing2','GateMode', ...
    'BudgetRepresentative','NominalStatus','ClientRawDV_mps','ClientMarginedDV_mps', ...
    'AtlanteRawDV_mps','AtlanteMarginedDV_mps', ...
    'CycleRawDV_mps','CycleMarginedDV_mps','ClientTOF_min','AtlanteTOF_min','CycleTOF_min', ...
    'NominalRejection','FiniteStatus','FiniteReason'});
cases = sortrows(cases,{'Trial','Architecture'});

[t,a,l,k] = ind2sub(size(phaseDV),find(isfinite(phaseDV)));
ids = sub2ind(size(phaseDV),t,a,l,k);
legs = ["Client";"Atlante"];
labels = reshape(string(R.phaseNames),[],1);
phases = table(t,a,legs(l),k,labels(k),phaseDV(ids),2*phaseDV(ids),phaseTOF(ids), ...
    'VariableNames',{'Trial','Architecture','Leg','PhaseNumber','Phase', ...
    'RawDV_mps','MarginedDV_mps','TOF_min'});
phases = sortrows(phases,{'Trial','Architecture','Leg','PhaseNumber'});
end

function configured = hasGeometry(geometry)
configured = ~isempty(fieldnames(geometry));
end
