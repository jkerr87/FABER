function out = main_FABER_round_trip_safety(options)
% Coupled Taxi round trip dispersions, directed faults and optional RPO CAM screening.
% The Phase A budget and this study share the same nominal RPO and Lambert code.
if nargin<1
    options = struct();
end

root = fileparts(mfilename('fullpath'));
addpath(fullfile(root,'ORBMECH'),fullfile(root,'ProximityOps'), ...
    fullfile(root,'MonteCarlo','FABER_safety_v3_nine_phase'),'-begin');
cfg = studyOptions(options);
phase = FABER_phaseA_defaults();
phase.captureWait_min = 45;
phase = override(phase,cfg.phaseA);
phase.N_RPO_MC = 1;
validateattributes(phase.captureWait_min,{'numeric'},{'scalar','finite','nonnegative'});
validateattributes(phase.farFieldIsp_s,{'numeric'},{'scalar','finite','positive'});
validateattributes(phase.massTaxi_kg,{'numeric'},{'scalar','finite','positive'});
validateattributes(phase.payloadMass_kg,{'numeric'},{'scalar','finite','positive'});
if exist('astroConstants','file')==2
    mu = astroConstants(13);
    Re = astroConstants(23);
else
    mu = 398600.4418;
    Re = 6378.1363;
end

% Solve the reference transfer once and reuse its outbound commands across trials.
nominal = FABER_outbound_lambert(mu,Re,phase);
randomState = rng;
restoreRng = onCleanup(@() rng(randomState)); %#ok<NASGU>
stream = RandStream('mt19937ar','Seed',cfg.seed);

% Reuse the same random draws when comparing cases with dispersions enabled.
draw = sampleInputs(stream,cfg);
zeroDraw = draw(1,:);
for name = string(draw.Properties.VariableNames)
    zeroDraw.(name) = zeros(size(zeroDraw.(name)));
end

labels = cfg.cases(:)';
total = cfg.N*numel(labels);
rows = cell(total,1);
records = cell(total,1);
clientEvents = [];
hostEvents = [];
nominalClientEvents = [];
nominalHostEvents = [];
returnPlan = []; returnArc = []; captureDrift = []; returnReason = "";
trial = 0;

for label = labels
    for i = 1:cfg.N
        trial = trial+1;
        d = draw(i,:);
        if label=="Nominal" || (label~="Dispersed" && ~cfg.directedUsesDispersions)
            d = zeroDraw;
        end

        returnPlan = []; returnArc = []; captureDrift = []; returnReason = "";

        % Apply timing and state errors before propagating the outbound coast.
        depEpoch_s = nominal.departureEpoch_s+d.DepartureTime_s;
        baseTaxi = coast(nominal.departurePreBurnECI,d.DepartureTime_s);
        baseModule = coast(nominal.moduleDepartureECI,d.DepartureTime_s);
        atlante = coast(nominal.atlanteDepartureECI,d.DepartureTime_s);
        taxi = perturb(baseTaxi,d.InjectionLVLH');
        module = perturb(baseModule,d.ModuleLVLH');
        eOut = execution(d,1,label,"Outbound",cfg.partialBurnFraction);
        arc = executeArc(nominal,taxi,module,phase.x0_client,eOut);
        arc.departureEpoch_s = depEpoch_s;
        arc.arrivalEpoch_s = depEpoch_s+arc.duration_s;
        client = []; host = [];
        status = "Outbound far field";
        reason = arc.reason;

        % Only an acceptable achieved gate state can initialize client RPO.
        if arc.accepted
            handoff = nominal;
            handoff.arrivalLVLH = arc.arrivalLVLH;
            handoff.departureEpoch_s = depEpoch_s;
            handoff.arrivalEpoch_s = arc.arrivalEpoch_s;
            handoff.moduleArrivalECI = arc.targetArrivalECI;
            handoff.atlanteDepartureECI = atlante;

            handoff.departureBurn_mps = arc.departureBurn_mps;
            handoff.arrivalBurn_mps = arc.arrivalBurn_mps;
            handoff.truthCoast = @coast;
            handoff.planEphemeris = @(x,dt) ephemeris(mu,Re,x,dt,0);

            mc = FABER_phaseA_rpo_config(phase,mu,Re,handoff);
            mc.seed = cfg.rpoSeed+i-1;
            mc = override(mc,cfg.rpo);
            mc.safetySaveBurns = cfg.rpoFaultTrials>0 && ismember(label,["Dispersed","Nominal"]);
            mc.returnHandoff = @coupledReturn;
            R = [];
            evalc('R = FABER_RPO_9phase_MC(mc);');
            client = R.clientCandidates{1,1};
            if R.valid(1,1)
                record = R.records{1,1};
                host = record.ATLANTE;
                status = "Complete";
                reason = "";
            elseif isempty(client)
                status = "Client RPO";
                reason = R.rejectionReason(1,1);
            elseif strlength(returnReason)>0
                status = "Return far field";
                if label=="Capture failure"
                    status = "Capture";
                end

                reason = returnReason;
            else
                status = "Atlante RPO";
                reason = R.rejectionReason(1,1);
            end

            if mc.safetySaveBurns && label=="Dispersed"
                if ~isempty(client)
                    clientEvents = [clientEvents,client.safetyEvents]; %#ok<AGROW>
                end

                if ~isempty(host)
                    hostEvents = [hostEvents,host.safetyEvents]; %#ok<AGROW>
                end
            elseif mc.safetySaveBurns
                if ~isempty(client)
                    nominalClientEvents = [nominalClientEvents,client.safetyEvents]; %#ok<AGROW>
                end

                if ~isempty(host)
                    nominalHostEvents = [nominalHostEvents,host.safetyEvents]; %#ok<AGROW>
                end
            end
        end

        records{trial} = struct('outbound',arc,'client',client,'captureDrift',captureDrift, ...
            'returnPlan',returnPlan,'returnArc',returnArc,'host',host);
        rows{trial} = trialRow(trial,label,i,d,status,reason,records{trial});
    end
end

% Preserve full trial records alongside tabular gate and clearance metrics.
out = struct();
out.options = cfg;
out.phaseA = phase;
out.nominalOutbound = nominal;
out.trials = vertcat(rows{:});
out.records = records;
trials = out.trials;
farMin = min([trials.OutboundClosestApproach_m,trials.ReturnClosestApproach_m],[],2,'omitnan');
out.trials.FarProtectedEnvelopeEntered = farMin<cfg.farProtectedRadius_m;
out.trials.FarContactRadiusReached = farMin<=cfg.contactRadius_m;
out.trials.PassiveContactRadiusReached = trials.PassiveClosestApproach_m<=cfg.contactRadius_m;
out.summary = summarize(out.trials);
out.sensitivity = sensitivity(out.trials);
out.rpoSafetyClient = rpoSafety(clientEvents,nominalClientEvents, ...
    phase.h_target_km+phase.moduleAltitudeOffset_km,"clientGeometry");
out.rpoSafetyAtlante = rpoSafety(hostEvents,nominalHostEvents,phase.h_target_km,"hostGeometry");
out.figures = gobjects(0);
if cfg.makePlots
    out.figures = FABER_plot_round_trip_safety(out);
end

if strlength(cfg.outputFolder)>0
    saveOutputs(out);
end

disp(out.summary);

    function x = coast(x0,dt)
        atTime = ephemeris(mu,Re,x0,dt,cfg.truthJ2);
        x = atTime(dt);
    end

    function [transfer,why] = coupledReturn(clientState,massCaptured)
    % After client RPO, target the return from the captured stack state.
        transfer = [];
        if label=="Capture failure"
            chief = coast(handoff.moduleArrivalECI,clientState.elapsed_s);
            chaser = FABER_relative_state(chief,clientState.terminalLVLH,'lvlhToEci');
            captureDrift = passiveCoast(chaser,chief);
            returnReason = "Capture failed. No payload attached or return attempted.";
        else
            [returnPlan,returnReason] = FABER_return_lambert(mu,Re,phase,handoff, ...
                clientState,massCaptured);
            if ~isempty(returnPlan)
                eReturn = execution(d,3,label,"Return",cfg.partialBurnFraction);
                returnArc = executeArc(returnPlan,returnPlan.departurePreBurnECI, ...
                    returnPlan.atlanteDepartureECI,phase.x0_atlante,eReturn);
                returnArc.departureEpoch_s = returnPlan.departureEpoch_s;
                returnArc.arrivalEpoch_s = returnPlan.departureEpoch_s+returnArc.duration_s;
                returnArc.massArrival_kg = massCaptured*exp(-returnArc.dv_mps/ ...
                    (9.80665*phase.farFieldIsp_s));
                returnReason = returnArc.reason;
                if returnArc.accepted
                    transfer = returnPlan;
                    transfer.arrivalLVLH = returnArc.arrivalLVLH;
                    transfer.arrivalPostBurnECI = returnArc.arrivalECI;
                    transfer.referenceArrivalECI = returnArc.targetArrivalECI;
                    transfer.arrivalBurn_mps = returnArc.arrivalBurn_mps;
                    transfer.departureBurn_mps = returnArc.departureBurn_mps;
                    transfer.massArrival_kg = returnArc.massArrival_kg;
                end
            end
        end

        why = returnReason;
    end

    function result = executeArc(plan,start,chief,gate,errors)
        departure = errors.departureScale*plan.departureBurn_mps + errors.departureError;
        arrival = errors.arrivalScale*plan.arrivalBurn_mps + errors.arrivalError;
        if errors.departureMissed
            departure = zeros(3,1);
            arrival = zeros(3,1);
        elseif errors.arrivalMissed
            arrival = zeros(3,1);
        end

        duration = plan.transferTime_s;
        start(4:6) = start(4:6)+departure;
        chaserAt = ephemeris(mu,Re,start,duration,cfg.truthJ2);
        chiefAt = ephemeris(mu,Re,chief,duration,cfg.truthJ2);
        result = coastMetrics(chaserAt,chiefAt,duration,Re,cfg);
        target = chiefAt(duration);
        arrivalPreBurn = chaserAt(duration);
        arrival = arrival(:);
        % A protected-envelope or Earth-altitude breach stops the arrival command.
        if result.minRadius_m < cfg.farProtectedRadius_m || ...
                result.minAltitude_km < phase.outboundMinimumPerigee_km
            arrival = zeros(3,1);
        end

        afterBurn = arrivalPreBurn;
        afterBurn(4:6) = afterBurn(4:6)+arrival;
        relative = FABER_relative_state(target,afterBurn,'eciToLvlh');
        residual = relative-gate(:);
        result.positionError_m = norm(residual(1:3));
        result.velocityError_mps = norm(residual(4:6));
        result.arrivalLVLH = relative;
        result.arrivalECI = afterBurn;
        result.targetArrivalECI = target;

        result.departureBurn_mps = departure;
        result.arrivalBurn_mps = arrival;
        result.dv_mps = norm(departure)+norm(arrival);
        result.duration_s = duration;
        result.reason = "";
        if result.minAltitude_km < phase.outboundMinimumPerigee_km
            result.reason = "Far-field altitude below limit";
        elseif result.minRadius_m < cfg.farProtectedRadius_m
            result.reason = "Far-field protected-envelope entry";
        elseif errors.departureMissed || errors.arrivalMissed
            result.reason = "Commanded far-field burn missed";
        elseif result.positionError_m>cfg.gatePositionTolerance_m || ...
                result.velocityError_mps>cfg.gateVelocityTolerance_mps
            result.reason = "RPO entry outside gate tolerances";
        end

        result.accepted = strlength(result.reason)==0;
        result.passiveDrift = [];
        if ~result.accepted && result.minAltitude_km>=phase.outboundMinimumPerigee_km ...
                && cfg.passiveHorizon_s>0
            result.passiveDrift = passiveCoast(afterBurn,target);
        end
    end

    function result = passiveCoast(chaser,chief)
    % Measure post-fault drift over the configured observation horizon.
        horizon = cfg.passiveHorizon_s;
        chaserAt = ephemeris(mu,Re,chaser,horizon,cfg.truthJ2);
        chiefAt = ephemeris(mu,Re,chief,horizon,cfg.truthJ2);
        result = coastMetrics(chaserAt,chiefAt,horizon,Re,cfg);
    end

    function result = rpoSafety(events,nominalEvents,h_km,geometryField)
    % Use dispersed events when available and fall back to nominal events.
        result = [];
        source = "Dispersed";
        if isempty(events)
            events = nominalEvents;
            source = "Nominal";
        end
        if isempty(events)
            return;
        end

        settings = cfg.rpoSafety;
        settings.N = cfg.rpoFaultTrials;
        settings.vehicleMass_kg = phase.massTaxi_kg;
        hasGroups = isfield(settings,'camForceVectors_LVLH_N') && ...
            ~isempty(settings.camForceVectors_LVLH_N);
        hasVehicleGroups = isfield(settings,'camForceVectorsByVehicle') && ...
            isfield(settings.camForceVectorsByVehicle,'Taxi') && ...
            ~isempty(settings.camForceVectorsByVehicle.Taxi);
        if ~isfield(settings,'enableActiveCAM')
            settings.enableActiveCAM = hasGroups || hasVehicleGroups;
        end

        assert(~settings.enableActiveCAM || hasGroups || hasVehicleGroups, ...
            'Active CAM screening requires verified Taxi firing groups.');
        if ~isfield(settings,'hazards') && isfield(cfg.rpo,geometryField) && ...
                ~isempty(fieldnames(cfg.rpo.(geometryField)))
            settings.hazards = cfg.rpo.(geometryField);
        end

        result = FABER_safety_MC_taxi(events,mu,Re,h_km,settings);
        result.summary.SourceCase = source;
    end
end

function cfg = studyOptions(options)
supported = ["Nominal","Dispersed","Outbound departure missed","Outbound arrival missed", ...
    "Outbound departure partial","Capture failure","Return departure missed", ...
    "Return arrival missed","Return departure partial"];
cfg = struct('N',30,'seed',4242,'rpoSeed',42,'truthJ2',0, ...
    'sigmaInjectionLVLH',[10 10 5 .001 .001 .0005], ...
    'sigmaModuleLVLH',[20 50 10 .001 .001 .0005], ...
    'sigmaBurnScale',0.00005,'sigmaBurnVector_mps',0.0005, ...
    'sigmaDepartureTime_s',0.1,'partialBurnFraction',0.8,'directedUsesDispersions',false, ...
    'gatePositionTolerance_m',500,'gateVelocityTolerance_mps',0.1, ...
    'farProtectedRadius_m',300,'contactRadius_m',2,'passiveHorizon_s',7200, ...
    'sampleStep_s',60,'makePlots',true,'plotLimit',40,'outputFolder',"", ...
    'phaseA',struct(),'rpo',struct(),'rpoFaultTrials',30,'rpoSafety',struct());
cfg.cases = supported([1 2 3 4 6 7 8]);
cfg = override(cfg,options);
cfg.cases = unique(string(cfg.cases),'stable');
assert(~isempty(cfg.cases) && all(ismember(cfg.cases,supported),'all'), ...
    'Unknown or empty round trip fault case list.');
validateattributes(cfg.directedUsesDispersions,{'logical'},{'scalar'});
validateattributes(cfg.N,{'numeric'},{'scalar','integer','positive','finite'});
validateattributes(cfg.rpoFaultTrials,{'numeric'},{'scalar','integer','nonnegative','finite'});
validateattributes(cfg.plotLimit,{'numeric'},{'scalar','integer','positive','finite'});
validateattributes(cfg.sigmaInjectionLVLH,{'numeric'},{'vector','numel',6,'finite','nonnegative'});
validateattributes(cfg.sigmaModuleLVLH,{'numeric'},{'vector','numel',6,'finite','nonnegative'});
for field = ["sigmaBurnScale","sigmaBurnVector_mps","sigmaDepartureTime_s","truthJ2","passiveHorizon_s"]
    validateattributes(cfg.(field),{'numeric'},{'scalar','real','finite','nonnegative'});
end

for field = ["sampleStep_s","gatePositionTolerance_m","gateVelocityTolerance_mps","contactRadius_m"]
    validateattributes(cfg.(field),{'numeric'},{'scalar','real','finite','positive'});
end

validateattributes(cfg.farProtectedRadius_m,{'numeric'}, ...
    {'scalar','real','finite','>',cfg.contactRadius_m});
validateattributes(cfg.partialBurnFraction,{'numeric'},{'scalar','real','finite','>',0,'<',1});
% Preserve the one-candidate paired interface and the actual propagated state.
protected = ["N","seed","x0_client","x0_host","arrivalStatesClient","arrivalStatesHost", ...
    "arrivalEpochClient_s","clientH_km","h_km","mu_km3s2","Re_km","returnHandoff", ...
    "massTaxi_kg","payloadMass_kg","enablePayloadMassChange","methodsH2","methodsC1", ...
    "methodsC2","gateModes","makePlots","saveCSV"];
assert(~any(ismember(string(fieldnames(cfg.rpo)),protected)), ...
    'Set orbit and timing in phaseA. RPO overrides are for feasibility inputs.');
end

function cfg = override(cfg,values)
for key = string(fieldnames(values))'
    cfg.(key) = values.(key);
end
end

function d = sampleInputs(stream,cfg)
% Draw initial-state errors and separate execution errors for each far-field burn.
N = cfg.N;
injection = randn(stream,N,6).*cfg.sigmaInjectionLVLH(:)';
module = randn(stream,N,6).*cfg.sigmaModuleLVLH(:)';
scales = cfg.sigmaBurnScale*randn(stream,N,4);
vectors = cfg.sigmaBurnVector_mps*randn(stream,N,12);
time = cfg.sigmaDepartureTime_s*randn(stream,N,1);
d = table(injection,module,scales,vectors,time,'VariableNames', ...
    {'InjectionLVLH','ModuleLVLH','BurnScaleError','BurnVectorError_mps','DepartureTime_s'});
end

function e = execution(d,index,label,leg,partial)
e = struct('departureScale',max(0,1+d.BurnScaleError(index)), ...
    'arrivalScale',max(0,1+d.BurnScaleError(index+1)), ...
    'departureError',d.BurnVectorError_mps(3*index-2:3*index)', ...
    'arrivalError',d.BurnVectorError_mps(3*index+1:3*index+3)', ...
    'departureMissed',label==leg+" departure missed", ...
    'arrivalMissed',label==leg+" arrival missed");
if label==leg+" departure partial"
    e.departureScale = partial*e.departureScale;
end
end

function x = perturb(base,delta)
x = FABER_relative_state(base,delta,'lvlhToEci');
end

function atTime = ephemeris(mu,Re,x,dt,J2)
if dt==0
    atTime = @(t) repmat(x(:),1,numel(t));
    return;
end

opts = odeset('RelTol',2e-12,'AbsTol',1e-9);
sol = ode113(@(~,y) gravity(y,1e9*mu,1000*Re,J2),[0 dt],x(:),opts);
atTime = @(t) deval(sol,t);
end

function dx = gravity(x,mu,Re,J2)
p = x(1:3);
r = norm(p);
a = -mu*p/r^3;
if J2~=0
    z2r2 = (p(3)/r)^2;
    a = a + (1.5*J2*mu*Re^2/r^5)*p.*([5;5;5]*z2r2-[1;1;3]);
end

dx = [x(4:6);a];
end

function result = coastMetrics(chaserAt,chiefAt,duration,Re,cfg)
t = linspace(0,duration,max(2,ceil(duration/cfg.sampleStep_s)+1));
chaser = chaserAt(t);
chief = chiefAt(t);
distance = vecnorm(chaser(1:3,:)-chief(1:3,:));
altitude = vecnorm(chaser(1:3,:))/1000-Re;
[minRadius,tMin] = refinedMin(t,distance,@(q) positionDistance(chaserAt(q),chiefAt(q)));
minAltitude = refinedMin(t,altitude,@(q) positionDistance(chaserAt(q),zeros(6,1))/1000-Re);
% Keep a thinner trajectory for plotting while screening on the full time grid.
keep = unique(round(linspace(1,numel(t),min(numel(t),500))));
plotT = sort(unique([t(keep),tMin]));
chaser = chaserAt(plotT);
chief = chiefAt(plotT);
relative = zeros(6,numel(plotT));
for k = 1:numel(plotT)
    relative(:,k) = FABER_relative_state(chief(:,k),chaser(:,k),'eciToLvlh');
end

result = struct('time_s',plotT,'trajectoryLVLH',relative, ...
    'minRadius_m',minRadius,'timeOfClosestApproach_s',tMin,'minAltitude_km',minAltitude);
end

function distance = positionDistance(x,y)
distance = norm(x(1:3)-y(1:3));
end

function [best,when] = refinedMin(t,x,objective)
% Refine sampled minima to estimate closest approach and minimum altitude.
[best,index] = min(x);
when = t(index);
local = find(x(2:end-1)<=x(1:end-2) & x(2:end-1)<=x(3:end))+1;
for k = unique([1 local numel(t)])
    a = t(max(1,k-1));
    b = t(min(numel(t),k+1));
    if b>a
        [tk,xk] = fminbnd(objective,a,b);
        if xk<best
            best = xk;
            when = tk;
        end
    end
end
end

function row = trialRow(trial,label,drawId,d,status,reason,result)
% Keep missing values for downstream stages the trial did not reach.
outbound = result.outbound;
client = result.client;
plan = result.returnPlan;
inbound = result.returnArc;
host = result.host;
captureEpoch = NaN; returnDeparture = NaN; returnArrival = NaN;
captureMass = NaN; arrivalMass = NaN; returnDV = NaN;
returnPos = NaN; returnVel = NaN; returnMin = NaN; returnAlt = NaN;
clientDV = NaN; hostDV = NaN; totalDV = NaN; elapsed_hr = NaN;
driftDistances = [];
if ~isempty(outbound.passiveDrift)
    driftDistances(end+1) = outbound.passiveDrift.minRadius_m;
end

if ~isempty(result.captureDrift)
    driftDistances(end+1) = result.captureDrift.minRadius_m;
end

if ~isempty(client)
    captureEpoch = outbound.arrivalEpoch_s+client.elapsed_s;
    clientDV = client.dv;
end

if ~isempty(plan)
    returnDeparture = plan.departureEpoch_s;
    captureMass = plan.massCaptured_kg;
end

if ~isempty(inbound)
    returnArrival = inbound.arrivalEpoch_s;
    returnDV = inbound.dv_mps;
    returnPos = inbound.positionError_m;
    returnVel = inbound.velocityError_mps;
    returnMin = inbound.minRadius_m;
    returnAlt = inbound.minAltitude_km;
    arrivalMass = inbound.massArrival_kg;
    if ~isempty(inbound.passiveDrift)
        driftDistances(end+1) = inbound.passiveDrift.minRadius_m;
    end
end

if ~isempty(host)
    hostDV = host.dv;
    totalDV = outbound.dv_mps+clientDV+returnDV+hostDV;
    elapsed_hr = (returnArrival+host.elapsed_s-outbound.departureEpoch_s)/3600;
end

passiveMin = NaN;
if ~isempty(driftDistances)
    passiveMin = min(driftDistances);
end

row = table(trial,label,drawId,status,reason,status=="Complete", ...
    outbound.accepted,~isempty(client),~isempty(plan) || ...
    (~isempty(client) && label~="Capture failure"),~isempty(inbound) && inbound.accepted, ...
    outbound.departureEpoch_s,outbound.arrivalEpoch_s,captureEpoch,returnDeparture,returnArrival, ...
    outbound.positionError_m,outbound.velocityError_mps,returnPos,returnVel, ...
    outbound.minRadius_m,returnMin,passiveMin,outbound.minAltitude_km,returnAlt, ...
    captureMass,arrivalMass,outbound.dv_mps,clientDV,returnDV,hostDV,totalDV,elapsed_hr, ...
    'VariableNames',{'Trial','Case','Draw','Stage','Reason','Complete', ...
    'OutboundGateAccepted','ClientRPOAccepted','CaptureSucceeded','ReturnGateAccepted', ...
    'OutboundDepartureEpoch_s','ClientArrivalEpoch_s','CaptureEpoch_s', ...
    'ReturnDepartureEpoch_s','AtlanteArrivalEpoch_s','ClientGatePositionError_m', ...
    'ClientGateVelocityError_mps','AtlanteGatePositionError_m','AtlanteGateVelocityError_mps', ...
    'OutboundClosestApproach_m','ReturnClosestApproach_m','PassiveClosestApproach_m', ...
    'OutboundMinAltitude_km','ReturnMinAltitude_km','CapturedMass_kg','AtlanteGateMass_kg', ...
    'OutboundFarDV_mps','ClientRPODV_mps','ReturnFarDV_mps','AtlanteRPODV_mps', ...
    'CompleteDV_mps','CompleteDuration_hr'});
row = [row,d];
end

function S = summarize(T)
% Summarize case outcomes and compute P95 costs for completed round trips.
S = table();
for label = unique(T.Case,'stable')'
    x = T(T.Case==label,:);
    n = height(x);
    completed = nnz(x.Complete);
    returnReached = nnz(isfinite(x.ReturnClosestApproach_m));
    passiveScreened = nnz(isfinite(x.PassiveClosestApproach_m));
    p = completed/n;
    z = 1.95996398454;
    center = (p+z^2/(2*n))/(1+z^2/n);
    halfWidth = z*sqrt(p*(1-p)/n+z^2/(4*n^2))/(1+z^2/n);
    S = [S; table(label,n,nnz(x.OutboundGateAccepted),nnz(x.ClientRPOAccepted), ...
        nnz(x.CaptureSucceeded),returnReached,nnz(x.ReturnGateAccepted),completed, ...
        p,max(0,center-halfWidth),min(1,center+halfWidth),passiveScreened, ...
        nnz(x.FarProtectedEnvelopeEntered),nnz(x.FarContactRadiusReached), ...
        nnz(x.PassiveContactRadiusReached), ...
        percentile95(x.CompleteDV_mps),percentile95(x.CompleteDuration_hr), ...
        'VariableNames',{'Case','Trials','OutboundGateAccepted','ClientRPOAccepted', ...
        'Captures','ReturnArcsExecuted','ReturnGateAccepted','Completed','CompletionFraction', ...
         'Wilson95Lower','Wilson95Upper','PassiveCoastsScreened', ...
        'FarProtectedEnvelopeEntries','FarContactRadiusCases','PassiveContactRadiusCases', ...
        'P95CompleteDV_mps', ...
        'P95CompleteDuration_hr'})]; %#ok<AGROW>
end
end

function p = percentile95(x)
x = x(isfinite(x));
p = NaN;
if ~isempty(x)
    p = prctile(x,95);
end
end

function S = sensitivity(T)
% Descriptive rank correlations, using dispersed trials only.
T = T(T.Case=="Dispersed",:);
inputs = ["Injection position";"Injection velocity";"Module position";"Module velocity"; ...
    "Departure epoch";"Outbound departure burn scale";"Outbound arrival burn scale"; ...
    "Return departure burn scale";"Return arrival burn scale"; ...
    "Outbound departure burn vector";"Outbound arrival burn vector"; ...
    "Return departure burn vector";"Return arrival burn vector"];
x = [vecnorm(T.InjectionLVLH(:,1:3),2,2),vecnorm(T.InjectionLVLH(:,4:6),2,2), ...
    vecnorm(T.ModuleLVLH(:,1:3),2,2),vecnorm(T.ModuleLVLH(:,4:6),2,2), ...
    abs(T.DepartureTime_s),abs(T.BurnScaleError), ...
    vecnorm(T.BurnVectorError_mps(:,1:3),2,2),vecnorm(T.BurnVectorError_mps(:,4:6),2,2), ...
    vecnorm(T.BurnVectorError_mps(:,7:9),2,2),vecnorm(T.BurnVectorError_mps(:,10:12),2,2)];
metrics = ["Client gate position error","Atlante gate position error","Complete round trip DV","Completion"];
y = [T.ClientGatePositionError_m,T.AtlanteGatePositionError_m,T.CompleteDV_mps,double(T.Complete)];
S = table();
for j=1:numel(metrics)
    for i=1:numel(inputs)
        keep = isfinite(x(:,i)) & isfinite(y(:,j));
        a = x(keep,i);
        b = y(keep,j);
        rho = NaN;
        if numel(a)>=3 && rangeLocal(a)>0 && rangeLocal(b)>0
            cc = corrcoef(ranks(a),ranks(b));
            rho = cc(1,2);
        end

        S = [S; table(inputs(i),metrics(j),nnz(keep),rho, ...
            'VariableNames',{'Input','Metric','Trials','SpearmanRho'})]; %#ok<AGROW>
    end
end
end

function w = rangeLocal(x)
w = max(x)-min(x);
end

function r = ranks(x)
[~,~,group] = unique(x);
[~,order] = sort(x);
r = accumarray(group(order),(1:numel(x))',[],@mean);
r = r(group);
end

function saveOutputs(out)
folder = out.options.outputFolder;
if ~isfolder(folder)
    mkdir(folder);
end

writetable(out.trials,fullfile(folder,'Round_Trip_Trials.csv'));
writetable(out.summary,fullfile(folder,'Round_Trip_Summary.csv'));
writetable(out.sensitivity,fullfile(folder,'Round_Trip_Sensitivity.csv'));
if ~isempty(out.rpoSafetyClient)
    writetable(out.rpoSafetyClient.trials,fullfile(folder,'Client_RPO_Faults.csv'));
    writetable(out.rpoSafetyClient.summary,fullfile(folder,'Client_RPO_Fault_Summary.csv'));
end

if ~isempty(out.rpoSafetyAtlante)
    writetable(out.rpoSafetyAtlante.trials,fullfile(folder,'Atlante_RPO_Faults.csv'));
    writetable(out.rpoSafetyAtlante.summary,fullfile(folder,'Atlante_RPO_Fault_Summary.csv'));
end

for k=1:numel(out.figures)
    f = out.figures(k);
    exportgraphics(f,fullfile(folder,string(f.Tag)+".png"),'Resolution',160);
end
end
