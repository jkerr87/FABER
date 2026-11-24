function plan = FABER_environment_schedule(budget,cfg,solar,mu,Re,enforce)
% Replay the selected Phase A paths, moving operations to consecutive sunlight.
% Times are seconds after the UTC epoch. Subphases use the original LVLH paths.
phase = budget.cfg;
baseline = budget.outboundLambert;
origin = baseline.departureEpoch_s;
captureActive_s = phase.captureWait_min*60;
if isfinite(cfg.captureActive_min), captureActive_s = cfg.captureActive_min*60; end
atlante = @(t) FABER_environment_circular(mu,baseline.atlanteDepartureECI,t);
module = @(t) FABER_environment_circular(mu,baseline.moduleDepartureECI,t);
park = @(t) fixedLVLH(atlante,t,phase.taxiParkingDepartureState);
events = table(strings(0,1),strings(0,1),strings(0,1),zeros(0,1),zeros(0,1), ...
    false(0,1),zeros(0,1),zeros(0,1),'VariableNames', ...
    {'Vehicle','Activity','Mode','Start_s','End_s','SunRequired','DeltaV_mps','Retrieval'});
links = table(strings(0,1),strings(0,1),zeros(0,1),zeros(0,1),zeros(0,1),strings(0,1), ...
    'VariableNames',{'Sender','Receiver','Start_s','End_s','Received_kW','Kind'});
issues = table(strings(0,1),strings(0,1),zeros(0,1),strings(0,1), ...
    'VariableNames',{'Vehicle','Activity','Epoch_s','Reason'});
holds = table(strings(0,1),strings(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    zeros(0,1),false(0,1),'VariableNames', ...
    {'Vehicle','BeforeActivity','Start_s','End_s','MaxDrift_m','MinClearance_m','Passed'});
paths = cell(0,1);
away = {zeros(0,2),zeros(0,2)};
recoveryDV = 0;
idx = budget.RPO_MC.SelectedCycle;
if isfinite(idx)
    record = budget.RPO_MC.records{idx,1};
else
    record = [];
end
explicit = ~isempty(cfg.retrievals);
retrievals = cfg.retrievals;
if ~explicit
    retrievals = struct('moduleId',"Module 1",'phaseAngle_deg',phase.phaseAngle_deg, ...
        'cargoCount',cfg.modulesPerRetrieval);
end

retrievalSummary = table((1:numel(retrievals))',reshape(string({retrievals.moduleId}),[],1), ...
    reshape([retrievals.phaseAngle_deg],[],1),reshape([retrievals.cargoCount],[],1), ...
    'VariableNames',{'Retrieval','Module','PhaseAngle_deg','CargoCount'});
for name = ["Start_s","Departure_s","Capture_s","Handover_s","StartMass_kg", ...
        "CapturedMass_kg","EndMass_kg","Propellant_kg","DeltaV_mps"]
    retrievalSummary.(name) = nan(height(retrievalSummary),1);
end
retrievalSummary.Status = repmat("Pending",height(retrievalSummary),1);
retrievalArcs = cell(numel(retrievals),2);
activeTrip = 0;
clientChief = module;
t = 0;
mass = phase.massTaxi_kg;
taxiDocked = false;
taxiBerth = park;
returnPlan = [];
moduleState = baseline.moduleDepartureECI;
moduleRate = sqrt(mu/norm(moduleState(1:3))^3);
for trip = 1:numel(retrievals)
    activeTrip = trip;
    initial = module((retrievals(trip).phaseAngle_deg-phase.phaseAngle_deg)*pi/180/moduleRate);
    clientChief = @(q) FABER_environment_circular(mu,initial,q);
    retrievalSummary.Start_s(trip) = t;
    retrievalSummary.StartMass_kg(trip) = mass;
    retrievalSummary.Status(trip) = "Blocked";

    departureState = park(t);
    if trip==1 && enforce && cfg.taxiInitialWait_s>0
        [t,departureState,ok] = waitAtParking(t,cfg.taxiInitialWait_s);
        if ~ok, break; end
    end
    if trip>1
        [t,mass,ok,departureState] = leaveBerth(t,mass,record.ATLANTE.terminalLVLH);
        if ~ok, break; end
        taxiDocked = false;
    end

    outbound = baseline;
    if t>0 || (explicit && retrievals(trip).phaseAngle_deg~=phase.phaseAngle_deg)
        outbound = targetModule(t,clientChief,departureState);
    end
    if isempty(outbound)
        block("Taxi","Outbound Lambert",t,"No arc passed the outbound perigee and endpoint checks.");
        break;
    end

    retrievalArcs{trip,1} = outbound;
    [t,endMass,berth,returnPlan] = retrieve(outbound,clientChief,t,mass,explicit);
    retrievalArcs{trip,2} = returnPlan;
    taxiDocked = isfinite(endMass);
    if ~taxiDocked, break; end

    taxiBerth = berth;
    retrievalSummary.Handover_s(trip) = t;
    retrievalSummary.EndMass_kg(trip) = endMass;
    retrievalSummary.Propellant_kg(trip) = retrievalSummary.StartMass_kg(trip)-endMass;
    retrievalSummary.Status(trip) = "Delivered";
    mass = endMass;
end
for trip = 1:height(retrievalSummary)
    if isfinite(retrievalSummary.Start_s(trip))
        retrievalSummary.DeltaV_mps(trip) = sum(events.DeltaV_mps(events.Retrieval==trip));
    end
end
activeTrip = 0;
taxiEnd = t;
delivered = retrievalSummary(retrievalSummary.Status=="Delivered",:);
allDelivered = height(delivered)==height(retrievalSummary);

% One Kiki/Rango pair works at a time. Every selected module receives its own
% Phase A loaded and unloaded path and explicit service operations.
indices = cfg.assemblyModuleIndices(:)';
completed = zeros(0,1);
safeStops = cell(0,1);
sorties = table('Size',[0 5],'VariableTypes',repmat("double",1,5), ...
    'VariableNames',{'Queue','ModuleIndex','Pair','LoadStart_s','BerthEnd_s'});
assemblyTime = 0;
deliveries = cfg.preloadedModules+cumsum(delivered.CargoCount);
for j=1:numel(indices)
    moduleIndex = indices(j);
    pair = cfg.assemblyPair(1+mod(j-1,numel(cfg.assemblyPair)));
    kiki = "Kiki "+pair; rango = "Rango "+pair;
    r = budget.assemblyRelocations.bestRecords{moduleIndex};
    if j>cfg.preloadedModules
        delivery = find(deliveries>=j,1);
        if isempty(delivery)
            block(kiki,"Assembly relocation "+moduleIndex,max(assemblyTime,taxiEnd), ...
                "Cargo has not been handed over.");
            break;
        end
        assemblyTime = max(assemblyTime,delivered.Handover_s(delivery));
    end
    if isempty(r)
        block(kiki,"Assembly relocation "+moduleIndex,assemblyTime,"Phase A found no feasible round trip.");
        break;
    end
    tag = "Assembly "+moduleIndex;
    home = @(q) fixedLVLH(atlante,q,[r.outbound.source_m(:);zeros(3,1)]);
    site = @(q) fixedLVLH(atlante,q,[r.outbound.target_m(:);zeros(3,1)]);
    if cfg.chargeAsNeeded
        [assemblyTime,ok] = prepareCharge(kiki,rango,tag,assemblyTime,home,cfg.assemblyCharge_min(j,:));
        if ~ok, break; end
    end
    loadStart = assemblyTime;
    [assemblyTime,ok] = operation(kiki,tag+" load cargo","K3",assemblyTime,cfg.kikiLoad_min*60,home,tag,false,true);
    if ~ok, break; end
    firstPath = numel(paths)+1;
    [assemblyTime,ok,departure] = reservedRelocation(kiki,tag+" loaded","K4l",r.outbound,assemblyTime);
    if ~ok
        safeStops{end+1,1} = struct('kiki',kiki,'rango',rango,'t',assemblyTime,'state',home,'tag',tag);
        break;
    end
    away{pair}(end+1,:) = [departure Inf];
    if cfg.retainRangoForRelocations
        carryRango(firstPath,rango,tag+" retained delivery");
    end
    if cfg.retainRangoForRelocations
        [assemblyTime,ok] = rangoTransfer(kiki,rango,tag+" deploy Rango",assemblyTime,site,tag);
        if ~ok, break; end
    end
    if ~cfg.chargeAsNeeded
        [assemblyTime,ok] = operation(kiki,tag+" charge Rango","K6",assemblyTime, ...
            cfg.rangoCharge_min*60,site,tag,false,true);
        if ~ok, break; end
        if cfg.rangoCharge_min>0
            chargeStart = assemblyTime-cfg.rangoCharge_min*60;
            append(rango,tag+" receive charge","R5",chargeStart,assemblyTime,false,site,tag,0);
            addLink(kiki,rango,chargeStart,assemblyTime,cfg.rangoPort_kW);
        end
    end
    handoverStart = assemblyTime;
    [assemblyTime,ok] = operation(kiki,tag+" handover to Rango","K5",assemblyTime, ...
        cfg.kikiHandover_min*60,site,tag,true,true);
    if ~ok, break; end
    if cfg.kikiHandover_min>0
        append(rango,tag+" receive cargo","R2",assemblyTime-cfg.kikiHandover_min*60, ...
            assemblyTime,true,site,tag,0);
    elseif assemblyTime>handoverStart
        error('Unexpected zero-duration handover delay.');
    end
    [workEnd,ok] = operation(rango,tag+" installation","R4",assemblyTime, ...
        cfg.rangoAssembly_min*60,site,tag,cfg.installationRequiresSun,true);
    if ~ok, break; end
    if cfg.retainRangoForRelocations
        append(kiki,tag+" attached standby during installation","K11",assemblyTime,workEnd,false,site,tag,0);
        [assemblyTime,ok] = rangoTransfer(kiki,rango,tag+" recover Rango",workEnd,site,tag);
        if ~ok, break; end
    end
    % Without the retained pair assumption, Kiki can return during installation.
    firstPath = numel(paths)+1;
    [assemblyTime,ok] = reservedRelocation(kiki,tag+" unloaded","K4e",r.returnTrip,assemblyTime);
    if ~ok
        safeStops{end+1,1} = struct('kiki',kiki,'rango',rango,'t',assemblyTime,'state',site,'tag',tag);
        break;
    end
    if cfg.retainRangoForRelocations
        carryRango(firstPath,rango,tag+" retained return");
    end
    if cfg.chargeAsNeeded
        [assemblyTime,ok] = operation(kiki,tag+" connect at Atlante","K2",assemblyTime, ...
            cfg.assemblyBerth_min*60,home,tag,false,true);
        if ~ok, break; end
    end
    away{pair}(end,2) = assemblyTime;
    sorties(end+1,:) = {j,moduleIndex,pair,loadStart,assemblyTime};
    if ~cfg.chargeAsNeeded
        chargeStart = assemblyTime;
        [assemblyTime,ok] = operation(kiki,tag+" recharge at Atlante","K7",assemblyTime, ...
            cfg.kikiCharge_min*60,home,tag,false,true);
        if ~ok, break; end
        addLink("Atlante",kiki,chargeStart,assemblyTime,cfg.kikiPort_kW);
    end
    assemblyTime = max(assemblyTime,workEnd);
    completed(end+1,1) = moduleIndex; %#ok<AGROW>
end
stop = max([taxiEnd,assemblyTime,events.End_s']);
if cfg.includeRefueling
    if allDelivered && numel(completed)==numel(indices)
        t = stop;
        dock = @(q) fixedLVLH(atlante,q,[0;0;0;0;0;0]);
        berthTime = 5.42*60;
        [t,ok] = operation("Taxi","Refuel berth","T9b",t,berthTime,taxiBerth,"Service",true,true);
        if ok
            a = t; t = t+cfg.taxiRefuel_min*60-2*berthTime;
            append("Taxi","Taxi refuel","T9",a,t,false,taxiBerth,"Service",0);
            append("Atlante","Supply Taxi propellant","A4",a,t,false,dock,"Service",0);
            [t,ok] = operation("Taxi","Refuel release","T9b",t,berthTime,taxiBerth,"Service",true,true);
        end
        if ok
            t = postFluidTopup(t,dock);
            a=t; t=t+cfg.atlanteTanker_min*60;
            append("Atlante","External tanker","A4",a,t,false,dock,"Service",0);
            t = postFluidTopup(t,dock);
            for pair=unique(cfg.assemblyPair(1:min(numel(indices),numel(cfg.assemblyPair))))
                a=t; t=t+cfg.kikiRefuel_min*60;
                append("Kiki "+pair,"Kiki refuel","K8",a,t,false,dock,"Service",0);
                append("Atlante","Supply Kiki "+pair+" propellant","A4",a,t,false,dock,"Service",0);
                t = postFluidTopup(t,dock);
            end
        end
        stop = max(stop,t);
    else
        block("Atlante","Refueling sequence",stop,"Refueling awaits completed retrieval and assembly operations.");
    end
end
stop = stop+cfg.standbyAfter_min*60;
if stop==0
    % Retain one EPS interval when the first departure is blocked at epoch.
    stop = min(cfg.sampleStep_s,cfg.maxDays*86400);
end
if ~taxiDocked && stop>taxiEnd
    append("Taxi","Post-stop state open","T8",taxiEnd,stop,false, ...
        @(q) nan(6,numel(q)),"Open",0);
end
knownStops = strings(0,1);
for j=1:numel(safeStops)
    held = safeStops{j};
    knownStops(end+1,1) = held.kiki; %#ok<AGROW>
    if stop>held.t
        append(held.kiki,"Attached wait at "+held.tag,"K11",held.t,stop,false,held.state,held.tag,0);
        if cfg.retainRangoForRelocations
            append(held.rango,"Attached wait at "+held.tag,"R8",held.t,stop,false,held.state,held.tag,0);
        end
    end
end
for v = unique(issues.Vehicle(startsWith(issues.Vehicle,"Kiki") | startsWith(issues.Vehicle,"Rango")))'
    if any(knownStops==v), continue; end
    a = min(issues.Epoch_s(issues.Vehicle==v));
    if stop>a
        mode = "K9";
        if startsWith(v,"Rango"), mode="R6"; end
        append(v,"Post-stop state open",mode,a,stop,false,@(q) nan(6,numel(q)),"Open",0);
    end
end
if cfg.dockedPairSupport, dockedPower(stop); end
assert(stop<=cfg.maxDays*86400,'Study exceeds the configured solar ephemeris horizon.');
plan = struct('events',sortrows(events,'Start_s'),'links',links,'holds',holds, ...
    'issues',issues,'paths',{paths},'stop_s',stop,'atlante',atlante,'taxiPark',park, ...
    'taxiEnd_s',taxiEnd,'taxiDocked',taxiDocked,'taxiBerth',taxiBerth,'completedModules',completed,'assemblySorties',sorties, ...
    'returnLambert',returnPlan,'recoveryDV_mps',recoveryDV, ...
    'retrievalSummary',retrievalSummary,'retrievalArcs',{retrievalArcs});

    function [finish,state,ok] = waitAtParking(time,dt)
        x0 = phase.taxiParkingDepartureState(:);
        chief = atlante(time);
        n = sqrt(mu/norm(chief(1:3))^3);
        hazards = obstacle([0;0;0],cfg.rpoHoldKeepOut_m,cfg.atlanteHoldGeometry);
        [drift,relative] = holdScreen(x0,dt,n,hazards);
        ok = drift.ok;
        finish = time;
        state = [];
        label = "Parking hold before outbound Lambert";
        holds(end+1,:) = {"Taxi",label,time,time+dt,drift.maxDrift,drift.clearance,ok};
        if ~ok
            block("Taxi",label,time,"Retimed departure requires a parking wait that fails drift or clearance.");
            return;
        end
        finish = time+dt;
        path = @(q) FABER_environment_eci(atlante(q),hcwStates(x0,q-time,n));
        state = FABER_environment_eci(atlante(finish),relative(:,end));
        append("Taxi",label,"T10",time,finish,false,path,"Atlante",0);
    end

    function arc = targetModule(time,chief,departureState)
        state = chief(time);
        period = 2*pi*sqrt(norm(state(1:3))^3/mu);
        rev = phase.outboundLambertRevolutions;
        if isempty(rev), rev = phase.phasingRevs; end
        arc = FABER_lambert_handoff(mu/1e9,Re/1000,departureState,@(dt) chief(time+dt), ...
            phase.x0_client,origin+time,period,rev,phase.outboundLambertTOF_s, ...
            phase.outboundMinimumPerigee_km);
        if ~isempty(arc)
            arc.moduleDepartureECI = state;
            arc.moduleArrivalECI = arc.referenceArrivalECI;
            arc.atlanteDepartureECI = atlante(time);
        end
    end

    function [finish,mass_kg,ok,departureState] = leaveBerth(time,mass_kg,startState)
        finish = time;
        departureState = [];
        goal = phase.taxiParkingDepartureState(:);
        chiefState = atlante(time);
        n = sqrt(mu/norm(chiefState(1:3))^3);
        dt = cfg.taxiDeparture_min*60;
        retreat = HCW_solver(startState(1:3),startState(4:6),goal(1:3),goal(4:6), ...
            dt,n,max(3,ceil(dt/min(5,cfg.sampleStep_s))+1));
        path = retreat.trajectory(1:3,:);
        range = vecnorm(path);
        inside = find(range<cfg.rpoHoldKeepOut_m);
        % Start at the serviced berth. Permit outward exit through the KOZ,
        % and reject inward motion or re-entry before the parking handoff.
        invalidExit = any(diff(inside)>1) || any(diff(range(inside))<-1e-6) || ...
            (~isempty(inside) && (inside(1)~=1 || inside(end)==numel(range)));
        label = "Taxi departure alignment";
        if invalidExit || FABER_path_clearance(path,cfg.atlanteHoldGeometry)<=0
            block("Taxi",label,time,"Berth-to-parking trajectory fails the departure clearance screen.");
            ok = false;
            return;
        end
        nextMass = burnMass(mass_kg,retreat.dvTotal,cfg.recoveryIsp_s);
        ok = massOK(nextMass,0,label,time);
        if ~ok, return; end
        if enforce
            [finish,ok] = movement("Taxi",label,"T5",time,dt,retreat.trajectory, ...
                startState,atlante,n,cfg.atlanteHoldGeometry,"Atlante",true,true,retreat.dvTotal,0);
        else
            relative = @(q) interp1(retreat.t,retreat.trajectory',min(dt,max(0,q(:)-time)),'linear')';
            finish = time+dt;
            append("Taxi",label,"T5",time,finish,true, ...
                @(q) FABER_environment_eci(atlante(q),relative(q)),"Atlante",retreat.dvTotal);
        end
        if ok
            mass_kg = nextMass;
            departureState = FABER_environment_eci(atlante(finish),retreat.finalState);
        end
    end

    function [time,endMass,berth,returnArc] = retrieve(arc,chief,time,mass_kg,retarget)
        endMass = NaN;
        berth = park;
        returnArc = [];
        outMass = burnMass(mass_kg,burnCost(arc),phase.farFieldIsp_s);
        if ~massOK(outMass,0,"Outbound Lambert",time), return; end

        arrival = arc.arrivalEpoch_s-origin;
        retrievalSummary.Departure_s(activeTrip) = time;
        append("Taxi","Outbound Lambert","T1",time,arrival,false,arcPath(arc),"Far field",burnCost(arc));
        time = arrival;
        if isempty(record)
            block("Taxi","Client RPO",time,"Phase A has no complete selected RPO pair.");
            return;
        end
        client = record.Client;
        if ~matchesGate(arc,client,"Client",time), return; end
        client = setGate(client,arc.arrivalLVLH);
        geometry = obstacle([0;0;0],cfg.rpoHoldKeepOut_m,cfg.clientHoldGeometry);
        before = recoveryDV;
        [time,ok] = approach("Taxi","Client","T2",client,chief,time,geometry);
        if ~ok, return; end
        captured = rpoMass(outMass,client,record.massOutboundStart_kg,before)+phase.payloadMass_kg;
        if ~massOK(captured,phase.payloadMass_kg,"Capture operations",time), return; end

        captureHold = @(q) fixedLVLH(chief,q,client.terminalLVLH);
        [time,ok] = operation("Taxi","Capture operations","T3",time,captureActive_s, ...
            captureHold,"Client",cfg.captureRequiresSun,false);
        if ~ok, return; end
        retrievalSummary.Capture_s(activeTrip) = time;
        retrievalSummary.CapturedMass_kg(activeTrip) = captured;

        remainder = phase.captureWait_min*60-captureActive_s;
        if remainder>1e-6
            [time,ok] = operation("Taxi","Attached capture dwell","T10",time,remainder, ...
                captureHold,"Client",false,true);
            if ~ok, return; end
        end

        if enforce || retarget
            client.elapsed_s = time-arrival-phase.captureWait_min*60;
            [returnArc,why] = FABER_return_lambert(mu/1e9,Re/1000,phase,arc,client,captured);
        else
            returnArc = budget.returnLambert;
            why = "Phase A return transfer is open.";
        end
        if isempty(returnArc)
            block("Taxi","Return Lambert",time,why);
            return;
        end
        returnMass = burnMass(captured,burnCost(returnArc),phase.farFieldIsp_s);
        if ~massOK(returnMass,phase.payloadMass_kg,"Return Lambert",time), return; end

        start = returnArc.departureEpoch_s-origin;
        finish = returnArc.arrivalEpoch_s-origin;
        append("Taxi","Return Lambert","T4",start,finish,false,arcPath(returnArc),"Far field",burnCost(returnArc));
        time = finish;
        host = record.ATLANTE;
        if ~matchesGate(returnArc,host,"Atlante",time), return; end
        host = setGate(host,returnArc.arrivalLVLH);
        geometry = obstacle([0;0;0],cfg.rpoHoldKeepOut_m,cfg.atlanteHoldGeometry);
        before = recoveryDV;
        [time,ok] = approach("Taxi","Atlante","T5",host,atlante,time,geometry);
        if ~ok, return; end
        unloaded = rpoMass(returnMass,host,record.massReturnStart_kg,before)-phase.payloadMass_kg;
        if ~massOK(unloaded,0,"Taxi to Atlante handover",time), return; end

        berth = @(q) fixedLVLH(atlante,q,host.terminalLVLH);
        [time,ok] = operation("Taxi","Taxi to Atlante handover","T6",time, ...
            cfg.taxiHandover_min*60,berth,"Atlante",true,false);
        if ok, endMass = unloaded; end
    end

    function next = rpoMass(current,approachRecord,screenedMass,previousRecovery)
        ratio = 1-approachRecord.propellant_kg/screenedMass;
        next = burnMass(current*ratio,recoveryDV-previousRecovery,cfg.recoveryIsp_s);
    end

    function ok = massOK(current,payload,activity,time)
        ok = isfinite(current) && current>payload && ...
            (isnan(cfg.taxiDryMass_kg) || current>=cfg.taxiDryMass_kg+payload);
        if ~ok
            block("Taxi",activity,time,"Taxi mass falls below the configured usable propellant limit.");
        end
    end

    function ok = matchesGate(arc,approachRecord,leg,time)
        residual = arc.arrivalLVLH-approachRecord.transferDetails{1}.startState;
        ok = norm(residual(1:3))<=1 && norm(residual(4:6))<=0.001;
        if ~ok
            block("Taxi",leg+" RPO",time,"Retargeted gate does not match the screened RPO entry.");
        end
    end

    function r = setGate(r,state)
        r.trace{1}(:,1) = state;
        r.transferDetails{1}.startState = state;
    end

    function [finish,ok] = approach(veh,leg,mode,r,chief,time,hazards)
        finish = time; ok = true;
        chiefAtStart = chief(time);
        n = sqrt(mu/norm(chiefAtStart(1:3))^3);
        for k=1:numel(r.trace)
            activity = leg+" "+string(budget.RPO_MC.phaseNames{k});
            dt = r.tofPhase(k)*60;
            if dt<=0, continue; end
            trace = r.trace{k};
            if enforce && k>=7
                % Protect the uninterrupted closing chain through the HP4
                % fly-through and include the immediately following service.
                begin = r.transferDetails{4}.startState;
                bounds = rpoBounds(hazards,begin);
                [trace,dt] = terminalChain(r,cfg.sampleStep_s);
                tail = cfg.taxiHandover_min*60;
                if leg=="Client"
                    tail = captureActive_s*double(cfg.captureRequiresSun);
                end
                if ~isfinite(tail), tail=0; end
                [finish,ok] = movement(veh,leg+" closing and final approach",mode, ...
                    finish,dt,trace,begin,chief,n,bounds,leg,false,true,sum(r.dvPhase(7:9)),tail);
                return;
            end
            if mod(k,2)==0
                waitState = trace(:,1);
                bounds = rpoBounds(hazards,waitState);
                [finish,ok] = hold(veh,activity,mode,finish,dt,waitState,chief,n,bounds,leg);
            else
                begin = r.transferDetails{(k+1)/2}.startState;
                bounds = rpoBounds(hazards,begin);
                [finish,ok] = movement(veh,activity,mode,finish,dt,trace, ...
                    begin,chief,n,bounds,leg,false,true,r.dvPhase(k),0);
            end
            if ~ok, return; end
        end
    end

    function [finish,ok,departure] = reservedRelocation(veh,tag,mode,r,time)
        % Reserve all three legs and any intermediate holds before leaving
        % the attached location. Discard the ledger of a failed preview.
        counts = [height(events),numel(paths),height(holds),height(issues)];
        beforeDV = recoveryDV;
        departure = time;
        offsets = 0;
        if enforce
            offsets = 0:cfg.assemblyDepartureSearchStep_s:cfg.maxWait_s;
        end
        failure = "";
        for delay = offsets
            if time+delay+(r.exitTOF_min+r.transferTOF_min+r.finalTOF_min)*60>cfg.maxDays*86400
                break;
            end
            if delay>0
                x0 = [r.source_m(:);zeros(3,1)];
                append(veh,"Attached wait before "+tag,"K11",time,time+delay,false, ...
                    @(q) fixedLVLH(atlante,q,x0),tag,0);
            end
            [finish,ok] = relocation(veh,tag,mode,r,time+delay);
            if ok
                first = find(events.Activity==tag+" departure" & events.Vehicle==veh,1,'last');
                departure = events.Start_s(first);
                if ~enforce || departure-time<=cfg.maxWait_s, return; end
                failure = "The departure lies beyond the configured wait limit.";
            end
            if failure=="" && height(issues)>counts(4)
                failure = issues.Activity(end)+" - "+issues.Reason(end);
            end
            events = events(1:counts(1),:);
            paths = paths(1:counts(2),1);
            holds = holds(1:counts(3),:);
            issues = issues(1:counts(4),:);
            recoveryDV = beforeDV;
        end
        finish = time;
        ok = false;
        block(veh,tag,time,"Stayed attached. No safe relocation schedule within the search horizon. "+failure);
    end

    function [finish,ok] = relocation(veh,tag,transitMode,r,time)
        finish = time; ok = true;
        n = sqrt(mu/(Re+1000*phase.h_target_km)^3);
        hazards = obstacle([r.source_m(:),r.target_m(:)],cfg.assemblyHoldKeepOut_m,cfg.assemblyHoldGeometry);
        pathsHere = {r.exitPath,r.transferPath,r.finalPath};
        minutes = [r.exitTOF_min,r.transferTOF_min,r.finalTOF_min];
        cost = [r.exitDV_mps,r.transferDV_mps,r.finalDV_mps];
        labels = [" departure"," gate transfer"," terminal arrival"];
        for k=1:3
            trace = pathsHere{k};
            begin = [trace(1:3,1);zeros(3,1)];
            [finish,ok] = movement(veh,tag+labels(k),transitMode,finish,minutes(k)*60, ...
                trace,begin,atlante,n,hazards,tag,k==1,true,cost(k),0);
            if ~ok, return; end
        end
    end

    function [finish,ok] = movement(veh,label,mode,time,dt,trace,x0,chief,n,hazards,leg,attached,canWait,dv,tail)
        if size(trace,1)==3, trace=[trace;zeros(size(trace))]; end
        step = linspace(0,dt,max(2,ceil(dt/cfg.sampleStep_s)+1));
        relative = interp1(linspace(0,dt,size(trace,2)),trace',step,'linear')';
        pathAt = @(a,t) FABER_environment_eci(chief(t), ...
            interp1(step,relative',min(dt,max(0,t(:)-a)),'linear')');
        coverage = step;
        if tail>0
            coverage = unique([step,dt+linspace(0,tail,max(2,ceil(tail/cfg.sampleStep_s)+1))]);
        end
        start = daylight(time,dt+tail,pathAt,coverage,chief);
        finish = time; ok = isfinite(start);
        if ~ok
            required_min = (dt+tail+2*cfg.sunGuard_s)/60;
            reason = string(sprintf("Needs %.1f min of continuous sunlight, including %.1f min for " + ...
                "the following operation and %.1f min of guard time. No window fits within %.1f h.", ...
                required_min,tail/60,2*cfg.sunGuard_s/60,cfg.maxWait_s/3600));
            block(veh,label,time,reason);
            return;
        end
        if start>time
            if ~canWait
                block(veh,label,time,"Sunlight unavailable at a fly-through gate. No hold is authorized.");
                ok = false; return;
            end
            if attached
                waitMode = holdMode(veh);
                append(veh,"Attached wait before "+label,waitMode,time,start,false, ...
                    @(q) fixedLVLH(chief,q,x0),leg,0);
            else
                [start,ok] = waitAndRecover(veh,label,mode,time,start,dt,relative,step,coverage, ...
                    x0,chief,n,hazards,leg);
                if ~ok, return; end
            end
        end
        finish = start+dt;
        x = @(t) FABER_environment_eci(chief(t), ...
            interp1(step,relative',min(dt,max(0,t(:)-start)),'linear')');
        append(veh,label,mode,start,finish,true,x,leg,dv);
    end

    function [start,ok] = waitAndRecover(veh,label,mode,time,start,dt,relative,step,coverage,x0,chief,n,hazards,leg)
        [drift,rel] = holdScreen(x0,start-time,n,hazards);
        holds(end+1,:) = {veh,label,time,start,drift.maxDrift,drift.clearance,drift.ok};
        ok = drift.ok;
        if ~ok
            block(veh,label,time,"Eclipse wait fails the passive drift or clearance screen.");
            return;
        end
        endState = rel(:,end);
        append(veh,"Eclipse wait before "+label,holdMode(veh),time,start,false, ...
            @(t) FABER_environment_eci(chief(t),hcwStates(x0,t-time,n)),leg,0);
        if norm(endState-x0)<1e-5
            return;
        end
        trim_s = cfg.recoveryTime_s;
        trim = HCW_solver(endState(1:3),endState(4:6),x0(1:3),x0(4:6),trim_s,n, ...
            max(3,ceil(trim_s/cfg.sampleStep_s)+1));
        check = FABER_path_clearance(trim.trajectory(1:3,:),hazards)>0;
        if isfield(hazards,'maxRange_m')
            check = check && max(vecnorm(trim.trajectory(1:3,:)))<hazards.maxRange_m;
        end
        trimTimes = start+trim.t;
        trimChief = chief(trimTimes);
        trimLit = lit(FABER_environment_eci(trimChief,trim.trajectory),trimTimes) && ...
            lit(trimChief,trimTimes);
        guard = coverage(end)+linspace(0,cfg.sunGuard_s, ...
            max(2,ceil(cfg.sunGuard_s/cfg.sampleStep_s)+1));
        moveOffsets = unique([coverage,guard]);
        moveTimes = start+trim_s+moveOffsets;
        moveRelative = interp1(step,relative',min(dt,moveOffsets),'linear')';
        moveChief = chief(moveTimes);
        moveLit = lit(FABER_environment_eci(moveChief,moveRelative),moveTimes) && ...
            lit(moveChief,moveTimes);
        if moveTimes(end)>cfg.maxDays*86400 || ~check || ~trimLit || ~moveLit || ...
                trim.dvTotal>cfg.maxRecoveryDV_mps || ...
                max(trim.relativeSpeed)>cfg.maxRecoverySpeed_mps
            block(veh,label,start,"Safe drift, but sunlit reacquisition and its next transfer did not pass.");
            ok = false; return;
        end
        trimStart = start;
        append(veh,"Reacquire hold before "+label,mode,start,start+trim_s,true, ...
            @(q) FABER_environment_eci(chief(q),interp1(trim.t,trim.trajectory', ...
            min(trim_s,max(0,q(:)-trimStart)),'linear')'),leg,trim.dvTotal);
        recoveryDV = recoveryDV+trim.dvTotal;
        start = start+trim_s;
    end

    function [finish,ok] = hold(veh,label,mode,time,dt,x0,chief,n,hazards,leg)
        finish = time+dt;
        q = linspace(time,finish,max(2,ceil(dt/cfg.sampleStep_s)+1));
        shadow = ~lit(fixedLVLH(chief,q,x0),q);
        ok = true;
        if enforce && shadow
            [drift,rel] = holdScreen(x0,dt,n,hazards);
            restore = norm(rel(1:3,end)-x0(1:3))<=cfg.holdContinuationPosition_m && ...
                norm(rel(4:6,end)-x0(4:6))<=cfg.holdContinuationVelocity_mps;
            holds(end+1,:) = {veh,label,time,finish,drift.maxDrift,drift.clearance,drift.ok && restore};
            if ~(drift.ok && restore)
                block(veh,label,time,"Phase A hold cannot drift through eclipse and resume within tolerances.");
                finish = time; ok = false; return;
            end
            path = @(q) FABER_environment_eci(chief(q),hcwStates(x0,q-time,n));
        else
            path = @(q) fixedLVLH(chief,q,x0);
        end
        append(veh,label,holdMode(veh),time,finish,false,path,leg,0);
    end

    function [result,relative] = holdScreen(x0,dt,n,hazards)
        q = linspace(0,dt,max(2,ceil(dt/cfg.sampleStep_s)+1));
        [relative,phi] = hcwStates(x0,q,n);
        maxDrift = 0; clearance = Inf;
        for k=1:numel(q)
            % Circumscribed error radius of the propagated initial-state box.
            radius = norm(abs(phi(1:3,:,k))*cfg.holdErrorLVLH(:));
            if isfield(hazards,'extraRadius_m')
                radius = radius+hazards.extraRadius_m;
            end
            maxDrift = max(maxDrift,norm(relative(1:3,k)-x0(1:3))+radius);
            clearance = min(clearance,FABER_clearance_m(relative(1:3,k),hazards,radius));
            if isfield(hazards,'maxRange_m')
                clearance = min(clearance,hazards.maxRange_m-norm(relative(1:3,k))-radius);
            end
        end
        result = struct('maxDrift',maxDrift,'clearance',clearance, ...
            'ok',clearance>0 && maxDrift<=cfg.maxHoldDrift_m);
    end

    function bounds = rpoBounds(bounds,x0)
        az = budget.RPO_MC.cfg.AZ_m;
        if norm(x0(1:3))>=az
            bounds.radii_m(1) = max(bounds.radii_m(1),az);
        else
            bounds.maxRange_m = az;
        end
    end

    function [finish,ok] = operation(veh,label,mode,time,dt,x,leg,sunRequired,attached)
        finish = time;
        ok = isfinite(dt) && dt>=0;
        if ~ok
            append(veh,label+" (duration open)",mode,time,time,sunRequired,x,leg,NaN);
            block(veh,label,time,"Duration is open. Set the handover or service duration explicitly.");
            return;
        end
        if dt==0
            append(veh,label,mode,time,time,sunRequired,x,leg,0);
            return;
        end
        start = time;
        if sunRequired
            offset = linspace(0,dt,max(2,ceil(dt/cfg.sampleStep_s)+1));
            chief = atlante;
            if leg=="Client", chief=clientChief; end
            start = daylight(time,dt,@(~,q) x(q),offset,chief);
        end
        if ~isfinite(start) || (start>time && ~attached)
            block(veh,label,time,"No sunlight window with a confirmed attached hold.");
            ok = false; return;
        end
        if start>time
            append(veh,"Attached wait before "+label,holdMode(veh),time,start,false,x,leg,0);
        end
        finish = start+dt;
        append(veh,label,mode,start,finish,sunRequired,x,leg,0);
    end

    function start = daylight(time,dt,x,offset,target)
        start = time;
        if ~enforce, return; end
        for t0 = time+(0:cfg.searchStep_s:cfg.maxWait_s)
            if t0+dt+cfg.sunGuard_s>cfg.maxDays*86400, break; end
            q = t0+offset;
            % Hold each LVLH endpoint fixed during the adjacent lighting guard.
            guardSamples = max(2,ceil(cfg.sunGuard_s/cfg.sampleStep_s)+1);
            guardedTime = unique([linspace(max(0,t0-cfg.sunGuard_s),t0,guardSamples), ...
                q,linspace(q(end),q(end)+cfg.sunGuard_s,guardSamples)]);
            states = x(t0,guardedTime);
            if lit(states,guardedTime) && lit(target(guardedTime),guardedTime)
                start = t0; return;
            end
        end
        start = NaN;
    end

    function ok = lit(states,t)
        light = FABER_environment_light(solar,t,states,cfg.aspectAxisLVLH,Re);
        ok = all(light>=cfg.sunRequiredFraction);
    end

    function fn = arcPath(arc)
        a = arc.departureEpoch_s-origin;
        dt = arc.transferTime_s;
        opts = odeset('RelTol',1e-11,'AbsTol',1e-8);
        sol = ode113(@(~,x) [x(4:6);-mu*x(1:3)/norm(x(1:3))^3], ...
            [0 dt],arc.departurePostBurnECI,opts);
        fn = @(t) deval(sol,min(dt,max(0,t(:)'-a)));
    end

    function append(veh,label,mode,a,b,sunRequired,x,leg,dv)
        assert(b>=a && b<=cfg.maxDays*86400,'Activity outside the study horizon.');
        events(end+1,:) = {veh,label,mode,a,b,sunRequired,dv,activeTrip};
        paths{end+1,1} = struct('vehicle',veh,'label',label,'leg',leg,'start',a,'stop',b, ...
            'state',x,'sunRequired',sunRequired,'clientChief',clientChief);
    end

    function block(veh,label,t,why)
        issues(end+1,:) = {veh,label,t,string(why)};
    end

    function addLink(source,destination,a,b,power,kind)
        if nargin<6, kind="Charge"; end
        if b>a && power>0
            links(end+1,:) = {source,destination,a,b,power,kind};
        end
    end

    function carryRango(first,rango,label)
        last = numel(paths);
        for i=first:last
            p = paths{i};
            append(rango,label,"R8",p.start,p.stop,p.sunRequired,p.state,p.leg,0);
        end
    end

    function [finish,ok] = prepareCharge(kiki,rango,tag,time,home,minutes)
        finish = time;
        ok = true;
        if minutes(2)>0
            [finish,ok] = operation(kiki,tag+" connect Rango charge link","K6c",finish, ...
                cfg.rangoChargeConnect_min*60,home,tag,false,true);
            if ~ok, return; end
            a = finish;
            [finish,ok] = operation(kiki,tag+" charge Rango at Atlante","K7",a, ...
                minutes(2)*60,home,tag,false,true);
            if ~ok, return; end
            append(rango,tag+" receive charge","R5",a,finish,false,home,tag,0);
            addLink(kiki,rango,a,finish,cfg.rangoPort_kW);
            addLink("Atlante",kiki,a,finish,cfg.kikiPort_kW);
        end
        if minutes(1)>0
            a = finish;
            [finish,ok] = operation(kiki,tag+" recharge at Atlante","K7",a, ...
                minutes(1)*60,home,tag,false,true);
            if ~ok, return; end
            addLink("Atlante",kiki,a,finish,cfg.kikiPort_kW);
        end
    end

    function [finish,ok] = rangoTransfer(kiki,rango,label,start,site,tag)
        [finish,ok] = operation(kiki,label,"K5",start,cfg.rangoTransfer_min*60,site,tag,true,true);
        if ok && cfg.rangoTransfer_min>0
            append(rango,label,"R1",finish-cfg.rangoTransfer_min*60,finish,true,site,tag,0);
        end
    end

    function finish = postFluidTopup(start,dock)
        finish = start+cfg.postFluidTopup_min*60;
        if finish==start, return; end
        for pair=1:2
            kiki = "Kiki "+pair; rango = "Rango "+pair;
            append(kiki,"Post-refuel top-up","K7",start,finish,false,dock,"Service",0);
            addLink("Atlante",kiki,start,finish,cfg.kikiPort_kW);
            if cfg.retainRangoForRelocations
                append(rango,"Retained top-up","R5",start,finish,false,dock,"Service",0);
                addLink(kiki,rango,start,finish,cfg.rangoPort_kW);
            end
        end
    end

    function dockedPower(stop)
        % Cover standby loads through the PM4 retained electrical topology.
        % Hold the shared fluid interface clear during Atlante A4 operations.
        input = fullfile(fileparts(mfilename('fullpath')),'input', ...
            "power_modes_"+cfg.cluster+".csv");
        power = readtable(input,'TextType','string');
        e = events;
        scheduledLinks = links;
        edges = unique([0;stop;e.Start_s;e.End_s;scheduledLinks.Start_s;scheduledLinks.End_s]);
        for j=1:numel(edges)-1
            a = edges(j); b = edges(j+1); m = (a+b)/2;
            active = e.Start_s<=m & e.End_s>m;
            if any(active & e.Vehicle=="Atlante" & e.Mode=="A4")
                continue;
            end
            for pair=1:2
                offDock = away{pair};
                if any(offDock(:,1)<=m & offDock(:,2)>m)
                    continue;
                end
                kiki = "Kiki "+pair; rango = "Rango "+pair;
                if any(active & (e.Vehicle==kiki | e.Vehicle==rango) & e.Activity=="Post-stop state open")
                    continue;
                end
                receiving = scheduledLinks.Start_s<=m & scheduledLinks.End_s>m;
                rangoPower = 0;
                if cfg.retainRangoForRelocations && ...
                        ~any(receiving & scheduledLinks.Receiver==rango)
                    q = e.Mode(active & e.Vehicle==rango);
                    if isempty(q), q="R8"; end
                    rangoPower = power.Mean_W(power.Mode==q(1) & power.Vehicle=="Rango")/1000;
                    addLink(kiki,rango,a,b,rangoPower,"Docked support");
                end
                if ~any(receiving & scheduledLinks.Receiver==kiki)
                    q = e.Mode(active & e.Vehicle==kiki);
                    if isempty(q), q="KB"; end
                    own = power.Mean_W(power.Mode==q(1) & power.Vehicle=="Kiki")/1000;
                    addLink("Atlante",kiki,a,b,own+rangoPower/0.90,"Docked support");
                end
            end
        end
    end
end

function [path,duration] = terminalChain(r,sampleStep)
times = [];
states = [];
duration = 0;
for k=7:9
    dt = 60*r.tofPhase(k);
    if dt<=0, continue; end
    q = duration+linspace(0,dt,size(r.trace{k},2));
    use = true(size(q));
    if ~isempty(times), use(1)=false; end
    times = [times,q(use)]; %#ok<AGROW>
    states = [states,r.trace{k}(:,use)]; %#ok<AGROW>
    duration = duration+dt;
end
q = linspace(0,duration,max(2,ceil(duration/sampleStep)+1));
path = interp1(times,states',q,'linear')';
end

function x = fixedLVLH(chief,t,state)
x = FABER_environment_eci(chief(t),repmat(state(:),1,numel(t)));
end

function geometry = obstacle(points,radius,geometry)
centers = points';
radii = repmat(radius,size(points,2),1);
if isfield(geometry,'centers_m')
    centers = [centers;geometry.centers_m];
    if isfield(geometry,'radii_m')
        radii = [radii;geometry.radii_m(:)];
    else
        radii = [radii;zeros(size(geometry.centers_m,1),1)];
    end
end
geometry.centers_m = centers;
geometry.radii_m = radii;
end

function mode = holdMode(veh)
mode = "K11";
if veh=="Taxi", mode = "T10"; end
if startsWith(veh,"Rango"), mode = "R7"; end
end

function [x,phi] = hcwStates(x0,t,n)
A = [zeros(3) eye(3); 3*n^2 0 0 0 2*n 0;0 0 0 -2*n 0 0;0 0 -n^2 0 0 0];
x = zeros(6,numel(t)); phi = zeros(6,6,numel(t));
for i=1:numel(t)
    phi(:,:,i) = expm(A*t(i));
    x(:,i) = phi(:,:,i)*x0;
end
end

function v = burnCost(arc)
v = norm(arc.departureBurn_mps)+norm(arc.arrivalBurn_mps);
end

function mass = burnMass(initial,dv,isp)
mass = initial*exp(-dv/(9.80665*isp));
end
