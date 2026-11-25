function firing = FABER_transfer_firing(out)
% Report non-RPO firing times; keep coast durations and burn times separate.
cfg = out.cfg;
T = out.phaseTable;
segments = string(T.Segment);
mc = out.RPO_MC.cfg;
main = mc.propulsion.main;
rcs = mc.propulsion.rcs;
main.Isp_s = cfg.farFieldIsp_s;
rows = cell(0,13);

for ship = ["Taxi","ATLANTE"]
    if ship=="Taxi"
        mass = cfg.massTaxi_kg;
        propulsion = main;
        name = "MAIN";
    else
        mass = cfg.massAtlante_kg+cfg.massKiki_kg+cfg.massRango_kg;
        propulsion = cfg.atlanteTransferPropulsion;
        name = "Configured";
    end
    for segment = ["Deployment","Inter-client"]
        selected = T.Vehicle==ship & segments==segment & string(T.Phase)~="RAAN reserve";
        cost = sum(T.DV_perEvent_mps(selected));
        label = "Injection correction";
        if segment=="Inter-client", label = "Client orbit transfer"; end
        add(ship,label,"Combined nominal burns",cost,mass,propulsion,name,NaN,false);
    end

    for phase = ["Parking acquisition","Departure state alignment correction"]
        selected = T.Vehicle==ship & T.Architecture==ship+" Vbar" & string(T.Phase)==phase;
        if ~any(selected), continue; end
        parkingPropulsion = propulsion;
        parkingType = name;
        if ship=="Taxi"
            parkingPropulsion = rcs;
            parkingType = "RCS";
        else
            parkingPropulsion = cfg.atlanteParkingPropulsion;
        end
        add(ship,phase,"V-bar parking",T.DV_perEvent_mps(selected),mass, ...
            parkingPropulsion,parkingType,NaN,false);
    end

    selected = T.Vehicle==ship & segments=="End of life" & string(T.Phase)=="Deorbit";
    disposalMass = cfg.disposalAtlanteMass_kg;
    if ship=="Taxi", disposalMass = cfg.disposalTaxiMass_kg; end
    if any(selected)
        add(ship,"Disposal","Insertion burn only",T.DV_perEvent_mps(selected), ...
            disposalMass,propulsion,name,NaN,false);
    end
end

d = out.taxiInitialDeparture;
dv = d.DV_mps;
if ~d.feasible, dv = NaN; end
add("Taxi","Atlante Departure","Continuous release",dv,cfg.massTaxi_kg, ...
    rcs,"RCS",d.burnTime_s,true);
lambert("Far field to Client",out.outboundLambert,cfg.massTaxi_kg);
returnMass = NaN;
if ~isempty(out.returnLambert), returnMass = out.returnLambert.massCaptured_kg; end
lambert("Far field to Atlante",out.returnLambert,returnMass);

for phase = ["Loaded relocation","Unloaded return"]
    selected = T.Vehicle=="Kiki" & segments=="Assembly" & string(T.Phase)==phase;
    if ~any(selected), continue; end
    mass = cfg.massKiki_kg+cfg.massRango_kg;
    if phase=="Loaded relocation", mass = mass+cfg.payloadMass_kg; end
    add("Kiki",phase,"Reported mean relocation with Rango attached", ...
        T.DV_perEvent_mps(selected),mass,cfg.kikiTransferPropulsion,"Configured",NaN,false);
end

firing = cell2table(rows,'VariableNames',{'Vehicle','Transfer','BurnOrScope','RawDV_mps', ...
    'BurnStartMass_kg','ThrusterType','InstalledThrusters','ConfiguredActiveThrusters','RatedThrustEach_N', ...
    'Isp_s','ThrustEfficiency','FiringTime_s','TimeBasis'});

    function lambert(label,arc,mass)
        dv = [NaN NaN];
        if ~isempty(arc), dv = [norm(arc.departureBurn_mps),norm(arc.arrivalBurn_mps)]; end
        add("Taxi",label,"Departure burn",dv(1),mass,main,"MAIN",NaN,false);
        mass = mass*exp(-dv(1)/(9.80665*main.Isp_s));
        add("Taxi",label,"RPO gate velocity match",dv(2),mass,main,"MAIN",NaN,false);
    end

    function add(ship,label,scope,dv,mass,propulsion,type,measured,isContinuous)
        installed = NaN; count = NaN; thrust = NaN; isp = NaN; efficiency = NaN;
        time = NaN;
        basis = "Transfer open";
        if isempty(fieldnames(propulsion))
            type = "Not configured";
        else
            if isfield(propulsion,'type'), type = string(propulsion.type); end
            if isfield(propulsion,'nInstalled'), installed = propulsion.nInstalled; end
            count = propulsion.nActiveThrusters;
            thrust = propulsion.thrust_N;
            isp = propulsion.Isp_s;
            efficiency = propulsion.thrustEfficiency;
            validateattributes(count,{'numeric'},{'scalar','integer','positive'});
            validateattributes(thrust,{'numeric'},{'scalar','real','finite','positive'});
            validateattributes(isp,{'numeric'},{'scalar','real','finite','positive'});
            validateattributes(efficiency,{'numeric'},{'scalar','real','positive','<=',1});
        end
        if isfinite(dv)
            if dv==0
                time = 0;
                basis = "No burn";
            elseif isContinuous
                time = measured;
                basis = "Continuous-thrust simulation";
            elseif isfinite(count)
                exhaust = 9.80665*isp;
                time = -mass*expm1(-dv/exhaust)*exhaust/(count*thrust*efficiency);
                basis = "Constant-thrust estimate from raw impulsive delta-v";
            else
                basis = "Propulsion not configured";
            end
        end
        rows(end+1,:) = {ship,label,scope,dv,mass,type,installed,count,thrust,isp,efficiency,time,basis};
    end
end
