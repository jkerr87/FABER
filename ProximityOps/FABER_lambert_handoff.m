function transfer = FABER_lambert_handoff(mu,Re,departureECI,chiefAtTime, ...
    gateLVLH,epoch_s,period,rev,times,minPerigee_km)
% Search two-body Lambert arcs, then independently verify the chosen endpoint.
% Inputs and outputs use m, m/s and s except mu, Re and perigee (km units).
taxi = departureECI(:)/1000;
% Match the departure orbit's direction, including retrograde SSO cases.
angularMomentum = cross(taxi(1:3),taxi(4:6));
orbitType = double(angularMomentum(3)<0);
validateattributes(rev,{'numeric'},{'scalar','integer','nonnegative','finite'});
if isempty(times)
    times = period*(rev+(0.15:0.05:0.85));
end

validateattributes(times,{'numeric'},{'vector','positive','finite','real'});
validateattributes(epoch_s,{'numeric'},{'scalar','finite','real'});
validateattributes(minPerigee_km,{'numeric'},{'scalar','finite','real','nonnegative'});

% Rank candidates by both burns after screening their perigee altitude.
bestDV = Inf;
candidates = nan(numel(times)*(1+double(rev>0)),6);
candidate = 0;
for dt = times(:)'
    chief = chiefAtTime(dt);
    gate = FABER_relative_state(chief,gateLVLH,'lvlhToEci')/1000;
    for branch = 0:double(rev>0)
        [pAxis,p,e,err,vi,vf] = lambertMR(taxi(1:3)',gate(1:3)',dt,mu,orbitType,rev,branch,0);
        if err~=0 || ~all(isfinite([pAxis p e vi(:)' vf(:)']))
            continue
        end

        perigee = p/(1+e)-Re;
        departureBurn = vi(:)-taxi(4:6);
        arrivalBurn = gate(4:6)-vf(:);
        dv = norm(departureBurn)+norm(arrivalBurn);
        candidate = candidate+1;
        candidates(candidate,:) = [dt,branch,perigee, ...
            1000*norm(departureBurn),1000*norm(arrivalBurn),1000*dv];
        if perigee<minPerigee_km
            continue
        end

        if dv<bestDV
            bestDV = dv;
            transfer = struct('transferTime_s',dt,'departureEpoch_s',epoch_s, ...
                'arrivalEpoch_s',epoch_s+dt,'revolutions',rev, ...
                'branch',branch,'perigeeAltitude_km',perigee, ...
                'departureBurn_mps',1000*departureBurn,'arrivalBurn_mps',1000*arrivalBurn, ...
                'departurePreBurnECI',1000*taxi,'departurePostBurnECI',1000*[taxi(1:3);vi(:)], ...
                'referenceArrivalECI',chief,'gateTargetECI',1000*gate);
        end
    end
end

if ~isfinite(bestDV)
    transfer = [];
    return;
end

% Integrate the selected arc independently before handing it to RPO.
state0 = transfer.departurePostBurnECI/1000;
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[~,states] = ode113(@(~,x) [x(4:6);-mu*x(1:3)/norm(x(1:3))^3], ...
    [0 transfer.transferTime_s],state0,options);
arrival = 1000*states(end,:)';
if norm(arrival(1:3)-transfer.gateTargetECI(1:3)) >= 1 || ...
        norm(arrival(4:6)+transfer.arrivalBurn_mps-transfer.gateTargetECI(4:6)) >= 0.001
    transfer = [];
    return;
end

% Retain the candidate grid for the transfer-time trade plot.
transfer.search = array2table(candidates(1:candidate,:), ...
    'VariableNames',{'Duration_s','Branch','PerigeeAltitude_km', ...
    'DepartureDV_mps','ArrivalDV_mps','TotalDV_mps'});
transfer.search.PassesPerigee = transfer.search.PerigeeAltitude_km>=minPerigee_km;
transfer.arrivalPreBurnECI = arrival;
arrival(4:6) = arrival(4:6)+transfer.arrivalBurn_mps;
transfer.arrivalPostBurnECI = arrival;
transfer.arrivalLVLH = FABER_relative_state(transfer.referenceArrivalECI,arrival,'eciToLvlh');
end
