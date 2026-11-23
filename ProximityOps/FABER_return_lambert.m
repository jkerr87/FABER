function [transfer,reason] = FABER_return_lambert(mu,Re,cfg,outbound,client,massCaptured_kg)
% Model ideal instantaneous attachment at the outbound RPO terminal state.
transfer = [];
reason = "Capture dwell is unspecified";
if isnan(cfg.captureWait_min)
    return;
end
reason = "No Lambert arc passed the return perigee and endpoint checks";
n = sqrt(mu/(Re+cfg.h_target_km)^3);
captureEpoch_s = outbound.arrivalEpoch_s+client.elapsed_s;
module = circularCoast(mu,outbound.moduleArrivalECI,client.elapsed_s);
taxi = FABER_relative_state(module,client.terminalLVLH,'lvlhToEci');
% Use the position and momentum of the combined center of mass.
payloadFraction = cfg.payloadMass_kg/massCaptured_kg;
stackAtCapture = taxi + payloadFraction*(module-taxi);
dwell_s = 60*cfg.captureWait_min;
stackAtDeparture = coast(mu,stackAtCapture,dwell_s);
departureEpoch_s = captureEpoch_s+dwell_s;
atlanteAtDeparture = circularCoast(mu,outbound.atlanteDepartureECI, ...
    departureEpoch_s-outbound.departureEpoch_s);
target = @(dt) circularCoast(mu,atlanteAtDeparture,dt);
revolutions = cfg.returnLambertRevolutions;
if isempty(revolutions)
    revolutions = cfg.phasingRevs;
end
transfer = FABER_lambert_handoff(mu,Re,stackAtDeparture,target,cfg.x0_atlante, ...
    departureEpoch_s,2*pi/n,revolutions,cfg.returnLambertTOF_s, ...
    cfg.outboundMinimumPerigee_km);
if isempty(transfer)
    return;
end
transfer.captureEpoch_s = captureEpoch_s;
transfer.captureDwell_s = dwell_s;
transfer.captureTaxiLVLH = client.terminalLVLH;
transfer.captureTaxiECI = taxi;
transfer.captureModuleECI = module;
transfer.captureStackECI = stackAtCapture;
transfer.atlanteDepartureECI = atlanteAtDeparture;
transfer.massCaptured_kg = massCaptured_kg;
dv = norm(transfer.departureBurn_mps)+norm(transfer.arrivalBurn_mps);
transfer.massArrival_kg = massCaptured_kg*exp(-dv/(9.80665*cfg.farFieldIsp_s));
reason = "";
end

function x = circularCoast(mu,x0,dt)
r = x0(1:3);
v = x0(4:6);
n = sqrt(1e9*mu/norm(r)^3);
a = n*dt;
x = [cos(a)*r+sin(a)/n*v; cos(a)*v-n*sin(a)*r];
end

function x = coast(mu,x0,dt)
x = x0;
if dt>0
    opts = odeset('RelTol',1e-11,'AbsTol',1e-12);
    [~,states] = ode113(@(~,y) [y(4:6);-1e9*mu*y(1:3)/norm(y(1:3))^3], ...
        [0 dt],x0,opts);
    x = states(end,:)';
end
end
