function transfer = FABER_outbound_lambert(mu,Re,cfg)
% Connect a nominal parking handoff to the module's RPO entry gate.
% Inertial calculation in km and s. Public state vectors use m and m/s.
rAtlante = Re+cfg.h_target_km;
rModule = rAtlante+cfg.moduleAltitudeOffset_km;
nModule = sqrt(mu/rModule^3);
period = 2*pi/nModule;
u0 = cfg.outboundDepartureArgumentOfLatitude_deg;
raan = cfg.outboundAtlanteRAAN_deg;
atlante = circularState(mu,rAtlante,cfg.i_target_deg,raan,u0);
taxi = FABER_relative_state(1000*atlante,cfg.taxiParkingDepartureState,'lvlhToEci');
rev = cfg.outboundLambertRevolutions;
if isempty(rev)
    rev = cfg.phasingRevs;
end

moduleAtTime = @(dt) 1000*circularState(mu,rModule, ...
    cfg.i_target_deg+cfg.moduleInclinationOffset_deg,raan+cfg.moduleRAANOffset_deg, ...
    u0+cfg.phaseAngle_deg+rad2deg(nModule*dt));
transfer = FABER_lambert_handoff(mu,Re,taxi,moduleAtTime,cfg.x0_client, ...
    cfg.outboundDepartureEpoch_s,period,rev,cfg.outboundLambertTOF_s, ...
    cfg.outboundMinimumPerigee_km);
if isempty(transfer)
    error('FABER:OutboundLambertNoSolution', ...
        'No outbound Lambert arc passed the perigee and endpoint checks. Adjust duration, revolutions or geometry.');
end

transfer.moduleArrivalECI = transfer.referenceArrivalECI;
transfer.moduleDepartureECI = moduleAtTime(0);
transfer.atlanteDepartureECI = 1000*atlante;
end

function x = circularState(mu,r,inc,raan,u)
radial = [cosd(raan);sind(raan);0];
transverse = [-sind(raan)*cosd(inc); cosd(raan)*cosd(inc); sind(inc)];
x = [r*(cosd(u)*radial+sind(u)*transverse); ...
     sqrt(mu/r)*(-sind(u)*radial+cosd(u)*transverse)];
end
