
function sk = station_keeping(mass_kg, rho, area_m2, Cd, ...
                              altitude_km, lifetime_years)
% Preliminary orbital drag-compensation Delta-V.
%
% This is NOT a relative-formation-keeping controller.
% Assumes constant atmospheric density and orbital conditions.
%
% Units:
% mass_kg        [kg]
% rho            [kg/m^3]
% area_m2        [m^2]
% altitude_km    [km]
% lifetime_years [years]
% Delta-V        [m/s]

mu = 398600.4418;   % [km^3/s^2]
Re = 6378.1363;     % [km]

v_rel = sqrt(mu/(Re + altitude_km))*1000;

a_drag = 0.5*rho*v_rel^2*Cd*area_m2/mass_kg;

DV_annual = a_drag*365.25*24*3600;
DV_lifetime = DV_annual*lifetime_years;

sk = struct();

sk.DV_SK = DV_lifetime;    % Existing output preserved
sk.DV_annual_mps = DV_annual;
sk.DV_lifetime_mps = DV_lifetime;
sk.dragAcceleration_mps2 = a_drag;
sk.orbitalVelocity_mps = v_rel;

end
