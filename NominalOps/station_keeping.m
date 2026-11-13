function sk_results = station_keeping(mass_dry, rho, area_drag, Cd, altitude_km, lifetime_years)

% Constants
v_rel = sqrt(398600 / (6371 + altitude_km)) * 1000;

% Drag DV Calculations
acc_drag_nominal = 0.5 * rho * v_rel^2 * Cd * area_drag / mass_dry;
deltaV_drag_annual = acc_drag_nominal * 365.25 * 24 * 3600;
deltaV_drag_lifetime = deltaV_drag_annual * lifetime_years;

% Result struct
sk_results.DV_SK = deltaV_drag_lifetime;

end