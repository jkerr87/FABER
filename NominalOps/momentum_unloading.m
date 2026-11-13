function DV_mom_unloading = momentum_unloading(mass_idle, rho, altitude_km, area_drag, Cd, offset_cg_cp, Cr, delta_I, lifetime_years, arm_thrusters)

% Constants
Irr = 1367;
c = 299792458;
mu_earth = astroConstants(13);
R_earth = astroConstants(23);
v_rel = sqrt(398600 / (6371 + altitude_km)) * 1000;

% Disturbances
% Drag
F_drag = 0.5 * rho * v_rel^2 * Cd * area_drag;
DV_mom_unloading.T_drag = F_drag * offset_cg_cp;            

% SRP
F_srp = (Irr / c * Cr * area_drag);
DV_mom_unloading.T_srp = F_srp * offset_cg_cp;              

% GG
r_orbit = (R_earth + altitude_km);
n_orbit = sqrt(mu_earth / r_orbit^3);
DV_mom_unloading.T_gg = 3/2 * n_orbit^2 * delta_I;          

% Magnetic desturbance
D_mag = 0.01 * mass_idle;                    % Dipolo residuo: 0.01 Am^2 per kg
B_field = (2 * 7.96e15) / ((r_orbit * 1000)^3);       % Campo magnetico terrestre (modello base)
DV_mom_unloading.T_mag = D_mag * B_field;                   

% Total torque and momentum
DV_mom_unloading.T_total_disturbance = DV_mom_unloading.T_drag + DV_mom_unloading.T_srp + DV_mom_unloading.T_gg + DV_mom_unloading.T_mag;            % Worst-case: sum of all contribute
H_dump_annual = DV_mom_unloading.T_total_disturbance * 365.25 * 24 * 3600; 
DV_mom_unloading.H_dump_lifetime = H_dump_annual * lifetime_years;

% Total impulse
Impulse_tot_lifetime = DV_mom_unloading.H_dump_lifetime / arm_thrusters;

% Delta_V
DV_mom_unloading.DV = Impulse_tot_lifetime / mass_idle;      % Approximation for small variation of mass

end

