function PH = phasing_manoeuvre(h_base, h_launcher_nominal, dh_err, phase_deg, N_rev, N_modules)
    % Constants
    mu = astroConstants(13);
    R_E = astroConstants(23);
    
    % Corner altitudes
    h_corners = h_launcher_nominal + dh_err * [1; 1; -1; -1];
    
    % Initialization
    DV_phasing_per_cycle = zeros(4, 1);
    DV_phase_1 = zeros(4, 1);
    DV_phase_2 = zeros(4, 1);
    
    % Time Initialization
    time_phase_1 = zeros(4, 1);
    time_phase_2 = zeros(4, 1);
    time_cycle = zeros(4, 1);
    
    % Corner variation
    for j = 1:4
        h_launcher = h_corners(j);
    
        % Phasing 1: rendezvous between Taxi and modules
        r_launcher = R_E + h_launcher;
        T_ref_1 = 2*pi * sqrt(r_launcher^3 / mu);
        T_phasing_1 = T_ref_1 * (1 - phase_deg / (360 * N_rev));
        a_phasing_1 = (mu * (T_phasing_1 / (2*pi))^2)^(1/3);
    
        v_ref_1 = sqrt(mu / r_launcher);
        v_phasing_1 = sqrt(mu * (2/r_launcher - 1/a_phasing_1));
    
        % Delta_V1
        DV_phase_1(j) = 2000 * abs(v_phasing_1 - v_ref_1);
        
        % Time spent in Phasing 1
        time_phase_1(j) = N_rev * T_phasing_1;
    
        % Phasing 2: rendezvous between Taxi and Structure
        r_base = R_E + h_base;
        T_ref_2 = 2*pi * sqrt(r_base^3 / mu);
        T_phasing_2 = T_ref_2 * (1 - phase_deg / (360 * N_rev));
        a_phasing_2 = (mu * (T_phasing_2 / (2*pi))^2)^(1/3);
    
        v_ref_2 = sqrt(mu / r_base);
        v_phasing_2 = sqrt(mu * (2/r_base - 1/a_phasing_2));
    
        % Delta_V2
        DV_phase_2(j) = 2000 * abs(v_phasing_2 - v_ref_2);
        
        % Time spent in Phasing 2
        time_phase_2(j) = N_rev * T_phasing_2;
    
        % Total phasing Delta_V for 1 complete cycle
        DV_phasing_per_cycle(j) = DV_phase_1(j) + DV_phase_2(j);
        
        % Total time for 1 complete cycle (Phasing 1 + Phasing 2)
        time_cycle(j) = time_phase_1(j) + time_phase_2(j);
    end
    
    % Select the worst corner (the one requiring the highest DeltaV)
    [max_DV_cycle, worst_idx] = max(DV_phasing_per_cycle);
    
    % Results (Delta V)
    PH.worst_h_launcher = h_corners(worst_idx);
    PH.DV_cycle_max = max_DV_cycle;
    PH.DV_phase1_module = DV_phase_1(worst_idx);
    PH.DV_phase2_base = DV_phase_2(worst_idx);
    PH.DV_total_lifetime = max_DV_cycle * N_modules;
    
    % Results (Time)
    PH.time_phase1_min = time_phase_1(worst_idx) / 60;
    PH.time_phase2_min = time_phase_2(worst_idx) / 60;
    PH.time_cycle_min = time_cycle(worst_idx) / 60;
    
    % Tempo totale su 7 anni per tutti i moduli
    PH.time_total_lifetime_hours = PH.time_cycle_min * N_modules / 60;
end