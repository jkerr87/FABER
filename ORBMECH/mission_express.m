function [t1_refined_mjd, t2_refined_mjd, dv_opt, DV, depT_mjd, arrT_mjd] = mission_express(date1, date2, idx_dep, idx_arr, ndep, narr, dv_MIN)

% in MJD2000 days -- Time, modified Julian day since 01/01/2000, 12:00 noon
mjd_init_d  = date2mjd2000m(date1.dep1(1), date1.dep1(2), date1.dep1(3));
mjd_final_d = date2mjd2000m(date1.dep2(1), date1.dep2(2), date1.dep2(3));

mjd_init_a  = date2mjd2000m(date2.arr1(1), date2.arr1(2), date2.arr1(3));
mjd_final_a = date2mjd2000m(date2.arr2(1), date2.arr2(2), date2.arr2(3));

% in seconds - DEPARTURE
init_d  = mjd_init_d  * 86400;
final_d = mjd_final_d * 86400;

% in seconds - ARRIVAL
init_a  = mjd_init_a  * 86400;
final_a = mjd_final_a * 86400;

%% for loop to calculate all combinations
n_dep = ndep;
n_arr = narr;

% Create a grid of departure and arrival times
depT_s   = linspace(init_d,  final_d,  n_dep);
depT_mjd = linspace(mjd_init_d, mjd_final_d, n_dep);

arrT_s   = linspace(init_a,  final_a,  n_arr);
arrT_mjd = linspace(mjd_init_a, mjd_final_a, n_arr);

% IMPORTANT: rows = arrival, cols = departure
DV       = NaN(n_arr, n_dep);
TOF_days = NaN(n_arr, n_dep);
DV1      = NaN(n_arr, n_dep);   % NEW: injection Δv1

for i = 1:n_dep
    
    % Set initial times
    t1_s   = depT_s(i);
    t1_mjd = depT_mjd(i);

    % get eph at departure
    [kepE, mu_S] = uplanet(t1_mjd, idx_dep);

    % get cartesian elements (kep already in rad)
    [rE, vE] = kep2car(kepE(1), kepE(2), kepE(3), kepE(4), ...
                      kepE(5), kepE(6), mu_S);
    rE = rE(:);  vE = vE(:);

    for j = 1:n_arr
        
        t2_mjd = arrT_mjd(j);
        t2_s   = arrT_s(j);

        % arrival must be after departure
        if t2_mjd <= t1_mjd
            continue;
        end

        [kepM, mu_S] = uplanet(t2_mjd, idx_arr);

        [rM, vM] = kep2car(kepM(1), kepM(2), kepM(3), kepM(4), ...
                          kepM(5), kepM(6), mu_S);
        rM = rM(:);  vM = vM(:);
        
        ToF = (t2_s - t1_s);               % [s]
        TOF_days(j,i) = (t2_mjd - t1_mjd); % [days]

        % calling lambert (Sun-centered!)
        [~,~,~,ERROR,VI,VF] = lambertMR(rE, rM, ToF, mu_S, 0, 0, 0);

        if ERROR ~= 0
            continue;
        end

        VI = VI(:); VF = VF(:);

        % Saving results (rows=arrival j, cols=departure i)
        % Δv1 and Δv2
        dv1 = norm(VI - vE);   % departure impulse
        dv2 = norm(vM - VF);   % arrival impulse

        DV1(j,i) = dv1;        % store Δv1
        DV(j,i) = dv1 + dv2;  % store total Δv;

    end
end

% Analyzing data - getting best delta v
[dv_min, idx2] = min(DV(:), [], 'omitnan');
[j_best, i_best] = ind2sub(size(DV), idx2);  % j_best=arrival index, i_best=departure index

% best times
t1_best_mjd = depT_mjd(i_best);
t1_best_s   = depT_s(i_best);

t2_best_mjd = arrT_mjd(j_best);
t2_best_s   = arrT_s(j_best);

% time of flight
tof_best_mjd = t2_best_mjd - t1_best_mjd;
tof_best_s   = t2_best_s   - t1_best_s;

%% DELTA V CONSTRAINT========================================================
% Launcher constraint mask - ensuring the launcher can put you into orbit
if dv_MIN ~= 0
    mask_ok = DV1 <= dv_MIN;
    
    DV_constrained = DV;
    DV_constrained(~mask_ok) = NaN;   % discard infeasible points
    
    [dv_min_constr, idxc] = min(DV_constrained(:), [], 'omitnan');
    
    if isempty(dv_min_constr) || isnan(dv_min_constr)
        error('No feasible transfers satisfying Δv1 <= %.2f km/s in this window.', dv_MIN);
    end
    
    [jc, ic] = ind2sub(size(DV_constrained), idxc);
    
    t1_constr_mjd = depT_mjd(ic);
    t2_constr_mjd = arrT_mjd(jc);
    tof_constr_days = t2_constr_mjd - t1_constr_mjd;
    
% Resolve planet name from index
    planetNames = {"Mercury","Venus","Earth","Mars","Jupiter","Saturn","Uranus","Neptune"};
    thisPlanetName = planetNames{idx_arr};

    fprintf('\n===== Step 4: With launcher constraint Δv1 ≤ %.2f km/s =====\n', dv_MIN);
    fprintf('Target planet: %s  (index %d)\n', thisPlanetName, idx_arr);
    fprintf('Constrained Δv_min: %.6f km/s\n', dv_min_constr);
    fprintf('TOF:                %.3f days\n', tof_constr_days);

    [YYd,MMd,DDd,hhD,mmD,ssD] = mjd2000_to_date(t1_constr_mjd);
    [YYa,MMa,DDa,hhA,mmA,ssA] = mjd2000_to_date(t2_constr_mjd);

    fprintf('Departure: %04d-%02d-%02d %02d:%02d:%05.2f\n', ...
            YYd,MMd,floor(DDd),hhD,mmD,ssD);
    fprintf('Arrival:   %04d-%02d-%02d %02d:%02d:%05.2f\n', ...
            YYa,MMa,floor(DDa),hhA,mmA,ssA);

    fprintf('============================================================\n\n');
end


%% FMINCON ===============================================================
% seeing if launcher constriant exists
if dv_MIN ~= 0

    TOF_min_days = 60;                 % keep Lambert away from short-TOF junk
    TOF_min_s    = TOF_min_days*86400;

    t1_constr_s  = t1_constr_mjd * 86400;
    t2_constr_s  = t2_constr_mjd * 86400;
    tof_constr_s = t2_constr_s - t1_constr_s;

    x0 = [t1_constr_s; tof_constr_s];

    % bounds on t1 and TOF
    lb = [init_d; TOF_min_s];
    ub = [final_d; final_a - init_d];  % biggest possible TOF in your windows
    
    % objective in [t1_s, TOF_s]
    obj = @(x) dv_total_from_t1tof_idx(x(1), x(2), init_a, final_a, idx_arr, idx_dep);
    nonlcon = @(x) launcher_nonlcon(x, init_a, final_a, dv_MIN, idx_arr, idx_dep);

    opts = optimoptions('fmincon', ...
    'Display','iter', ...
    'Algorithm','sqp'); 
    %opts = optimoptions('fmincon', ...
    %    'Display','iter', ...
    %    'Algorithm','interior-point', ...
    %    'FiniteDifferenceType','central', ...
    %    'FiniteDifferenceStepSize', 8000, ...   % ~2.2 hr step stabilizes gradients
    %    'TolFun',1e-8, 'TolX',1e-10, ...
    %    'MaxFunctionEvaluations', 8000);
    
    [x_opt, dv_opt] = fmincon(obj, x0, [], [], [], [], lb, ub, nonlcon, opts);

    t1_refined_s  = x_opt(1);
    tof_refined_s = x_opt(2);
    t2_refined_s  = t1_refined_s + tof_refined_s;
    
    t1_refined_mjd = t1_refined_s/86400;
    t2_refined_mjd = t2_refined_s/86400;

else
    % no launcher constraint case
    TOF_min_days = 60;                 % keep Lambert away from short-TOF junk
    TOF_min_s    = TOF_min_days*86400;
    
    % initial guess from grid
    x0 = [t1_best_s; tof_best_s];
    
    % bounds on t1 and TOF
    lb = [init_d; TOF_min_s];
    ub = [final_d; final_a - init_d];  % biggest possible TOF in your windows
    
    % objective in [t1_s, TOF_s]
    obj = @(x) dv_total_from_t1tof_idx(x(1), x(2), init_a, final_a, idx_arr, idx_dep);
    
    opts = optimoptions('fmincon', ...
        'Display','iter', ...
        'Algorithm','interior-point', ...
        'FiniteDifferenceType','central', ...
        'FiniteDifferenceStepSize', 2000, ...   % ~2.2 hr step stabilizes gradients
        'TolFun',1e-8, 'TolX',1e-10, ...
        'MaxFunctionEvaluations', 8000);
    
    [x_opt, dv_opt] = fmincon(obj, x0, [], [], [], [], lb, ub, [], opts);
    
    t1_refined_s  = x_opt(1);
    tof_refined_s = x_opt(2);
    t2_refined_s  = t1_refined_s + tof_refined_s;
    
    t1_refined_mjd = t1_refined_s/86400;
    t2_refined_mjd = t2_refined_s/86400;

end
    %% RETURN RESULTS =====================================================
fprintf("\n===== Refined Solution (Step 6) =====\n");
fprintf("Grid Δv_min:   %.6f km/s\n", dv_min);
fprintf("Refined Δv:    %.6f km/s\n", dv_opt);
fprintf("Refined TOF:   %.3f days\n", tof_refined_s/86400);

[depY,depM,depD,depH,depMin,depSec] = mjd2000_to_date(t1_refined_mjd);
[arrY,arrM,arrD,arrH,arrMin,arrSec] = mjd2000_to_date(t2_refined_mjd);

fprintf("Refined Departure: %04d-%02d-%02d %02d:%02d:%05.2f\n", depY,depM,floor(depD),depH,depMin,depSec);
fprintf("Refined Arrival:   %04d-%02d-%02d %02d:%02d:%05.2f\n", arrY,arrM,floor(arrD),arrH,arrMin,arrSec);
fprintf("=====================================\n\n");

fprintf('Unconstrained grid Δv_min: %.6f km/s\n', dv_min);
fprintf('Constrained grid Δv_min:   %.6f km/s\n', dv_min_constr);


end


function dv = dv_total_from_t1tof_idx(t1_s, tof_s, init_a_s, final_a_s, targetIdx, init_idx)

    t2_s = t1_s + tof_s;

    % keep arrival inside window
    if t2_s < init_a_s || t2_s > final_a_s
        dv = 1e6 + 1e-2*abs(t2_s - min(max(t2_s,init_a_s),final_a_s));
        return;
    end

    t1_mjd = t1_s/86400;
    t2_mjd = t2_s/86400;

    % Earth at t1
    [kepE, mu_S] = uplanet(t1_mjd, init_idx);
    [rE, vE] = kep2car(kepE(1),kepE(2),kepE(3), ...
                       kepE(4),kepE(5),kepE(6), mu_S);
    rE = rE(:); vE = vE(:);

    % target planet at t2
    [kepT, ~] = uplanet(t2_mjd, targetIdx);
    [rT, vT] = kep2car(kepT(1),kepT(2),kepT(3), ...
                       kepT(4),kepT(5),kepT(6), mu_S);
    rT = rT(:); vT = vT(:);

    [~,~,~,ERROR,VI,VF] = lambertMR(rE, rT, tof_s, mu_S, 0,0,0);

    if ERROR ~= 0 || any(~isfinite(VI)) || any(~isfinite(VF))
        dv = 1e5 + 1e-4*tof_s;
        return;
    end

    VI = VI(:); VF = VF(:);

    dv = norm(VI - vE) + norm(vT - VF);
end


function [c, ceq] = launcher_nonlcon(x, init_a_s, final_a_s, dv1_max, idx_arr, idx_dep)
    t1_s  = x(1);
    tof_s = x(2);
    t2_s  = t1_s + tof_s;

    % keep arrival inside window (same as your objective penalty)
    if t2_s < init_a_s || t2_s > final_a_s
        c   = 1;    % infeasible
        ceq = [];
        return;
    end

    t1_mjd = t1_s/86400;
    t2_mjd = t2_s/86400;

    [kepD, mu_S] = uplanet(t1_mjd, idx_dep);
    [rD, vD] = kep2car(kepD(1),kepD(2),kepD(3), ...
                       kepD(4),kepD(5),kepD(6), mu_S);
    rD = rD(:); vD = vD(:);

    [kepA, ~] = uplanet(t2_mjd, idx_arr);
    [rA, vA] = kep2car(kepA(1),kepA(2),kepA(3), ...
                       kepA(4),kepA(5),kepA(6), mu_S);
    rA = rA(:); vA = vA(:);

    [~,~,~,ERR,VI,VF] = lambertMR(rD, rA, tof_s, mu_S, 0,0,0);
    if ERR ~= 0
        c   = 1;    % mark infeasible
        ceq = [];
        return;
    end

    VI = VI(:);

    dv1 = norm(VI - vD);

    % inequality constraint: dv1 - dv1_max ≤ 0
    c   = dv1 - dv1_max;
    ceq = [];
end
