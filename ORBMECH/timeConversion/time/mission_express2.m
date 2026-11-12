function [t1_refined_mjd, t2_refined_mjd, dv_opt, DV, depT_mjd, arrT_mjd] = ...
    mission_express2(date1, date2, idx_dep, idx_arr, ndep, narr, dv_MIN)

% mission_express2: grid search + optional refinement for interplanetary transfer
% Times are handled in MJD2000 [days] for ephemerides and optimization.
%
% INPUT:
%   date1.dep1 = [Y M D], earliest departure
%   date1.dep2 = [Y M D], latest   departure
%   date2.arr1 = [Y M D], earliest arrival
%   date2.arr2 = [Y M D], latest   arrival
%   idx_dep, idx_arr: planet indices for uplanet
%   ndep, narr: grid resolution
%   dv_MIN: max launcher Δv1 (0 => no launcher constraint)
%
% OUTPUT:
%   t1_refined_mjd, t2_refined_mjd: refined best times (MJD2000)
%   dv_opt: refined total Δv
%   DV: grid of total Δv [narr x ndep]
%   depT_mjd, arrT_mjd: departure/arrival grids (MJD2000)

%% 1. Time windows -> MJD2000 (days)
mjd_init_d  = date2mjd2000([date1.dep1 0 0 0]); % [Y M D h m s]
mjd_final_d = date2mjd2000([date1.dep2 0 0 0]);

mjd_init_a  = date2mjd2000([date2.arr1 0 0 0]);
mjd_final_a = date2mjd2000([date2.arr2 0 0 0]);

%% 2. Build grids in MJD2000 days
n_dep = ndep;
n_arr = narr;

depT_mjd = linspace(mjd_init_d, mjd_final_d, n_dep);  % 1 x n_dep
arrT_mjd = linspace(mjd_init_a, mjd_final_a, n_arr);  % 1 x n_arr

DV       = NaN(n_arr, n_dep);  % total Δv [km/s]
TOF_days = NaN(n_arr, n_dep);  % time of flight [days]
DV1      = NaN(n_arr, n_dep);  % departure Δv1 [km/s]

%% 3. Double loop over grid
for i = 1:n_dep
    t1_mjd = depT_mjd(i);

    % Ephemerides at departure
    [kepD, mu_S] = uplanet(t1_mjd, idx_dep);
    [rD, vD] = kep2car(kepD(1), kepD(2), kepD(3), ...
                       kepD(4), kepD(5), kepD(6), mu_S);
    rD = rD(:);  vD = vD(:);

    for j = 1:n_arr
        t2_mjd = arrT_mjd(j);

        % Arrival must be after departure
        if t2_mjd <= t1_mjd
            continue;
        end

        dt_days = t2_mjd - t1_mjd;
        dt_s    = dt_days * 86400;  % Lambert needs seconds

        % Ephemerides at arrival
        [kepA, ~] = uplanet(t2_mjd, idx_arr);
        [rA, vA] = kep2car(kepA(1), kepA(2), kepA(3), ...
                           kepA(4), kepA(5), kepA(6), mu_S);
        rA = rA(:);  vA = vA(:);

        % Lambert (direct, 0 rev)
        [~,~,~,ERROR,VI,VF] = lambertMR(rD, rA, dt_s, mu_S, 0, 0, 0);
        if ERROR ~= 0
            continue;
        end

        VI = VI(:); VF = VF(:);

        dv1 = norm(VI - vD);
        dv2 = norm(vA - VF);

        DV1(j,i)      = dv1;
        DV(j,i)       = dv1 + dv2;
        TOF_days(j,i) = dt_days;
    end
end

%% 4. Unconstrained grid minimum
[dv_min, idx2] = min(DV(:), [], 'omitnan');
[j_best, i_best] = ind2sub(size(DV), idx2);

t1_best_mjd   = depT_mjd(i_best);
t2_best_mjd   = arrT_mjd(j_best);
tof_best_days = t2_best_mjd - t1_best_mjd;

% For possible unconstrained refinement in days:
t1_best_s = t1_best_mjd * 86400;
t2_best_s = t2_best_mjd * 86400;
tof_best_s = t2_best_s - t1_best_s; %#ok<NASGU> (kept only if you want it)

%% 5. Launcher constraint on Δv1 (grid level)
if dv_MIN ~= 0
    mask_ok = DV1 <= dv_MIN;
    DV_constrained = DV;
    DV_constrained(~mask_ok) = NaN;

    [dv_min_constr, idxc] = min(DV_constrained(:), [], 'omitnan');
    if isempty(dv_min_constr) || isnan(dv_min_constr)
        error('No feasible transfers satisfying Δv1 <= %.2f km/s in this window.', dv_MIN);
    end

    [jc, ic] = ind2sub(size(DV_constrained), idxc);

    t1_constr_mjd   = depT_mjd(ic);
    t2_constr_mjd   = arrT_mjd(jc);
    tof_constr_days = t2_constr_mjd - t1_constr_mjd;

    planetNames = {"Mercury","Venus","Earth","Mars","Jupiter","Saturn","Uranus","Neptune"};
    thisPlanetName = planetNames{idx_arr};

    fprintf('\n===== Step 4: With launcher constraint Δv1 ≤ %.2f km/s =====\n', dv_MIN);
    fprintf('Target planet: %s  (index %d)\n', thisPlanetName, idx_arr);
    fprintf('Constrained Δv_min: %.6f km/s\n', dv_min_constr);
    fprintf('TOF:                %.3f days\n', tof_constr_days);

    depDate = mjd20002date(t1_constr_mjd);
    arrDate = mjd20002date(t2_constr_mjd);

    fprintf('Departure: %04d-%02d-%02d %02d:%02d:%05.2f\n', ...
        depDate(1), depDate(2), depDate(3), depDate(4), depDate(5), depDate(6));
    fprintf('Arrival:   %04d-%02d-%02d %02d:%02d:%05.2f\n', ...
        arrDate(1), arrDate(2), arrDate(3), arrDate(4), arrDate(5), arrDate(6));
    fprintf('============================================================\n\n');
else
    dv_min_constr = NaN;
    t1_constr_mjd = t1_best_mjd;
    t2_constr_mjd = t2_best_mjd;
    tof_constr_days = tof_best_days;
end

%% 6. FMINCON refinement in DAYS: x = [t1_mjd; t2_mjd]
if dv_MIN ~= 0
    % Constrained refinement
    x0 = [t1_constr_mjd; t2_constr_mjd];

    lb = [mjd_init_d; mjd_init_a];
    ub = [mjd_final_d; mjd_final_a];

    obj     = @(x) dv_total_from_t1t2_days(x(1), x(2), mjd_init_a, mjd_final_a, idx_arr, idx_dep);
    nonlcon = @(x) launcher_nonlcon_days(x, dv_MIN, mjd_init_a, mjd_final_a, idx_arr, idx_dep);

    opts = optimoptions('fmincon', ...
        'Display','iter', ...
        'Algorithm','sqp');

    [x_opt, dv_opt] = fmincon(obj, x0, [], [], [], [], lb, ub, nonlcon, opts);

    t1_refined_mjd = x_opt(1);
    t2_refined_mjd = x_opt(2);
else
    % Unconstrained refinement (if you want it)
    x0 = [t1_best_mjd; t2_best_mjd];

    lb = [mjd_init_d; mjd_init_a];
    ub = [mjd_final_d; mjd_final_a];

    obj = @(x) dv_total_from_t1t2_days(x(1), x(2), mjd_init_a, mjd_final_a, idx_arr, idx_dep);

    opts = optimoptions('fmincon', ...
        'Display','iter', ...
        'Algorithm','sqp');

    [x_opt, dv_opt] = fmincon(obj, x0, [], [], [], [], lb, ub, [], opts);

    t1_refined_mjd = x_opt(1);
    t2_refined_mjd = x_opt(2);
end

tof_refined_days = t2_refined_mjd - t1_refined_mjd;
tof_refined_s    = tof_refined_days * 86400;

%% 7. PRINT SUMMARY & RETURN
fprintf("\n===== Refined Solution (Step 6) =====\n");
fprintf("Grid Δv_min:   %.6f km/s\n", dv_min);
fprintf("Refined Δv:    %.6f km/s\n", dv_opt);
fprintf("Refined TOF:   %.3f days\n", tof_refined_days);

depDateRef = mjd20002date(t1_refined_mjd);
arrDateRef = mjd20002date(t2_refined_mjd);

fprintf("Refined Departure: %04d-%02d-%02d %02d:%02d:%05.2f\n", ...
    depDateRef(1), depDateRef(2), depDateRef(3), depDateRef(4), depDateRef(5), depDateRef(6));
fprintf("Refined Arrival:   %04d-%02d-%02d %02d:%02d:%05.2f\n", ...
    arrDateRef(1), arrDateRef(2), arrDateRef(3), arrDateRef(4), arrDateRef(5), arrDateRef(6));
fprintf("=====================================\n\n");

fprintf('Unconstrained grid Δv_min: %.6f km/s\n', dv_min);
if ~isnan(dv_min_constr)
    fprintf('Constrained grid Δv_min:   %.6f km/s\n', dv_min_constr);
end

% OUTPUTS
DV = DV; %#ok<NASGU>   % (DV already in workspace; kept to match signature)

end  % end of mission_express2


%% ===== Helper: DV total from [t1_mjd, t2_mjd] in DAYS ====================
function dv = dv_total_from_t1t2_days(t1_mjd, t2_mjd, mjd_init_a, mjd_final_a, targetIdx, init_idx)

    % Basic feasibility checks (arrival after dep and inside window)
    if t2_mjd <= t1_mjd || t2_mjd < mjd_init_a || t2_mjd > mjd_final_a
        dv = 1e6;
        return;
    end

    dt_days = t2_mjd - t1_mjd;
    dt_s    = dt_days * 86400;

    % Departure planet
    [kepD, mu_S] = uplanet(t1_mjd, init_idx);
    [rD, vD] = kep2car(kepD(1),kepD(2),kepD(3), ...
                       kepD(4),kepD(5),kepD(6), mu_S);
    rD = rD(:); vD = vD(:);

    % Arrival planet
    [kepA, ~] = uplanet(t2_mjd, targetIdx);
    [rA, vA] = kep2car(kepA(1),kepA(2),kepA(3), ...
                       kepA(4),kepA(5),kepA(6), mu_S);
    rA = rA(:); vA = vA(:);

    [~,~,~,ERROR,VI,VF] = lambertMR(rD, rA, dt_s, mu_S, 0, 0, 0);

    if ERROR ~= 0 || any(~isfinite(VI)) || any(~isfinite(VF))
        dv = 1e5 + 1e-4*dt_s;
        return;
    end

    VI = VI(:); VF = VF(:);

    dv = norm(VI - vD) + norm(vA - VF);
end


%% ===== Helper: launcher constraint in DAYS ===============================
function [c, ceq] = launcher_nonlcon_days(x, dv1_max, mjd_init_a, mjd_final_a, idx_arr, idx_dep)

    t1_mjd = x(1);
    t2_mjd = x(2);

    % keep arrival inside window and after departure
    if t2_mjd <= t1_mjd || t2_mjd < mjd_init_a || t2_mjd > mjd_final_a
        c   = 1;
        ceq = [];
        return;
    end

    dt_days = t2_mjd - t1_mjd;
    dt_s    = dt_days * 86400;

    % Departure
    [kepD, mu_S] = uplanet(t1_mjd, idx_dep);
    [rD, vD] = kep2car(kepD(1),kepD(2),kepD(3), ...
                       kepD(4),kepD(5),kepD(6), mu_S);
    rD = rD(:); vD = vD(:);

    % Arrival (position only needed for Lambert)
    [kepA, ~] = uplanet(t2_mjd, idx_arr);
    [rA, ~] = kep2car(kepA(1),kepA(2),kepA(3), ...
                      kepA(4),kepA(5),kepA(6), mu_S);
    rA = rA(:);

    [~,~,~,ERR,VI,~] = lambertMR(rD, rA, dt_s, mu_S, 0, 0, 0);
    if ERR ~= 0 || any(~isfinite(VI))
        c   = 1;
        ceq = [];
        return;
    end

    VI = VI(:);
    dv1 = norm(VI - vD);

    % inequality constraint: dv1 <= dv1_max → dv1 - dv1_max ≤ 0
    c   = dv1 - dv1_max;
    ceq = [];
end
