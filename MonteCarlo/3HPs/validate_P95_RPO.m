function checks = validate_P95_RPO(R)
% VALIDATE_P95_RPO
%
% Propagate ONLY the selected P95 complete-retrieval
% design from ASV_RPO_MonteCarlo.
%
% Verifies:
% - Homing 1 speed <= 1.0 m/s
% - Homing 2 speed <= 1.0 m/s
% - Closing speed <= 0.3 m/s
% - Final approach speed <= 0.1 m/s
% - Homing remains outside KOZ
% - Closing enters KOZ within negative V-bar corridor
% - Minimum separation
% - Approximate AZ entry state
%
% Requires one HCW arc in each homing phase.
%
% This is sampled trajectory verification,
% not continuous-time certification.

cfg = R.cfg;
T = R.trialTable;

%% SELECT COMPLETE-CYCLE P95 DESIGN

[~,idx] = min(abs( ...
    T.CycleDV_mps - R.P95_DV_cycle_mps));

fprintf('\n====================================\n');
fprintf('P95 RPO TRAJECTORY VERIFICATION\n');
fprintf('====================================\n');
fprintf('Selected trial: %d\n',idx);
fprintf('Cycle DV: %.6f m/s\n',T.CycleDV_mps(idx));
fprintf('Cycle TOF: %.3f min\n',T.CycleTOF_min(idx));

n = sqrt(cfg.mu_km3s2 / ...
    (cfg.Re_km+cfg.h_km)^3);

N = 501;  % Verification sampling resolution

names = ["Homing 1";"Homing 2"; ...
         "Closing";"Final approach"];

speedLimit = [ ...
    cfg.maxSpeed_homing_mps; ...
    cfg.maxSpeed_closing_mps; ...
    cfg.maxSpeed_controlled_mps; ...
    cfg.maxSpeed_final_mps];

rows = {};
allPassed = true;

figure('Color','w','Name','P95 RPO Speed Verification');
tiledlayout(2,1);

for vehicle = 1:2

    if vehicle == 1

        vehicleName = 'Client';

        x0 = [ ...
            T.ClientX0_m(idx);
            T.ClientY0_m(idx);
            T.ClientZ0_m(idx);
            T.ClientVx0_mps(idx);
            T.ClientVy0_mps(idx);
            T.ClientVz0_mps(idx)];

        hp1 = T.ClientHP1_m(idx);
        hp2 = T.ClientHP2_m(idx);
        hp3 = T.ClientHP3_m(idx);

        t1 = T.ClientHoming1_min(idx)*60;
        t2 = T.ClientHoming2_min(idx)*60;
        tc = T.ClientClosing_min(idx)*60;
        tf = T.ClientFinal_min(idx)*60;

        Nh = T.ClientN_Homing1(idx);
        Nc = T.ClientN_Homing2(idx);

    else

        vehicleName = 'ATLANTE';

        x0 = [ ...
            T.HostX0_m(idx);
            T.HostY0_m(idx);
            T.HostZ0_m(idx);
            T.HostVx0_mps(idx);
            T.HostVy0_mps(idx);
            T.HostVz0_mps(idx)];

        hp1 = T.HostHP1_m(idx);
        hp2 = T.HostHP2_m(idx);
        hp3 = T.HostHP3_m(idx);

        t1 = T.HostHoming1_min(idx)*60;
        t2 = T.HostHoming2_min(idx)*60;
        tc = T.HostClosing_min(idx)*60;
        tf = T.HostFinal_min(idx)*60;

        Nh = T.HostN_Homing1(idx);
        Nc = T.HostN_Homing2(idx);
    end

    assert(Nh == 1 && Nc == 1, ...
        ['Verification requires one arc per ', ...
         'homing phase.']);

    rf = cfg.final_m(1);

    X1 = [0;-hp1;0;0;0;0];
    X2 = [0;-hp2;0;0;0;0];
    X3 = [0;-hp3;0;0;0;0];

    %% RECONSTRUCT HCW TRAJECTORIES

    [r1,v1] = hcwCoast(x0,X1,t1,n,N);
    [r2,v2] = hcwCoast(X1,X2,t2,n,N);
    [r3,v3] = hcwCoast(X2,X3,tc,n,N);

    %% FINAL QUINTIC APPROACH

    s = linspace(0,1,N);
    p = 10*s.^3-15*s.^4+6*s.^5;
    pd = (30*s.^2-60*s.^3+30*s.^4)/tf;

    dr = [0;-(rf-hp3);0];

    r4 = repmat(X3(1:3),1,N)+dr*p;
    v4 = dr*pd;

    RR = {r1,r2,r3,r4};
    VV = {v1,v2,v3,v4};
    durations = [t1,t2,tc,tf];

    %% VERIFY EACH PHASE

    fprintf('\n--- %s ---\n',vehicleName);

    nexttile(vehicle);
    hold on; grid on;

    offset = 0;

    for k = 1:4

        r = RR{k};
        v = VV{k};

        range = vecnorm(r,2,1);
        speed = vecnorm(v,2,1);

        vmax = max(speed);
        rmin = min(range);

        speedOK = vmax <= speedLimit(k)+1e-8;

        if k <= 2

            % Homing must remain outside KOZ
            zoneOK = all(range >= cfg.KOZ_m);

        elseif k == 3

            % Closing may enter KOZ only through
            % authorized negative V-bar corridor.

            inside = range < cfg.KOZ_m;

            lateral = hypot(r(1,inside),r(3,inside));

            zoneOK = all( ...
                lateral <= cfg.corridorHalfWidth_m) && ...
                all(r(2,inside)<0) && ...
                all(range >= rf);

        else

            % Nominal final approach follows V-bar.
            lateral = hypot(r(1,:),r(3,:));

            zoneOK = all( ...
                lateral <= cfg.corridorHalfWidth_m) && ...
                all(range >= rf-1e-8);
        end

        passed = speedOK && zoneOK;
        allPassed = allPassed && passed;

        fprintf( ...
            '%-15s vmax=%.4f  limit=%.2f  ', ...
            names(k),vmax,speedLimit(k));

        fprintf('minRange=%.2f m  ',rmin);

        if passed
            fprintf('PASS\n');
        else
            fprintf('FAIL');
            if ~speedOK
                fprintf(' [speed]');
            end

            if ~zoneOK
                fprintf(' [KOZ/corridor]');
            end

            fprintf('\n');
        end

        rows(end+1,:) = { ...
            vehicleName,char(names(k)), ...
            vmax,speedLimit(k),rmin, ...
            speedOK,zoneOK,passed};

        time = offset + ...
            linspace(0,durations(k)/60,N);

        plot(time,speed,'LineWidth',1.3, ...
            'DisplayName',char(names(k)));

        offset = time(end);
    end

    xlabel('Transfer time [min] (holds excluded)');
    ylabel('Relative speed [m/s]');
    title([vehicleName ' - selected P95 design']);
    legend('Location','best');

    %% APPROXIMATE AZ CROSSING DURING CLOSING

    rangesClosing = vecnorm(r3,2,1);
    entry = find(rangesClosing <= cfg.AZ_m,1,'first');

    if ~isempty(entry)

        vEntry = norm(v3(:,entry));

        fprintf('AZ entry speed: %.4f m/s\n',vEntry);
        fprintf('AZ entry range: %.2f m\n', ...
            rangesClosing(entry));

    else

        fprintf('AZ entry not detected in closing samples.\n');
    end

end

%% RESULTS TABLE

checks = cell2table(rows, ...
    'VariableNames',{ ...
    'Vehicle','Phase','PeakSpeed_mps', ...
    'SpeedLimit_mps','MinRange_m', ...
    'SpeedOK','ZoneOK','Passed'});

fprintf('\n====================================\n');

if allPassed
    fprintf('P95 DESIGN PASSES SAMPLED CHECKS\n');
else
    fprintf('P95 DESIGN FAILS SAMPLED CHECKS\n');
end

fprintf('====================================\n');

disp(checks);

end

%% ================================================================
% HCW COAST PROPAGATION
% ================================================================

function [r,v] = hcwCoast(x0,xf,T,n,N)

Phi = HCW_STM(n,T);

Prr = Phi(1:3,1:3);
Prv = Phi(1:3,4:6);

% Required initial velocity after departure burn
v0plus = Prv \ ...
    (xf(1:3)-Prr*x0(1:3));

Xplus = [x0(1:3);v0plus];

r = zeros(3,N);
v = zeros(3,N);

for k = 1:N

    tk = T*(k-1)/(N-1);

    X = HCW_STM(n,tk)*Xplus;

    r(:,k) = X(1:3);
    v(:,k) = X(4:6);
end

end

%% ================================================================
% ANALYTICAL HCW STM
% ================================================================

function Phi = HCW_STM(n,t)

nt = n*t;
c = cos(nt);
s = sin(nt);

Phi = zeros(6,6);

Phi(1:3,1:3) = [ ...
    4-3*c,0,0; ...
    6*(s-nt),1,0; ...
    0,0,c];

Phi(1:3,4:6) = [ ...
    s/n,2*(1-c)/n,0; ...
    -2*(1-c)/n,(4*s-3*nt)/n,0; ...
    0,0,s/n];

Phi(4:6,1:3) = [ ...
    3*n*s,0,0; ...
    6*n*(c-1),0,0; ...
    0,0,-n*s];

Phi(4:6,4:6) = [ ...
    c,2*s,0; ...
    -2*s,4*c-3,0; ...
    0,0,c];

end
