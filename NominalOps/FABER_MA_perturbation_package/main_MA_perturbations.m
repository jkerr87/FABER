%% FABER MA - separate RPO perturbation and parking station-keeping studies
% Needs FABER_RPO_9phase_MC.m, ASV_RPO_MC_config.m, analytical_manos.m,
% rpo_continuous.m and these new MA functions on the MATLAB path.
clear; clc;

%% 1. Common environmental inputs (REPLACE WITH VEHICLE PROPERTIES)
env=FABER_MA_dynamics('defaults');
env.useJ2=true;
env.useDrag=true;
env.useSRP=false;  % SRP disabled until attitude / sun direction chosen
env.rhoRefAlt_m=680120;   % reference density altitude, m
env.rhoRef_kgm3=1e-12; % EXAMPLE ONLY - run density/solar-cycle sensitivity
env.scaleHeight_m=60000;
% Use actual ballistic properties of BOTH bodies for differential drag.
env.chief.mass_kg=1000;        % illustrative client or host
% NOTE: these values are placeholders to replace with current configurations
env.chief.area_m2=8;
env.chief.Cd=2.2;
env.deputy.mass_kg=3503.59; % Taxi wet mass
env.deputy.area_m2=12;      % illustrative

%% 2. RPO perturbation correction allocation (nominal MC unchanged)
cfg=ASV_RPO_MC_config();
cfg.N=50; cfg.makePlots=false; cfg.saveCSV=false;
R=FABER_RPO_9phase_MC(cfg);
p=struct();
p.env=env;
p.inclination_deg=98.12;
p.selection='selected';       % 'selected', 'pareto', 'all'
p.maxCandidates=50;
p.correctionFraction=0.5;  % prescribed impulse at 50% transfer

p.thrust_N=cfg.thrust_N*cfg.nActiveThrusters*cfg.thrustEfficiency;
p.mass_kg=cfg.massTaxi_kg;
p.Isp_s=cfg.Isp_s;        % update from actual thruster selection
p.maxCorrectionBurnFraction=0.05;
p.holdPositionTolerance_m=0.5;
p.holdVelocityTolerance_mps=0.001;
B=FABER_RPO_perturbation_budget(R,p);

%% 3. Idle holding ellipse for Taxi (separate from short RPO holds)
s=struct();
s.env=env;
s.h_km=680.12; s.inclination_deg=98.12;
s.architecture='ellipse';  % alternative 'vbar'
s.centerAlongTrack_m=-1000;
s.radialAmplitude_m=20;
s.duration_days=1;       % start with 1, then test 7, 30, 90 ...
s.checkInterval_s=5400;
s.correctionTOF_s=1800;
s.positionTolerance_m=10;
s.velocityTolerance_mps=0.01;
s.minClearance_m=300;
s.thrust_N=p.thrust_N;
s.Isp_s=p.Isp_s;
SK=FABER_relative_SK_budget(s);

fprintf('\n--- PHASE A SUMMARY ---\n');
if ~isempty(B.table)
    good=B.table.CorrectionFeasible;
    fprintf('RPO perturbation cases passed: %d / %d attempted\n',sum(good),height(B.table));
    if any(good)
        fprintf('Median additional RPO DV (accepted): %.6f m/s\n',median(B.table.ExtraDV_mps(good)));
    end
end

fprintf('Taxi parking SK: %.6f m/s in %.1f day(s)\n',SK.TotalDV_mps,s.duration_days);
fprintf('Taxi parking corrections: %d | minimum separation: %.1f m\n',...
    SK.Ncorrections,SK.MinSeparation_m);

figure('Color','w');
plot(SK.time_s/86400,SK.positionError_m,'LineWidth',1.2); grid on;
xlabel('Elapsed time [days]'); ylabel('Position error vs. nominal ellipse [m]');
title('FABER relative parking drift and discrete correction events');
