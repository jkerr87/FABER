% Run the three-holding-point V-bar study and check its P95 trajectory.
clear; clc; close all;

cfg = ASV_RPO_MC_config();

cfg.hp2_m = [400 400];
cfg.methods = {'Vbar'};

cfg.N_homing = [1 1];
cfg.N_closing = [1 1];

cfg.N = 1000;
cfg.maxAttemptsFactor = 30;

cfg.makePlots = false;
cfg.saveCSV = false;

R = ASV_RPO_MonteCarlo(cfg);

checks = validate_P95_RPO(R);

disp(R.summaryTable);
disp(R.phaseTable);
