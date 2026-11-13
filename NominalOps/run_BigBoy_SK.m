% FABER Big Boy station-keeping trade (illustrative assumptions).
% Run this script from the same folder as ASV_BigBoy_StationKeeping.m.
clear; clc;

% Override any default model setting here. All omitted settings retain
% documented illustrative values in the main function.
cfg = struct();
cfg.hostAltitude_m = 457250;   % [m] Reference host orbit
cfg.inclination_deg = 40;      % [deg]
cfg.distances_m = [100 250 500 1000];   % parking distances [m]
cfg.deadbands_m = [5 10 25];            % position deadbands [m]
cfg.missionDays = 7;          % validate first; then try 30 or 90
cfg.step_s = 60;              % halve step for numerical convergence
cfg.saveCSV = true;
cfg.makePlots = true;

R_BB = ASV_BigBoy_StationKeeping(cfg);

% View results:
% R_BB.summaryTable
% R_BB.caseResults{1}.history
% R_BB.caseResults{1}.burnDV_mps
