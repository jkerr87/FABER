% Main to run monte carlo
clear; clc; close all;


cfg = ASV_RPO_MC_config();

cfg.N = 1000;
cfg.maxAttemptsFactor = 10;
cfg.nTrajectorySamples = 101;

R = ASV_RPO_MonteCarlo(cfg);