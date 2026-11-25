function mcCfg = FABER_phaseA_rpo_config(cfg,mu,Re,outbound)
% Shared nine-phase setup used by Phase A and the sortie study.
mcCfg = ASV_RPO_MC_config();

% Pass the module-relative Lambert handoff through without resampling its state.
mcCfg.mu_km3s2 = mu;
mcCfg.Re_km = Re;
mcCfg.h_km = cfg.h_target_km;
mcCfg.x0_client = outbound.arrivalLVLH;
mcCfg.arrivalEpochClient_s = outbound.arrivalEpoch_s;

mcCfg.clientH_km = norm(outbound.moduleArrivalECI(1:3))/1000-Re;
mcCfg.arrivalStatesClient = repmat(outbound.arrivalLVLH,1,cfg.N_RPO_MC);
mcCfg.x0_host = cfg.x0_atlante(:);
mcCfg.N = cfg.N_RPO_MC;
mcCfg.seed = cfg.RPO_seed_client;

% Use the Phase A reference approach and impulsive burn-duration limit.
mcCfg.methodsH2 = {'Vbar'};
mcCfg.methodsC1 = {'Rbar'};
mcCfg.methodsC2 = {'GeneralHCW'};
mcCfg.gateModes = {'FlyThrough'};
mcCfg.homing1_min = [30 100];
mcCfg.maxBurnFraction = 0.10;
% Older saved Phase A configurations inherit the standalone screen settings.
input = ["rpoClientGeometry","rpoAtlanteGeometry","rpoTerminalCorridorHalfWidth_m", ...
    "rpoMaxPropellantCycle_kg","rpoMaxDVPerCycle_mps","rpoMaxBurnDuration_s","rpoMaxPeakAccel_mps2"];
setting = ["clientGeometry","hostGeometry","terminalCorridorHalfWidth_m", ...
    "maxPropellantCycle_kg","maxDV_perCycle_mps","maxBurnDuration_s","maxPeakAccel_mps2"];
for j = 1:numel(input)
    if isfield(cfg,input(j)), mcCfg.(setting(j)) = cfg.(input(j)); end
end

% Deduct outbound far-field propellant before sizing the client approach.
outboundDV = norm(outbound.departureBurn_mps)+norm(outbound.arrivalBurn_mps);
mcCfg.massTaxi_kg = cfg.massTaxi_kg*exp(-outboundDV/(9.80665*cfg.farFieldIsp_s));
mcCfg.returnHandoff = @(client,mass) FABER_return_lambert(mu,Re,cfg,outbound,client,mass);
mcCfg.payloadMass_kg = cfg.payloadMass_kg;
mcCfg.makePlots = false;
mcCfg.saveCSV = false;
end
