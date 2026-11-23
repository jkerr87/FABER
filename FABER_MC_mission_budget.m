function B = FABER_MC_mission_budget(baseline,R,rpoRows,architecture,N,seed)
% Propagate paired close-field cycle results into a lifetime budget.
% Each campaign resamples N_deliveries independent attempts with replacement.
% Unsuccessful campaigns remain in the denominator and have unknown cost.
validateattributes(N,{'numeric'},{'scalar','integer','positive'});
validateattributes(architecture,{'numeric'}, ...
    {'scalar','integer','>=',1,'<=',size(R.valid,2)});
assert(~isempty(rpoRows) && numel(unique(rpoRows))==numel(rpoRows));
assert(all(baseline.Count(rpoRows)==baseline.Count(rpoRows(1))), ...
    'RPO rows must share the same delivery count.');
deliveries = baseline.Count(rpoRows(1));
validateattributes(deliveries,{'numeric'},{'scalar','integer','nonnegative'});
keep = true(height(baseline),1);
keep(rpoRows) = false;
includedDV = baseline.RawLifetime_mps(keep);
baseDV = sum(includedDV,'omitnan');
baseMarginedDV = sum(baseline.MarginedLifetime_mps(keep),'omitnan');
baseTOF = sum(baseline.TOF_lifetime_min(keep),'omitnan');
valid = R.valid(:,architecture);
cycleDV = R.CycleDV_mps(:,architecture);
cycleTOF = R.CycleTOF_min(:,architecture);
assert(~isempty(valid) && all(isfinite(cycleDV(valid))) && ...
    all(isfinite(cycleTOF(valid))));
stream = RandStream('mt19937ar','Seed',seed);
draw = randi(stream,numel(valid),N,deliveries);
feasible = all(reshape(valid(draw),N,deliveries),2);
rawDV = nan(N,1);
marginedDV = nan(N,1);
knownTOF = nan(N,1);
accepted = draw(feasible,:);
acceptedCycleDV = sum(reshape(cycleDV(accepted),size(accepted)),2);
rawDV(feasible) = baseDV + acceptedCycleDV;
marginedDV(feasible) = baseMarginedDV + ...
    FABER_apply_DV_margin(acceptedCycleDV,true);
knownTOF(feasible) = baseTOF + sum(reshape(cycleTOF(accepted),size(accepted)),2);
B.architecture = R.methods(architecture,:);
B.unknownDVRows = baseline.Maneuver(keep & isnan(baseline.RawLifetime_mps));
B.unknownTOFRows = baseline.Maneuver(keep & isnan(baseline.TOF_lifetime_min));
B.cycleAcceptance = mean(valid);
B.campaignAcceptance = B.cycleAcceptance^deliveries;
B.resampledAcceptance = mean(feasible);
knownRawDV_mps = rawDV;
rawDV_mps = rawDV;
knownMarginedDV_mps = marginedDV;
marginedDV_mps = marginedDV;
if ~isempty(B.unknownDVRows)
    rawDV_mps(:) = NaN;
    marginedDV_mps(:) = NaN;
end
B.trials = table(feasible,rawDV_mps,knownRawDV_mps,marginedDV_mps, ...
    knownMarginedDV_mps,knownTOF,'VariableNames',{'RPOFeasible', ...
    'RawDV_mps','KnownRawDV_mps','MarginedDV_mps','KnownMarginedDV_mps','KnownTOF_min'});
if any(valid) || deliveries==0
    % Conditional on each delivery being feasible. Draw directly from
    % accepted cycles so rare campaign success does not erase this estimate.
    pool = find(valid);
    if deliveries==0
        B.P95KnownRawDV_mps = baseDV;
        B.P95KnownMarginedDV_mps = baseMarginedDV;
        B.P95KnownTOF_min = baseTOF;
    else
        conditional = pool(randi(stream,numel(pool),N,deliveries));
        rpoP95 = prctile(sum(reshape(cycleDV(conditional),N,deliveries),2),95);
        B.P95KnownRawDV_mps = baseDV + rpoP95;
        B.P95KnownMarginedDV_mps = baseMarginedDV + ...
            FABER_apply_DV_margin(rpoP95,true);
        B.P95KnownTOF_min = baseTOF + prctile(sum(reshape( ...
            cycleTOF(conditional),N,deliveries),2),95);
    end
else
    B.P95KnownRawDV_mps = NaN;
    B.P95KnownMarginedDV_mps = NaN;
    B.P95KnownTOF_min = NaN;
end
B.P95RawDV_mps = B.P95KnownRawDV_mps;
B.P95MarginedDV_mps = B.P95KnownMarginedDV_mps;
if ~isempty(B.unknownDVRows)
    B.P95RawDV_mps = NaN;
    B.P95MarginedDV_mps = NaN;
end
end
