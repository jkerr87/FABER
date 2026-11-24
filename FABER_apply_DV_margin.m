function margined = FABER_apply_DV_margin(raw,stochastic,fraction,floor_mps)
% Apply the delta-V allowance once per powered event.
% NaN is unresolved. An explicit zero has no propulsion allocation.
if nargin < 3
    fraction = 0.05;
end

if nargin < 4
    floor_mps = 10;
end

validateattributes(fraction,{'numeric'},{'scalar','finite','nonnegative'});
validateattributes(floor_mps,{'numeric'},{'scalar','finite','nonnegative'});
assert(isscalar(stochastic) || isequal(size(raw),size(stochastic)));
assert(all(isnan(raw(:)) | (isfinite(raw(:)) & raw(:)>=0)));
margined = raw;
deterministic = raw>0 & ~stochastic;
margined(deterministic) = raw(deterministic) + ...
    max(fraction*raw(deterministic),floor_mps);
sampled = raw>0 & stochastic;
margined(sampled) = 2*raw(sampled);
end
