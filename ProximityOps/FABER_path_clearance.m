function clearance_m = FABER_path_clearance(position_m, geometry)
% Sampled clearance of a 3-by-N LVLH trajectory from configured obstacles.
% Use existing sphere/capsule geometry and expand by vehicle and margin.
if isempty(fieldnames(geometry))
    clearance_m = Inf;
    return;
end

if ~isfield(geometry,'extraRadius_m')
    geometry.extraRadius_m = 0;
end

validateattributes(geometry.extraRadius_m,{'numeric'}, ...
    {'scalar','real','finite','nonnegative'});
clearance_m = Inf;
for k = 1:size(position_m,2)
    clearance_m = min(clearance_m,FABER_clearance_m( ...
        position_m(:,k),geometry,geometry.extraRadius_m));
end
end
