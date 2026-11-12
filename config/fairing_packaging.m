function pack = fairing_packaging(N_largeModules, N_smallElements, ...
    smallElementSize, packingBox, packingFactor)
% FAIRING_PACKAGING
% Estimates launches and small-element packing capacity.
%
% Units:
% smallElementSize [L;W;H] = m
% packingBox      [L;W;H] = m
%
% Assumptions:
% - Rectangular payload packing proxy
% - Identical small elements
% - Orthogonal packing with possible axis rotation
% - Packing factor reduces the ideal geometric count
% - One large module per launch
% - No mass or structural-load constraints

smallDims = sort(smallElementSize(:));
bayDims = sort(packingBox(:));

assert(numel(smallDims)==3 && numel(bayDims)==3, ...
    'Dimensions must be 3-element vectors.');

assert(all(smallDims>0) && all(bayDims>0), ...
    'All dimensions must be positive.');

assert(packingFactor>0 && packingFactor<=1, ...
    'Packing factor must be between 0 and 1.');

% Geometric packing capacity
N_fit = floor(bayDims ./ smallDims);
N_geometric = prod(N_fit);

% Apply assumed utilization
N_smallPerFairing = floor(packingFactor*N_geometric);

if N_smallElements>0 && N_smallPerFairing<1
    error('No small elements fit under the packing assumptions.');
end

% Number of launches
N_launches_large = N_largeModules;

if N_smallElements==0
    N_launches_small = 0;
else
    N_launches_small = ...
        ceil(N_smallElements/N_smallPerFairing);
end

N_launches = N_launches_large + N_launches_small;

% Assembly operations
N_assemblyOperations = N_largeModules + N_smallElements;

% Store output
pack.N_geometric = N_geometric;
pack.N_smallPerFairing = N_smallPerFairing;
pack.N_launches_large = N_launches_large;
pack.N_launches_small = N_launches_small;
pack.N_launches = N_launches;
pack.N_assemblyOperations = N_assemblyOperations;
pack.PackagedElementVolume_m3 = prod(smallDims);

end
