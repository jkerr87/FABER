%% ISS-BAND LARGE STRUCTURES TYPE 1 CONFIGURATION (ID 4)
% Scenario inputs transcribed from the nine-use-case table.
% Geometry, packaging dimensions and hostPoint are ASSUMPTIONS, not CAD.
% Run from a folder containing fairing_packaging.m and
% FABER_make_config_geometry.m. This script defines workspace variables.

%% Identification
configID = 4;
structureName = 'ISS-band Large Structures Type 1';
vehicles = 'Veicoletto';

%% Structure geometry
% Envelope refers to a structure bounding box, not fairing dimensions.
envelope = [40; 40; 20]; % [m]
structureVolume = prod(envelope); % bounding-box volume [m^3]

%% Assembly elements (table inputs)
N_largeModules = 7;
N_smallElements = 970; % labelled 'raw material' in source table
% WARNING: N_smallElements is treated as separately installed objects
% ONLY as a mission-analysis proxy; confirm actual manufacturing workflow.

%% Fairing packaging (illustrative, inherited from baseline)
packingBox = [3; 3; 8]; % [m]
packingFactor = 0.50;
smallElementSize = [4; 0.5; 0.5]; % [m], placeholder for ALL cases
pack = fairing_packaging(N_largeModules,N_smallElements, ...
    smallElementSize,packingBox,packingFactor);
N_smallPerFairing = pack.N_smallPerFairing;

%% Installation locations
% Host point is a placeholder; geometrical clearance must be verified.
hostPoint = [0; -0.75*envelope(2); 0]; % [m]
[moduleLarge,moduleSmall] = FABER_make_config_geometry( ...
    envelope,N_largeModules,N_smallElements);
modulePoints = [moduleLarge;moduleSmall]; % [N x 3] row-wise [x,y,z]
timesMin = [5 10 15 20]; % candidate travel times [min]

%% Deliveries / containers
% IMPORTANT: the original table calls N_containers ONLY the small-piece
% container count; large modules are additional separate delivery units.
N_containers = 27; % exact table entry, not overwritten by packing
N_smallContainers = N_containers;
N_deliveryUnits = N_largeModules + N_smallContainers;
N_deliveries = N_deliveryUnits; % one delivery per module / small container
% Estimate MINIMUM launches based on the illustrative packaging model.
% Not the same thing as user's declared N_containers.
N_launches = pack.N_launches;
N_requiredSmallContainers = pack.N_launches_small;
packingConsistent = N_smallContainers >= N_requiredSmallContainers;
if ~packingConsistent
    warning('FABER:PackagingMismatch', ...
     ['ID %d declares %d small containers but fairing_packaging ' ...
      'requires at least %d under assumed dimensions.'], ...
      configID,N_smallContainers,N_requiredSmallContainers);
end

%% Installation operations vs. reported vehicle relocations
N_assemblyOperations = N_largeModules + N_smallElements;
% These TWO values are preserved from the table. Their meanings are
% not established well enough to multiply into the trajectory budget.
N_relocationsPerDeliveryLarge = 14;
N_relocationsPerDeliveryContainer = 78;
% Do not set N_assemblyOperations equal to reported relocation events.
% In particular, Taxi + Veicoletto participation does not always imply
% two distinct assembly-transfer trajectories.

%% Assignment of elements to DELIVERY UNITS
% Legacy variable name elementsPerContainer is maintained for existing
% downstream code, but first entries are individual large modules.
elementsPerDelivery = [ones(N_largeModules,1); ...
                       zeros(N_smallContainers,1)];
containerModuleIndices = cell(N_deliveryUnits,1);
for k=1:N_largeModules
    containerModuleIndices{k} = k;
end
for k=1:N_smallContainers
    first = floor((k-1)*N_smallElements/N_smallContainers)+1;
    last  = floor(k*N_smallElements/N_smallContainers);
    idx = first:last;
    c = N_largeModules+k;
    elementsPerDelivery(c) = numel(idx);
    containerModuleIndices{c} = N_largeModules+idx;
end
elementsPerContainer = elementsPerDelivery; % compatibility alias

%% Integrity checks
assert(size(moduleLarge,1)==N_largeModules);
assert(size(moduleSmall,1)==N_smallElements);
assert(size(modulePoints,1)==N_assemblyOperations);
assert(sum(elementsPerDelivery)==N_assemblyOperations);
assert(numel(containerModuleIndices)==N_deliveryUnits);

fprintf('ID %d | %s | modules %d | small %d | deliveries %d | declared small containers %d | packing min %d | estimated launches %d\n', ...
 configID,structureName,N_largeModules,N_smallElements,N_deliveries, ...
 N_smallContainers,N_requiredSmallContainers,N_launches);
