
%% MODULAR SPACE STATION CONFIGURATION

structureName = 'Modular Space Station';

%% Structure geometry
envelope = [40;40;20];  % [m]
structureVolume = prod(envelope);

%% Assembly elements
N_largeModules = 8;
N_smallElements = 20;

%% Fairing packaging
packingBox = [3;3;8];
packingFactor = 0.50;
smallElementSize = [4;0.5;0.5];

pack = fairing_packaging(N_largeModules,N_smallElements, ...
    smallElementSize,packingBox,packingFactor);

N_smallPerFairing = pack.N_smallPerFairing;

%% Installation locations
hostPoint = [0;-30;0];

% Large module attachment locations
moduleLarge = [
    -15 -10  0
     -5 -10  0
      5 -10  0
     15 -10  0
    -15  10  0
     -5  10  0
      5  10  0
     15  10  0
];

% Small element attachment locations
[xGrid,yGrid] = meshgrid( ...
    linspace(-20,20,5), ...
    linspace(-20,20,4));

moduleSmall = [xGrid(:),yGrid(:), ...
    zeros(numel(xGrid),1)];

% Large elements first, small elements second
modulePoints = [moduleLarge;moduleSmall];

timesMin = [5 10 15 20];

%% Delivery containers - one per launch
N_smallContainers = ceil(N_smallElements/N_smallPerFairing);
N_containers = N_largeModules+N_smallContainers;
N_launches = N_containers;
N_deliveries = N_containers;

N_assemblyOperations = N_largeModules+N_smallElements;

elementsPerContainer = [ ...
    ones(N_largeModules,1); ...
    zeros(N_smallContainers,1)];

containerModuleIndices = cell(N_containers,1);

% Each large module occupies its own container
for k = 1:N_largeModules
    containerModuleIndices{k} = k;
end

% Pack small elements into the remaining containers
for k = 1:N_smallContainers
    first = (k-1)*N_smallPerFairing+1;
    last = min(k*N_smallPerFairing,N_smallElements);

    idx = first:last;
    c = N_largeModules+k;

    elementsPerContainer(c) = numel(idx);
    containerModuleIndices{c} = N_largeModules+idx;
end

assert(size(modulePoints,1)==N_assemblyOperations);
assert(sum(elementsPerContainer)==N_assemblyOperations);
