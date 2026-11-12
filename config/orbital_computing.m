
%% ORBITAL COMPUTING INFRASTRUCTURE CONFIGURATION

structureName = 'Orbital Computing Infrastructure';

%% Structure geometry
envelope = [30;20;10];
structureVolume = prod(envelope);

%% Assembly elements
N_largeModules = 4;
N_smallElements = 12;

%% Fairing packaging
packingBox = [3;3;8];
packingFactor = 0.50;

% Preliminary radiator/truss packaged-size proxy
smallElementSize = [4;0.5;0.5];

pack = fairing_packaging(N_largeModules,N_smallElements, ...
    smallElementSize,packingBox,packingFactor);

N_smallPerFairing = pack.N_smallPerFairing;

%% Installation locations
hostPoint = [0;-20;0];

[xGrid,yGrid] = meshgrid( ...
    linspace(-15,15,4), ...
    linspace(-10,10,4));

modulePoints = [xGrid(:),yGrid(:), ...
    zeros(numel(xGrid),1)];

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

% One large module per container
for k = 1:N_largeModules
    containerModuleIndices{k} = k;
end

% Small-element containers
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
