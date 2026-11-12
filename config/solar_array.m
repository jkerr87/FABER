
%% SPACE SOLAR ARRAY CONFIGURATION

structureName = 'Space Solar Array';

%% Structure geometry
envelope = [40;20;2];
structureVolume = prod(envelope);

%% Assembly elements
N_panels = 20;
N_trusses = 10;

N_largeModules = 0;
N_smallElements = N_panels+N_trusses;

%% Fairing packaging
packingBox = [3;3;8];
packingFactor = 0.50;

% Common packaged-size proxy for panels and trusses
smallElementSize = [4;0.5;0.5];

pack = fairing_packaging(N_largeModules,N_smallElements, ...
    smallElementSize,packingBox,packingFactor);

N_smallPerFairing = pack.N_smallPerFairing;

%% Installation locations
hostPoint = [0;-15;0];

Nx = 6;
Ny = 5;

[xGrid,yGrid] = meshgrid( ...
    linspace(-envelope(1)/2,envelope(1)/2,Nx), ...
    linspace(-envelope(2)/2,envelope(2)/2,Ny));

modulePoints = [xGrid(:),yGrid(:), ...
    zeros(numel(xGrid),1)];

timesMin = [5 10 15 20];

%% Delivery containers - one per launch
N_smallContainers = ceil(N_smallElements/N_smallPerFairing);
N_containers = N_largeModules+N_smallContainers;
N_launches = N_containers;
N_deliveries = N_containers;

N_assemblyOperations = N_largeModules+N_smallElements;

elementsPerContainer = zeros(N_containers,1);
containerModuleIndices = cell(N_containers,1);

for k = 1:N_smallContainers
    first = (k-1)*N_smallPerFairing+1;
    last = min(k*N_smallPerFairing,N_smallElements);

    idx = first:last;

    elementsPerContainer(k) = numel(idx);
    containerModuleIndices{k} = idx;
end

assert(size(modulePoints,1)==N_assemblyOperations);
assert(sum(elementsPerContainer)==N_assemblyOperations);
