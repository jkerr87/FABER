function result = assembly_relocations_v2(mu, Re, h, Npts, ...
    timesMin, hostPoint, modulePoints)
% ASSEMBLY RELOCATIONS
%
% Repeated host -> module -> host transfers.
% No keep-out zone or collision constraints.
%
% INPUTS
% mu           [km^3/s^2]
% Re           [km]
% h            [km]
% Npts         Samples per HCW trajectory
% timesMin     Candidate times per transfer leg [min]
% hostPoint    Host location [3x1] [m]
% modulePoints Module locations [Nmodules x 3] [m]
%
% Requires existing HCW_solver.m

%% Constants

n = sqrt(mu/(Re+h)^3);

hostPoint = hostPoint(:);
Nmodules = size(modulePoints,1);

vHold = zeros(3,1);

%% Initialize results

DV_perModule = NaN(Nmodules,1);
DV_outbound = NaN(Nmodules,1);
DV_return = NaN(Nmodules,1);

TOF_perModule = NaN(Nmodules,1);
timePerLeg = NaN(Nmodules,1);

pathsOut = cell(Nmodules,1);
pathsReturn = cell(Nmodules,1);

%% Loop through module locations

for i = 1:Nmodules

    modulePos = modulePoints(i,:)';

    bestDV = inf;

    %% Optimize transfer time

    for timeMin = timesMin

        T = timeMin*60;

        try
            % Host -> module
            solOut = HCW_solver( ...
                hostPoint, vHold, ...
                modulePos, vHold, ...
                T, n, Npts);

            % Module -> host
            solBack = HCW_solver( ...
                modulePos, vHold, ...
                hostPoint, vHold, ...
                T, n, Npts);

            % Total round-trip delta V
            dvOut = solOut.dvTotal;
            dvBack = solBack.dvTotal;

            dvRoundTrip = dvOut + dvBack;

            if ~isfinite(dvRoundTrip)
                continue
            end

            if dvRoundTrip < bestDV

                bestDV = dvRoundTrip;

                DV_outbound(i) = dvOut;
                DV_return(i) = dvBack;

                timePerLeg(i) = timeMin;
                TOF_perModule(i) = 2*timeMin;

                pathsOut{i} = solOut.trajectory(1:3,:);
                pathsReturn{i} = solBack.trajectory(1:3,:);

            end

        catch ME
            warning('Module %d, T=%.1f min: %s', ...
                i,timeMin,ME.message);
        end

    end

    if isfinite(bestDV)
        DV_perModule(i) = bestDV;
    else
        warning('No HCW solution found for module %d.',i);
    end

end

%% Total mission results

result.feasible = all(isfinite(DV_perModule));

result.Nmodules = Nmodules;
result.DV_perModule = DV_perModule;
result.DV_outbound = DV_outbound;
result.DV_return = DV_return;

result.TOF_perModule_min = TOF_perModule;
result.timePerLeg_min = timePerLeg;

result.pathsOut = pathsOut;
result.pathsReturn = pathsReturn;

if result.feasible
    result.DV = sum(DV_perModule);
    result.TOF_min = sum(TOF_perModule);
else
    result.DV = NaN;
    result.TOF_min = NaN;
end

%% Print results

fprintf('\n=============================================\n');
fprintf('ASV ASSEMBLY RELOCATIONS\n');
fprintf('=============================================\n');

fprintf('Orbit altitude: %.1f km\n',h);
fprintf('Number of modules: %d\n',Nmodules);

ResultsTable = table( ...
    (1:Nmodules)', ...
    DV_outbound, ...
    DV_return, ...
    DV_perModule, ...
    TOF_perModule, ...
    timePerLeg, ...
    'VariableNames', { ...
    'Module', ...
    'DV_Outbound_mps', ...
    'DV_Return_mps', ...
    'DV_Total_mps', ...
    'TOF_min', ...
    'TimePerLeg_min'});

disp(ResultsTable);

fprintf('\nTotal Delta V: %.4f m/s\n',result.DV);
fprintf('Total Transfer Time: %.2f min\n',result.TOF_min);

end


function [dvTotal,path,minClearance,feasible] = ...
    evaluateRoute(waypoints,Tseg,n,Npts,halfForbidden,vHold)

dvTotal = 0;
path = [];
minClearance = inf;
feasible = true;

for k = 1:size(waypoints,2)-1

    try

        sol = HCW_solver( ...
            waypoints(:,k),vHold, ...
            waypoints(:,k+1),vHold, ...
            Tseg,n,Npts);

        r = sol.trajectory(1:3,:);

        signed = signedBoxDistance(r,halfForbidden);

        minClearance = min(minClearance,min(signed));

        if any(signed<=1e-8)
            feasible = false;
            dvTotal = NaN;
            return
        end

        dvTotal = dvTotal+sol.dvTotal;

        path = [path r];

    catch

        feasible = false;
        dvTotal = NaN;
        return

    end
end

end

%% ================================================================
%  HELPER: SIGNED DISTANCE FROM BOX
% ================================================================

function signed = signedBoxDistance(r,halfSize)

q = abs(r)-halfSize(:);

signed = vecnorm(max(q,0),2,1) + ...
    min(max(q,[],1),0);

end