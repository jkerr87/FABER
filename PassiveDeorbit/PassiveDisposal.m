function [disposalResults, burnResults, fig] = PassiveDisposal(Cd, mass_kg, areaTest_m2)
%PASSIVEDISPOSAL Estimate drag-only and burn-assisted orbital disposal.
% Inputs:
%   Cd          - Positive drag coefficient [-].
%   mass_kg     - Positive disposal mass, including residuals [kg].
%   areaTest_m2 - Positive scalar or vector of projected drag areas [m^2].
%                 Each area is an alternative for the same mass and Cd.
% Outputs:
%   disposalResults - Minimum area for five-year decay and natural decay
%                     time for each orbit, area and density profile.
%   burnResults     - Minimum insertion DeltaV and post-burn decay time for
%                     each orbit, area and profile. DeltaV is zero when
%                     natural decay meets the five-year limit.
%   fig             - Two-panel figure of insertion DeltaV versus area.
%
% The model uses fixed mass, Cd and projected area, fixed density profiles
% and orbital speed for drag. The initial orbits are circular. The 160 km
% cutoff is applied to the circular altitude or to the orbit-averaged
% perigee. The burn budget covers the insertion burn only.

%% 1. Check inputs and set up the reference cases
validateattributes(Cd, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'}, mfilename, 'Cd', 1);
validateattributes(mass_kg, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'}, mfilename, 'mass_kg', 2);
validateattributes(areaTest_m2, {'numeric'}, ...
    {'nonempty', 'vector', 'real', 'finite', 'positive'}, mfilename, 'areaTest_m2', 3);

% Use full double arrays for the calculations and sort the areas for plotting.
Cd = full(double(Cd));
mass_kg = full(double(mass_kg));
areaTest_m2 = sort(full(double(areaTest_m2(:))));

% PM2 reference orbits.
orbit = ["ISS-like low"; "ISS-like nominal"; "ISS-like high"; ...
         "SSO low"; "SSO nominal"; "SSO high"];
hStart = [389.78; 457.25; 524.72; 607.28; 680.12; 752.96];  % [km]
profile = ["Moderate"; "Low"];

mu = 3.986004418e14;       % Earth's gravitational parameter [m^3/s^2]
Re = 6378137;              % Earth radius [m]
hStop = 160;               % Numerical decay cutoff [km]
limitYears = 5;            % Provisional time limit from release [years]
secPerYear = 365.25 * 86400;
limitSec = limitYears * secPerYear;

% ECSS-E-ST-10-04C Rev.1, Table G-4 (JB-2006).
% Density [kg/m^3]: moderate activity in column 1, low activity in column 2.
hTable = (160:20:760).';   % [km]
rhoTable = [
    1.36E-09  9.17E-10  % 160 km
    6.15E-10  3.45E-10  % 180 km
    3.17E-10  1.47E-10  % 200 km
    1.77E-10  6.96E-11  % 220 km
    1.05E-10  3.54E-11  % 240 km
    6.47E-11  1.88E-11  % 260 km
    4.12E-11  1.03E-11  % 280 km
    2.69E-11  5.86E-12  % 300 km
    1.80E-11  3.40E-12  % 320 km
    1.23E-11  2.02E-12  % 340 km
    8.48E-12  1.22E-12  % 360 km
    5.95E-12  7.46E-13  % 380 km
    4.22E-12  4.63E-13  % 400 km
    3.02E-12  2.92E-13  % 420 km
    2.18E-12  1.87E-13  % 440 km
    1.59E-12  1.21E-13  % 460 km
    1.17E-12  8.04E-14  % 480 km
    8.60E-13  5.44E-14  % 500 km
    6.39E-13  3.77E-14  % 520 km
    4.77E-13  2.68E-14  % 540 km
    3.58E-13  1.96E-14  % 560 km
    2.71E-13  1.47E-14  % 580 km
    2.06E-13  1.14E-14  % 600 km
    1.57E-13  9.10E-15  % 620 km
    1.20E-13  7.41E-15  % 640 km
    9.28E-14  6.16E-15  % 660 km
    7.19E-14  5.22E-15  % 680 km
    5.60E-14  4.50E-15  % 700 km
    4.40E-14  3.93E-15  % 720 km
    3.48E-14  3.48E-15  % 740 km
    2.79E-14  3.10E-15  % 760 km
];
logRho = log(rhoTable);
nOrbit = numel(orbit);
nArea = numel(areaTest_m2);
nProfile = numel(profile);

%% 2. Passive disposal from a circular orbit
minArea = zeros(nOrbit, nProfile);
naturalYears = zeros(nOrbit, nArea, nProfile);

for j = 1:nProfile
    % Interpolate log-density between the tabulated altitudes.
    rho = @(h) exp(interp1(hTable, logRho(:, j), h, 'linear'));

    % Integrate the inverse of the circular altitude decay rate.
    % The factor 1000 converts the integration variable from km to m.
    % The constant Cd*A/m is applied separately after integration.
    integrand = @(h) 1000 ./ (rho(h) .* sqrt(mu * (Re + 1000 * h)));

    for i = 1:nOrbit
        h0 = hStart(i);
        knots = hTable(hTable > hStop & hTable < h0);
        decayIntegral = integral(integrand, hStop, h0, ...
            'Waypoints', knots, 'RelTol', 1e-9);
        minArea(i, j) = mass_kg * decayIntegral / (Cd * limitSec);
    end

    % Natural decay time scales as 1/A and can exceed the time limit.
    naturalYears(:, :, j) = limitYears * minArea(:, j) ./ areaTest_m2.';
end

[rowOrbit, rowArea] = ndgrid(1:nOrbit, 1:nArea);
rowOrbit = rowOrbit(:);
rowArea = rowArea(:);
moderateYears = naturalYears(:, :, 1);
lowYears = naturalYears(:, :, 2);
disposalResults = table(orbit(rowOrbit), hStart(rowOrbit), ...
    minArea(rowOrbit, 1), minArea(rowOrbit, 2), areaTest_m2(rowArea), ...
    moderateYears(:), lowYears(:), ...
    'VariableNames', {'Orbit', 'Altitude_km', 'MinArea_moderate_m2', ...
    'MinArea_low_m2', 'Area_m2', 'NaturalDecay_moderate_years', ...
    'NaturalDecay_low_years'});
disp(disposalResults)

%% 3. Disposal after one tangential retrograde burn
tolKm = 0.1;  % Perigee search tolerance [km]
rStop = Re + 1000 * hStop;
E = ((0:511) + 0.5) * pi / 512;  % Eccentric anomaly on half an orbit
cosE = cos(E);
options = odeset('RelTol', 1e-8, 'AbsTol', [0.01, 1e-10], ...
    'Events', @(t, x) stopAtPerigee(t, x, rStop));

% Start with the no-burn results. Update only cases that need a burn.
perigee = repmat(hStart, 1, nArea, nProfile);
dv = zeros(nOrbit, nArea, nProfile);
decayYears = naturalYears;

for j = 1:nProfile
    % Clamp density queries from trial steps below the terminal altitude.
    rho = @(r) exp(interp1(hTable, logRho(:, j), ...
        max(hStop, (r - Re) / 1000), 'linear'));

    for i = 1:nOrbit
        h0 = hStart(i);
        rA = Re + 1000 * h0;   % The burn altitude becomes the apogee.
        vCirc = sqrt(mu / rA);

        for k = 1:nArea
            area = areaTest_m2(k);
            if area >= minArea(i, j)
                continue
            end

            drag = Cd * area / mass_kg;
            rates = @(t, x) orbitDecay(t, x, mu, drag, rho, cosE);
            [hP, decaySec] = findPerigee(hStop, h0, tolKm, Re, rA, ...
                limitSec, rates, options);

            % Use vis-viva to calculate the speed change at apogee.
            rP = Re + 1000 * hP;
            a = (rA + rP) / 2;
            vA = sqrt(mu * (2 / rA - 1 / a));
            perigee(i, k, j) = hP;
            dv(i, k, j) = vCirc - vA;
            decayYears(i, k, j) = decaySec / secPerYear;
        end
    end
end

[rowOrbit, rowArea, rowProfile] = ndgrid(1:nOrbit, 1:nArea, 1:nProfile);
burnResults = table(orbit(rowOrbit(:)), hStart(rowOrbit(:)), ...
    profile(rowProfile(:)), areaTest_m2(rowArea(:)), perigee(:), ...
    dv(:), decayYears(:), 'VariableNames', ...
    {'Orbit', 'Altitude_km', 'Activity', 'Area_m2', 'Perigee_km', ...
     'DeltaV_m_s', 'Decay_years'});
disp(burnResults)

%% 4. Plot the insertion DeltaV for both density profiles
fig = figure('Name', 'Disposal analysis');
tiledlayout(1, 2);
for j = 1:nProfile
    nexttile
    for i = 1:nOrbit
        semilogx(areaTest_m2, dv(i, :, j).', '-o');
        hold on
    end
    grid on
    xlabel('Effective drag area [m^2]')
    ylabel('Minimum insertion \DeltaV [m/s]')
    title(profile(j) + " activity")
end
legend(orbit, 'Location', 'best')
end

function [hP, decaySec] = findPerigee(hStop, hStart, tolKm, Re, rA, ...
    limitSec, rates, options)
% Find the highest initial perigee that meets the time limit.
hLow = hStop;    % Feasible: the initial perigee is already at the cutoff.
hHigh = hStart;  % Infeasible: natural decay exceeds the time limit.
decaySec = 0;    % Elapsed time for an initial perigee at the cutoff.

while hHigh - hLow > tolKm
    hTrial = (hLow + hHigh) / 2;
    rP = Re + 1000 * hTrial;
    a = (rA + rP) / 2;
    e = (rA - rP) / (rA + rP);

    [t, ~, tEvent] = ode45(rates, [0, limitSec], [a; e], options);
    if ~isempty(tEvent)
        hLow = hTrial;
        decaySec = tEvent(1);
    elseif t(end) < limitSec
        error('PassiveDisposal:IntegrationFailed', ...
            'Orbit decay integration stopped early for a trial perigee of %.3f km.', hTrial);
    else
        hHigh = hTrial;
    end
end

hP = hLow;
end

function dx = orbitDecay(~, x, mu, drag, rho, cosE)
% Average the energy and angular momentum losses over one Keplerian orbit.
a = x(1);  % Semimajor axis [m]
e = x(2);  % Eccentricity [-]
weight = 1 - e * cosE;  % Time averaging weight for eccentric anomaly.
r = a * weight;
v = sqrt(mu * (2 ./ r - 1 / a));
density = rho(r);

da = -drag * a^2 / mu * mean(density .* v.^3 .* weight);
de = -drag * (1 - e^2) * mean(density .* v .* cosE);
dx = [da; de];
end

function [value, isTerminal, direction] = stopAtPerigee(~, x, rStop)
% Stop when the perigee radius falls to the cutoff.
value = x(1) * (1 - x(2)) - rStop;
isTerminal = 1;
direction = -1;
end
