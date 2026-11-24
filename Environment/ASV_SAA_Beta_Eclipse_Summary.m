%% ASV: SAA of six LVLH faces, solar beta, and eclipse duration
% Input: GMAT ReportFile, 10 columns: elapsed days, rSC(3), vSC(3), rSun(3)
% Units from GMAT: days, km, km/s, km.
% Outputs: command-window summary, CSV tables and two figures.
clear; clc; close all;

inputFile = 'FABER_ISS_Eclipse_2038.txt';
muEarth = 398600.4418;     % km^3/s^2
RE = 6378.137;            % km, spherical Earth approximation

%% 1) Import: GMAT can repeat headers during repeated runs
opts = delimitedTextImportOptions('NumVariables',10);
opts.DataLines = [2 Inf];
opts.Delimiter = {' ', '\t'};
opts.ConsecutiveDelimitersRule = 'join';
opts.LeadingDelimitersRule = 'ignore';
opts.VariableTypes = repmat({'double'},1,10);
opts.ExtraColumnsRule = 'ignore';
opts.EmptyLineRule = 'skip';
D = readmatrix(inputFile,opts);
D = D(all(isfinite(D),2),:);            % remove repeated header rows
assert(size(D,2)==10 && size(D,1)>2, 'Expected 10 numeric GMAT columns.');
D = sortrows(D,1);
[~,ia] = unique(round(D(:,1)*86400,6),'stable'); % remove repeated timestamps
D = D(sort(ia),:);

t = D(:,1);                            % elapsed days
r = D(:,2:4);                          % Earth to spacecraft, km
v = D(:,5:7);                          % spacecraft velocity, km/s
rSun = D(:,8:10);                     % Earth to Sun, km
normRows = @(x) x ./ vecnorm(x,2,2);

%% 2) LVLH basis: +X along-track, +Y negative orbit normal, +Z nadir
z = normRows(-r);
h = cross(r,v,2);
y = normRows(-h);
x = normRows(cross(y,z,2));
s = normRows(rSun-r);                   % unit spacecraft-to-Sun vector

axesMatrix = {x,-x,y,-y,z,-z};
faceNames = {'+X';'-X';'+Y';'-Y';'+Z';'-Z'};
SAA = zeros(numel(t),6);
for j = 1:6
    cosang = sum(axesMatrix{j}.*s,2);
    SAA(:,j) = acosd(max(-1,min(1,cosang)));
end

%% 3) Beta angle (signed, relative to orbit normal)
hHat = normRows(h);
beta = asind(max(-1,min(1,sum(hHat.*s,2))));

%% 4) Idealized umbral eclipse duration for circular orbit at each epoch
% Geometric cylindrical-shadow, spherical Earth, point-source Sun model.
% For a circular orbit, this is the ideal eclipse duration per orbit,
% NOT a direct measurement of individual GMAT EclipseLocator events.
rmag = vecnorm(r,2,2);
v2 = sum(v.^2,2);
a = 1 ./ (2./rmag - v2/muEarth);       % osculating semimajor axis, km
n = sqrt(muEarth./a.^3);             % mean motion, rad/s
Torb_min = 2*pi./n/60;
betaCritical = asind(min(1,RE./rmag));
arg = sqrt(max(0,1-(RE./rmag).^2)) ./ max(cosd(beta),eps);
alpha = acos(max(-1,min(1,arg)));     % half eclipse angle, radians
Eclipse_min = 2*alpha./n/60;
Eclipse_min(abs(beta)>=betaCritical) = 0;

%% 5) Incidence fraction and eclipse-aware ideal illumination
% Approximate Earth occultation from geometry at each reported sample.
% Note: when sampling hourly, short eclipse passages can be missed.
sunTowardEarth = sum(r.*s,2)<0;
perpDist = vecnorm(cross(r,s,2),2,2);
inUmbraSample = sunTowardEarth & (perpDist < RE);
incidence = max(0,cosd(SAA));
illumination = incidence .* double(~inUmbraSample);

%% 6) Print numerical results
fprintf('\nASV SOLAR GEOMETRY SUMMARY\n');
fprintf('Analyzed interval: %.2f days (%.2f years)\n',t(end)-t(1),(t(end)-t(1))/365.25);
fprintf('Unique GMAT samples: %d\n',numel(t));
fprintf('Solar beta angle: min %.2f deg, max %.2f deg\n',min(beta),max(beta));
fprintf('Orbital period: %.2f to %.2f min\n',min(Torb_min),max(Torb_min));
fprintf('Maximum ideal eclipse duration: %.2f min/orbit\n',max(Eclipse_min));
fprintf('Minimum ideal eclipse duration: %.2f min/orbit\n',min(Eclipse_min));
fprintf('Maximum ideal eclipse fraction: %.1f %%\n',100*max(Eclipse_min./Torb_min));

faceSummary = table(faceNames,min(SAA,[],1)',max(SAA,[],1)',mean(SAA,1)', ...
    min(incidence,[],1)',max(incidence,[],1)', ...
    'VariableNames',{'Face','SAA_min_deg','SAA_max_deg','SAA_mean_deg', ...
                      'Incidence_min','Incidence_max'});
disp(faceSummary)
writetable(faceSummary,'ASV_Face_SAA_Summary.csv');

result = table(t,beta,Eclipse_min,Torb_min, ...
    'VariableNames',{'ElapsedDays','Beta_deg','IdealEclipse_min','OrbitPeriod_min'});
for j=1:6
    nm = strrep(faceNames{j},'+','plus');
    nm = strrep(nm,'-','minus');
    result.(['SAA_' nm '_deg']) = SAA(:,j);
    result.(['SunIncidence_' nm]) = illumination(:,j);
end
writetable(result,'ASV_Solar_Geometry_Results.csv');

%% 7) Plots: long-term beta and theoretical eclipse, short-term SAA
f1=figure('Color','w','Units','inches','Position',[1 1 6.2 4]);
plot(t,beta,'LineWidth',1); grid on
xlabel('Elapsed days'); ylabel('Solar beta angle [deg]');
title('ASV solar beta-angle evolution');
exportgraphics(f1,'ASV_Beta_Angle.png','Resolution',300);

f2=figure('Color','w','Units','inches','Position',[1 1 6.2 4]);
plot(t,Eclipse_min,'LineWidth',1); grid on
xlabel('Elapsed days'); ylabel('Ideal eclipse duration [min/orbit]');
title('ASV theoretical eclipse duration');
exportgraphics(f2,'ASV_Eclipse_Duration.png','Resolution',300);

short = t <= t(1) + 3*Torb_min(1)/1440; % about three orbits
f3=figure('Color','w','Units','inches','Position',[1 1 6.2 4]);
plot((t(short)-t(1))*24,SAA(short,:),'LineWidth',1.1); grid on
xlabel('Elapsed time [h]'); ylabel('Sun aspect angle [deg]');
title('ASV SAA on all six LVLH faces'); ylim([0 180]);
legend(faceNames,'Location','eastoutside');
exportgraphics(f3,'ASV_SAA_Six_Faces.png','Resolution',300);
