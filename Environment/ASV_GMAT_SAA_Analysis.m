%% ASV GMAT SUN ASPECT ANGLE ANALYSIS
% GMAT input file MUST have 10 numeric columns in this order:
% 1 elapsed days, 2:4 SC EarthMJ2000Eq position km,
% 5:7 SC EarthMJ2000Eq velocity km/s,
% 8:10 Sun EarthMJ2000Eq position km.
% The Sun epoch must be synchronized to the spacecraft epoch in GMAT.
% Body axes assumed fixed to the LVLH axes constructed below.
clear; clc; close all;

%% USER SETTINGS
inputFile = 'ASV_SAA_Results\FABER_ISS_SunOrbit_2038.txt'; % put in current MATLAB folder, or use full path
outDir = 'ASV_SAA_Results';
eclipseFile = 'ASV_SAA_Results\FABER_ISS_Eclipse_2038.txt';

figWidth = 6.2; figHeight = 3.9; % inches, ~half widescreen PowerPoint slide
Re = 6378.137;          % km, spherical Earth approximation
Rsun = 695700;          % km, approximate solar photosphere radius

if ~isfile(inputFile)
    error('GMAT report not found: %s\nPut the file beside this script or set inputFile to its full path.',inputFile);
end
if ~exist(outDir,'dir'), mkdir(outDir); end

%% READ / VALIDATE DATA

% Detect delimiter from first line
fid = fopen(inputFile,'r');
firstLine = fgetl(fid);
fclose(fid);

if contains(firstLine,'S')
    fprintf('Detected GMAT S-delimited format.\n');
    M = readmatrix(inputFile, ...
        'FileType','text','Delimiter','S');
else
    M = readmatrix(inputFile,'FileType','text');
end

M = M(~all(isnan(M),2),:);

if size(M,2) ~= 10
    error('Expected 10 columns, found %d.',size(M,2));
end

if any(~isfinite(M(:)))
    error('Invalid numeric values in GMAT report.');
end

% Verify orbital altitude
Re = 6378.137;
r0 = M(1,2:4);
h0 = norm(r0) - Re;

fprintf('Imported %d rows and %d columns.\n', ...
    size(M,1),size(M,2));
fprintf('Initial orbit altitude: %.2f km\n',h0);

if any(~isfinite(M(:)))
    error('Missing/non-numeric entries found. Check Sun position, report headers, and column order.');
end
M = sortrows(M,1);
[Mtime,iu] = unique(M(:,1),'stable');
M = M(iu,:);
tDays = Mtime-Mtime(1);
tHours = tDays*24;
r = M(:,2:4);  % km
v = M(:,5:7);  % km/s
rSun = M(:,8:10); % km, geocentric Sun position
if size(M,1)<3
    error('At least three samples are required.');
end
if median(vecnorm(rSun,2,2)) < 1e7
    warning('Sun position magnitude is unexpectedly small; verify Sun coordinates and epoch.');
end

%% SPACECRAFT-TO-SUN UNIT VECTOR
sunFromSC = rSun-r;
sHat = unitrows(sunFromSC);

%% LOCAL FRAME: +Z NADIR, +Y NEGATIVE ORBIT NORMAL, +X ALONG-TRACK
zBody = unitrows(-r);
hOrb = cross(r,v,2);
yBody = unitrows(-hOrb);
xBody = unitrows(cross(yBody,zBody,2));

%% SOLAR ASPECT ANGLES
SAA_X = acosd(clamp(dot(xBody,sHat,2)));
SAA_Y = acosd(clamp(dot(yBody,sHat,2)));
SAA_Z = acosd(clamp(dot(zBody,sHat,2)));

% Beta angle: signed Sun elevation above the orbital plane.
% This sign may differ from the convention used by a particular GMAT report.
beta = asind(clamp(dot(unitrows(hOrb),unitrows(rSun),2)));

%% SOLAR OCCULTATION (SPHERICAL EARTH/SUN APPROXIMATION)
% Angle between Earth-center and Sun-center as seen from spacecraft
rNorm = vecnorm(r,2,2);
sunDist = vecnorm(sunFromSC,2,2);
earthDir = unitrows(-r);
sep = acosd(clamp(dot(earthDir,sHat,2)));
earthAngularRadius = asind(clamp(Re./rNorm));
sunAngularRadius = asind(clamp(Rsun./sunDist));
% Full umbra: Sun disk entirely covered by Earth disk
isUmbra = sep <= (earthAngularRadius-sunAngularRadius);
% Any Earth/Sun disk overlap, including full umbra
isAnyEclipse = sep < (earthAngularRadius+sunAngularRadius);
isPenumbraOnly = isAnyEclipse & ~isUmbra;
isFullSun = ~isAnyEclipse;

% Direct sunlight cosine incidence on each +axis (no penumbra attenuation)
incX = max(0,cosd(SAA_X));
incY = max(0,cosd(SAA_Y));
incZ = max(0,cosd(SAA_Z));
incX(~isFullSun) = 0;
incY(~isFullSun) = 0;
incZ(~isFullSun) = 0;

%% REPORT: angles and eclipse statistics
fprintf('\n========== ASV SOLAR GEOMETRY REPORT ==========\n');
fprintf('Samples: %d | Spanned duration: %.4f days / %.3f hours\n',size(M,1),max(tDays),max(tHours));
fprintf('Median sampling interval: %.1f seconds\n',median(diff(Mtime))*86400);
fprintf('Estimated initial orbit altitude: %.2f km\n',rNorm(1)-Re);
fprintf('\nGeometric SAA (all sampled states, including eclipse):\n');
fprintf('%-12s %9s %9s %9s %9s %9s\n','Axis','Min [deg]','Mean','Max','P5','P95');
labels={'LVLH +X','LVLH +Y','LVLH +Z'};
angles={SAA_X,SAA_Y,SAA_Z};
for j=1:3
  A=angles{j};
  fprintf('%-12s %9.2f %9.2f %9.2f %9.2f %9.2f\n',labels{j},min(A),mean(A),max(A),pct(A,5),pct(A,95));
end
fprintf('\nBeta angle [deg]: min %.2f | mean %.2f | max %.2f\n',min(beta),mean(beta),max(beta));




% Time-weighted eclipse fractions using interval midpoint classification.
% Good for approximate fractions if time steps resolve eclipse transitions;
% use GMAT EclipseLocator for accurate event durations and maxima.
dtSeconds=diff(Mtime)*86400;
midUmbra=isUmbra(1:end-1) & isUmbra(2:end);
midEclipse=isAnyEclipse(1:end-1) & isAnyEclipse(2:end);
Tseconds=sum(dtSeconds);
fprintf('\nApprox. fraction of span in full umbra: %.1f %%\n',100*sum(dtSeconds(midUmbra))/Tseconds);
fprintf('Approx. fraction of span in any eclipse: %.1f %%\n',100*sum(dtSeconds(midEclipse))/Tseconds);
fprintf('(Fractions depend on sample spacing. Use GMAT EclipseLocator for eclipse duration.)\n');

% Count discrete sampled intervals classified as eclipse events
starts=find(diff([false;isAnyEclipse])==1);
ends=find(diff([isAnyEclipse;false])==-1);
if numel(starts)~=numel(ends)
    warning('Unexpected eclipse edge alignment; inspect samples.');
else
    completed=starts>1 & ends<numel(isAnyEclipse);
    durationsMin=(Mtime(ends(completed))-Mtime(starts(completed)))*1440;
    fprintf('Completed sampled eclipse events: %d\n',nnz(completed));
    if ~isempty(durationsMin)
      fprintf('Approx. longest sampled eclipse: %.1f minutes\n',max(durationsMin));
    end
end
fprintf('\nMean direct-sun incidence factors (+X,+Y,+Z) with eclipse mask: %.3f, %.3f, %.3f\n',mean(incX),mean(incY),mean(incZ));
fprintf('WARNING: These incidence factors neglect penumbra partial flux, Earth albedo, self-shadowing and solar array articulation.\n');
fprintf('================================================\n\n');


%% EXPORT SAMPLED ECLIPSE INTERVALS
% Independent estimate from the GMAT 10-column state history.
% GMAT EclipseLocator report is authoritative for precise
% entry/exit times and event durations.

epochUTC = datetime(2038,1,1,0,0,0,'TimeZone','UTC');

% Detect eclipse transitions
eclipseStartIdx = find(diff([false; isAnyEclipse]) == 1);
eclipseEndIdx   = find(diff([isAnyEclipse; false]) == -1);

startUTC = epochUTC + days(Mtime(eclipseStartIdx));
endUTC   = epochUTC + days(Mtime(eclipseEndIdx));

durationMin = minutes(endUTC - startUTC);

EclipseEvents = table( ...
    startUTC, endUTC, durationMin, ...
    'VariableNames', ...
    {'StartUTC','EndUTC','Duration_min'});

writetable(EclipseEvents, ...
    fullfile(outDir,'FABER_SampledEclipseEvents.csv'));

fprintf('\nSampled eclipse events: %d\n',height(EclipseEvents));

if ~isempty(EclipseEvents)
    fprintf('Maximum sampled eclipse: %.2f min\n', ...
        max(EclipseEvents.Duration_min));
end



%% SAVE DATA / SUMMARY
Results=table(Mtime,tDays,tHours,SAA_X,SAA_Y,SAA_Z,beta,isUmbra,isPenumbraOnly,isAnyEclipse,incX,incY,incZ,...
    'VariableNames',{'GMAT_ElapsedDays','TimeFromFirstSample_days','TimeFromFirstSample_hours',...
    'SAA_X_deg','SAA_Y_deg','SAA_Z_deg','Beta_deg','Umbra','PenumbraOnly','AnyEclipse',...
    'SunIncidence_PosX','SunIncidence_PosY','SunIncidence_PosZ'});
writetable(Results,fullfile(outDir,'ASV_SAA_TimeHistory.csv'));
Summary=table(labels',cellfun(@min,angles)',cellfun(@mean,angles)',cellfun(@max,angles)',...
    cellfun(@(a)pct(a,5),angles)',cellfun(@(a)pct(a,95),angles)',...
    'VariableNames',{'Axis','Min_deg','Mean_deg','Max_deg','P5_deg','P95_deg'});
writetable(Summary,fullfile(outDir,'ASV_SAA_Summary.csv'));

%% FIGURE 1: SAA vs time
f=figure('Color','w','Units','inches','Position',[1 1 figWidth figHeight]);
plot(tHours,SAA_X,'LineWidth',1.3);hold on
plot(tHours,SAA_Y,'LineWidth',1.3);
plot(tHours,SAA_Z,'LineWidth',1.3);
grid on; xlabel('Elapsed time from first sample [h]');ylabel('SAA [deg]');
title('ASV Sun aspect angles | LVLH-pointing');
legend('+X along-track','+Y negative orbit normal','+Z nadir','Location','best');
ylim([0 180]);
exportgraphics(f,fullfile(outDir,'ASV_SAA_Axes.png'),'Resolution',300);
exportgraphics(f,fullfile(outDir,'ASV_SAA_Axes.svg'),'ContentType','vector');

%% FIGURE 2: Solar beta angle
f=figure('Color','w','Units','inches','Position',[1 1 figWidth figHeight]);
plot(tDays,beta,'LineWidth',1.4);grid on
xlabel('Elapsed time from first sample [days]');ylabel('Beta angle [deg]');
title('ASV solar beta angle');
exportgraphics(f,fullfile(outDir,'ASV_Beta_Angle.png'),'Resolution',300);
exportgraphics(f,fullfile(outDir,'ASV_Beta_Angle.svg'),'ContentType','vector');

%% FIGURE 3: Solar exposure/incidence on +X/+Y/+Z (idealized)
f=figure('Color','w','Units','inches','Position',[1 1 figWidth figHeight]);
plot(tHours,incX,'LineWidth',1.2);hold on
plot(tHours,incY,'LineWidth',1.2);
plot(tHours,incZ,'LineWidth',1.2);
grid on; xlabel('Elapsed time from first sample [h]');
ylabel('Direct-sun cosine incidence factor [-]');
title('ASV idealized illumination of +body faces');
legend('+X','+Y','+Z','Location','best');ylim([0 1.05]);
exportgraphics(f,fullfile(outDir,'ASV_Sun_Incidence.png'),'Resolution',300);
exportgraphics(f,fullfile(outDir,'ASV_Sun_Incidence.svg'),'ContentType','vector');

fprintf('Output CSV and PNG/SVG plots saved in: %s\n',fullfile(pwd,outDir));


% Approximate direct solar-array power for a fixed +X-facing panel
AU = 149597870.7;       % km
solarFlux1AU = 1361;    % W/m^2; refine for your mission/model
panelArea = 1.0;        % m^2 -- replace with actual illuminated area
efficiency = 0.30;      % replace with array efficiency

% +X incidence; zero when the Sun is behind the panel
cosInc = max(0, dot(xBody, sHat, 2));

% Simple full-sun-only model; penumbra is treated as zero sunlight here
sunFraction = double(isFullSun);

solarFlux = solarFlux1AU .* (AU ./ sunDist).^2;
Parray_W = solarFlux .* panelArea .* efficiency .* cosInc .* sunFraction;


%% LOCAL FUNCTIONS
function U=unitrows(A)
  den=vecnorm(A,2,2);
  if any(den<eps),error('Cannot normalize a zero-length vector. Check GMAT data.');end
  U=A./den;
end
function x=clamp(x)
  x=max(-1,min(1,x));
end
function y=pct(x,p)
  x=sort(x(isfinite(x)));n=numel(x);
  if n==0,y=NaN;return;end
  q=1+(n-1)*p/100;lo=floor(q);hi=ceil(q);
  y=x(lo)+(q-lo)*(x(hi)-x(lo));
end
