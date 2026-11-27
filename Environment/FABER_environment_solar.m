function model = FABER_environment_solar(epoch,maxDays,useToolbox,gmatCluster)
% Earth-centered Sun ephemeris and sampled Earth occultation geometry.
model.gmatEclipseIntervals = table();
if nargin>=4 && strlength(gmatCluster)>0
    model = gmatSolar(epoch,maxDays,gmatCluster);
    return;
end
grid = unique([0:300:maxDays*86400,maxDays*86400])';
jd = juliandate(epoch+seconds(grid));
model.source = "Analytic Sun (approximate)";
if useToolbox
    try
        sun = 1000*planetEphemeris(jd,'Earth','Sun', ...
            StartDate=floor(jd(1))-1,EndDate=ceil(jd(end))+1);
        model.source = "Aerospace Toolbox planetary ephemeris";
    catch err
        model.ephemerisDiagnostic = string(err.message);
    end
end
if ~exist('sun','var')
    d = jd-2451545;
    meanAnomaly = 357.529+0.98560028*d;
    longitude = 280.459+0.98564736*d+1.915*sind(meanAnomaly)+0.020*sind(2*meanAnomaly);
    obliquity = 23.439-0.00000036*d;
    au = 149597870700*(1.00014-0.01671*cosd(meanAnomaly)-0.00014*cosd(2*meanAnomaly));
    sun = au.*[cosd(longitude),sind(longitude).*cosd(obliquity), ...
        sind(longitude).*sind(obliquity)];
end
model.sun = @(t) interp1(grid,sun,t(:),'linear')';
end

function model = gmatSolar(epoch,maxDays,cluster)
cluster = validatestring(cluster,{'ISS','SSO'});
folder = fullfile(fileparts(mfilename('fullpath')),'ASV_SAA_Results');
sunFile = fullfile(folder,"FABER_"+cluster+"_SunOrbit_2038.txt");
eclipseFile = fullfile(folder,"FABER_"+cluster+"_Eclipse_2038.txt");
reportEpoch = datetime(2038,1,1,'TimeZone','UTC');
endUTC = epoch+days(maxDays);
data = readmatrix(sunFile,'FileType','text','Delimiter','S');
assert(size(data,2)==10,'Expected 10 columns in %s.',sunFile);
valid = all(isfinite(data(:,[1 8:10])),2);
times = 86400*data(valid,1)-seconds(epoch-reportEpoch);
[times,order] = unique(times);
sun = 1000*data(valid,8:10);
assert(numel(times)>1 && times(1)<=0 && times(end)>=maxDays*86400, ...
    'Requested epoch and horizon must fall inside the GMAT SunOrbit export.');
first = find(times<=0,1,'last');
last = find(times>=maxDays*86400,1);
grid = times(first:last);
sun = sun(order(first:last),:);
model.sun = @(t) interp1(grid,sun,t(:),'linear')';
model.source = "GMAT "+cluster+" SunOrbit (EarthMJ2000Eq)";

% Locator events describe the GMAT reference spacecraft. They are not used
% to gate maneuvers on a different simulated spacecraft trajectory.
date = '(\d{2} [A-Za-z]{3} \d{4} \d{2}:\d{2}:\d{2}\.\d+)';
pattern = ['(?m)^\s*' date '\s+' date ...
    '\s+([\d.]+)\s+Earth\s+(Umbra|Penumbra|Antumbra)\s+(\d+)'];
rows = regexp(fileread(eclipseFile),pattern,'tokens');
assert(~isempty(rows),'No Earth eclipse events found in %s.',eclipseFile);
rows = vertcat(rows{:});
startUTC = datetime(rows(:,1),'InputFormat','dd MMM yyyy HH:mm:ss.SSS', ...
    'Locale','en_US','TimeZone','UTC');
stopUTC = datetime(rows(:,2),'InputFormat','dd MMM yyyy HH:mm:ss.SSS', ...
    'Locale','en_US','TimeZone','UTC');
inWindow = startUTC<endUTC & stopUTC>epoch;
model.gmatEclipseIntervals = table(startUTC(inWindow),stopUTC(inWindow), ...
    str2double(rows(inWindow,3)),string(rows(inWindow,4)),str2double(rows(inWindow,5)), ...
    'VariableNames',{'StartUTC','EndUTC','ReportedDuration_s','Type','EventNumber'});
end
