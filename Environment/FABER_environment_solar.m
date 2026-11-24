function model = FABER_environment_solar(epoch,maxDays,useToolbox)
% Earth-centered Sun ephemeris and sampled Earth occultation geometry.
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
