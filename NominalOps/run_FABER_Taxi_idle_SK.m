function SK = run_FABER_Taxi_idle_SK(follower,reference,opts)
% Periodic two-impulse HCW relative stationkeeping against differential drag.
% follower and reference each supply mass_kg, area_m2, Cd and optionally
% staticDistance_m, ellipseCenter_m, ellipseRadialAmplitude_m and
% ellipsePhase_rad, measured in the same circular chief LVLH frame.
% Positive distance is behind the chief. The follower executes corrections.
% opts requires altitude_m and rho_kgm3. Optional: duration_days,
% correctionInterval_hr, IspRCS_s, AZ_m, AZmargin_m.
% HCW assumes a circular chief; the supplied ellipses are relative orbits.
arguments
    follower (1,1) struct
    reference (1,1) struct
    opts (1,1) struct
end
follower = vehicleInputs(follower);
reference = vehicleInputs(reference);
opts = defaults(opts,struct('duration_days',365.25,'correctionInterval_hr',6, ...
    'IspRCS_s',293,'AZ_m',NaN,'AZmargin_m',0));
validateattributes(opts.altitude_m,{'numeric'},{'scalar','finite','real','nonnegative'});
validateattributes(opts.rho_kgm3,{'numeric'},{'scalar','finite','real','nonnegative'});
validateattributes(opts.duration_days,{'numeric'},{'scalar','finite','real','positive'});
validateattributes(opts.correctionInterval_hr,{'numeric'},{'scalar','finite','real','positive'});
validateattributes(opts.IspRCS_s,{'numeric'},{'scalar','finite','real','positive'});
validateattributes(opts.AZmargin_m,{'numeric'},{'scalar','finite','real','nonnegative'});
assert(isscalar(opts.AZ_m) && (isnan(opts.AZ_m) || ...
    (isfinite(opts.AZ_m) && opts.AZ_m>=0)),'Invalid approach-zone radius.');

mu = 398600.4418e9;
r = 6378.1363e3 + opts.altitude_m;
n = sqrt(mu/r^3);
dragFactor = @(v) v.Cd*v.area_m2/v.mass_kg;
aDiff = -0.5*opts.rho_kgm3*(mu/r)*(dragFactor(follower)-dragFactor(reference));
a = [0;aDiff;0];

A = zeros(6);
A(1:3,4:6) = eye(3);
A(4,1) = 3*n^2;
A(4,5) = 2*n;
A(5,4) = -2*n;
A(6,3) = -n^2;
augmented = [A,[zeros(3);eye(3)];zeros(3,9)];

T = opts.correctionInterval_hr*3600;
duration_s = opts.duration_days*86400;
nCycles = round(duration_s/T);
assert(nCycles>=1 && abs(nCycles*T-duration_s)<1e-6, ...
    'Duration must contain a whole number of correction intervals.');
E = expm(augmented*T);
positionFromVelocity = E(1:3,4:6);
assert(rcond(positionFromVelocity)>1e-10, ...
    'Correction interval is ill-conditioned for two-impulse HCW targeting.');
dv1 = -positionFromVelocity \ (E(1:3,7:9)*a);
dv2 = -(E(4:6,4:6)*dv1 + E(4:6,7:9)*a);
dvInterval = norm(dv1)+norm(dv2);
dv = nCycles*dvInterval;
propellant = follower.mass_kg*(-expm1(-dv/(opts.IspRCS_s*9.80665)));

tLocal = linspace(0,T,121);
errorPosition = zeros(3,numel(tLocal));
for j = 1:numel(tLocal)
    E = expm(augmented*tLocal(j));
    errorPosition(:,j) = E(1:3,4:6)*dv1 + E(1:3,7:9)*a;
end
maxError = max(vecnorm(errorPosition));
strategy = ["Vbar";"orbit"];
minRange = zeros(2,1);
maxRange = zeros(2,1);
violations = nan(2,1);
static = [0;reference.staticDistance_m-follower.staticDistance_m;0];
rho = vecnorm(static+errorPosition);
minRange(1) = min(rho);
maxRange(1) = max(rho);
if isfinite(opts.AZ_m)
    violations(1) = nCycles*nnz(rho<opts.AZ_m+opts.AZmargin_m);
    violations(2) = 0;
end
minRange(2) = Inf;
for k = 1:nCycles
    t = (k-1)*T+tLocal;
    actual = ellipsePosition(follower,n*t)-ellipsePosition(reference,n*t) ...
        + errorPosition;
    rho = vecnorm(actual);
    minRange(2) = min(minRange(2),min(rho));
    maxRange(2) = max(maxRange(2),max(rho));
    if isfinite(opts.AZ_m)
        violations(2) = violations(2)+nnz(rho<opts.AZ_m+opts.AZmargin_m);
    end
end
SK.summary = table(strategy,repmat(dv,2,1),repmat(propellant,2,1), ...
    minRange,maxRange,repmat(maxError,2,1),violations,'VariableNames', ...
    {'Architecture','SK_DV_mps','Propellant_kg','MinRange_m','MaxRange_m', ...
    'MaxPositionError_m','AZViolationSteps'});
SK.annualDV_mps = dv*365.25/opts.duration_days;
SK.differentialDrag_mps2 = aDiff;
SK.correctionDV_mps = dvInterval;
SK.burns_mps = [dv1,dv2];
SK.follower = follower;
SK.reference = reference;
SK.options = opts;
end

function r = ellipsePosition(vehicle,theta)
phase = theta+vehicle.ellipsePhase_rad;
amp = vehicle.ellipseRadialAmplitude_m;
r = [amp*cos(phase);-vehicle.ellipseCenter_m-2*amp*sin(phase);zeros(size(phase))];
end

function v = vehicleInputs(v)
validateattributes(v.mass_kg,{'numeric'},{'scalar','finite','real','positive'});
validateattributes(v.area_m2,{'numeric'},{'scalar','finite','real','nonnegative'});
validateattributes(v.Cd,{'numeric'},{'scalar','finite','real','nonnegative'});
v = defaults(v,struct('staticDistance_m',0,'ellipseCenter_m',0, ...
    'ellipseRadialAmplitude_m',0,'ellipsePhase_rad',0));
validateattributes(v.staticDistance_m,{'numeric'},{'scalar','finite','real'});
validateattributes(v.ellipseCenter_m,{'numeric'},{'scalar','finite','real'});
validateattributes(v.ellipseRadialAmplitude_m,{'numeric'},{'scalar','finite','real','nonnegative'});
validateattributes(v.ellipsePhase_rad,{'numeric'},{'scalar','finite','real'});
end

function s = defaults(s,values)
fields = fieldnames(values);
for k = 1:numel(fields)
    field = fields{k};
    if ~isfield(s,field)
        s.(field) = values.(field);
    end
end
end
