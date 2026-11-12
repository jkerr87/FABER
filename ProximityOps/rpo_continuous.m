
function result = rpo_continuous(X0, Xf, n, T_s, axisName, N)
% RPO_CONTINUOUS
%
% Continuous-control final approach using quintic trajectory planning
% and HCW feedforward acceleration.
%
% COORDINATES:
% x = radial OUTWARD
% y = along-track (V-bar)
% z = cross-track
%
% INPUTS:
% X0       Initial state [x;y;z;vx;vy;vz] [m, m/s]
% Xf       Desired final state [x;y;z;vx;vy;vz] [m, m/s]
% n        Mean orbital motion [rad/s]
% T_s      Final approach duration [s]
% axisName 'Vbar', 'Rbar', or 'general'
% N        Number of trajectory samples
%
% OUTPUTS:
% result.DV_mps
% result.TOF_s
% result.TOF_min
% result.trajectory
% result.finalState
% result.acceleration_mps2
% result.peakAcceleration_mps2
% result.maxSpeed_mps
% result.minRange_m
%
% ASSUMPTIONS:
% Circular chief orbit and unperturbed HCW dynamics.
% Ideal continuous acceleration with perfect trajectory tracking.
% No finite-thrust limitations, navigation errors, or collision avoidance.
%
% The approach is planned to end at the specified Xf.
% For capture holding, choose a suitable final position with vf = 0.

%% INPUTS

if nargin < 6 || isempty(N)
    N = 600;
end

X0 = X0(:);
Xf = Xf(:);

assert(numel(X0)==6 && numel(Xf)==6, ...
    'X0 and Xf must contain six state components.');

assert(n > 0, 'Mean motion must be positive.');
assert(T_s > 0, 'Transfer duration must be positive.');
assert(N >= 3, 'At least three samples are required.');

r0 = X0(1:3);
v0 = X0(4:6);

rf = Xf(1:3);
vf = Xf(4:6);

%% APPROACH GEOMETRY VALIDATION

name = lower(string(axisName));

tol = 1e-8;

switch name

    case "vbar"

        % Strict V-bar: radial and cross-track states must be zero.
        assert(all(abs(X0([1 3 4 6])) < tol) && ...
               all(abs(Xf([1 3 4 6])) < tol), ...
            'V-bar approach requires zero off-axis position and velocity.');

    case "rbar"

        % Strict R-bar: along-track and cross-track states must be zero.
        assert(all(abs(X0([2 3 5 6])) < tol) && ...
               all(abs(Xf([2 3 5 6])) < tol), ...
            'R-bar approach requires zero off-axis position and velocity.');

    case "general"

        % Arbitrary 3D position and velocity boundary conditions.

    otherwise

        error('axisName must be Vbar, Rbar, or general.');

end

%% QUINTIC TRAJECTORY

% Normalized time avoids large powers of T in the coefficient matrix.

t = linspace(0,T_s,N);
s = t/T_s;

r = zeros(3,N);
v = zeros(3,N);
a = zeros(3,N);

% Prescribed kinematic accelerations at the endpoints.
a0 = zeros(3,1);
af = zeros(3,1);

% Polynomial:
% r(s) = c0+c1*s+c2*s^2+c3*s^3+c4*s^4+c5*s^5
%
% Derivatives:
% v = (dr/ds)/T
% a = (d2r/ds2)/T^2

M = [
    1 0 0 0 0 0;
    0 1 0 0 0 0;
    0 0 2 0 0 0;
    1 1 1 1 1 1;
    0 1 2 3 4 5;
    0 0 2 6 12 20
];

for k = 1:3

    % Boundary conditions in normalized time
    b = [
        r0(k);
        v0(k)*T_s;
        a0(k)*T_s^2;
        rf(k);
        vf(k)*T_s;
        af(k)*T_s^2
    ];

    c = M\b;

    % Position
    r(k,:) = ...
        c(1) + c(2)*s + c(3)*s.^2 + ...
        c(4)*s.^3 + c(5)*s.^4 + c(6)*s.^5;

    % Velocity
    v(k,:) = ( ...
        c(2) + 2*c(3)*s + 3*c(4)*s.^2 + ...
        4*c(5)*s.^3 + 5*c(6)*s.^4 )/T_s;

    % Acceleration
    a(k,:) = ( ...
        2*c(3) + 6*c(4)*s + ...
        12*c(5)*s.^2 + 20*c(6)*s.^3 )/T_s^2;

end

%% HCW CONTINUOUS CONTROL ACCELERATION

% HCW equations:
%
% xddot =  3*n^2*x + 2*n*ydot + ux
% yddot = -2*n*xdot + uy
% zddot = -n^2*z + uz
%
% Rearranged to obtain the required control acceleration.

u = [
    a(1,:) - 3*n^2*r(1,:) - 2*n*v(2,:);
    a(2,:) + 2*n*v(1,:);
    a(3,:) + n^2*r(3,:)
];

%% DELTA V AND PERFORMANCE

accelMagnitude = vecnorm(u,2,1);
speedMagnitude = vecnorm(v,2,1);
range = vecnorm(r,2,1);

% Equivalent continuous-control delta V
DV = trapz(t,accelMagnitude);

% Peak required acceleration
peakAcceleration = max(accelMagnitude);

% Maximum relative speed
maxSpeed = max(speedMagnitude);

% Minimum target separation
minRange = min(range);

%% TERMINAL STATE VERIFICATION

Xend = [r(:,end);v(:,end)];

posError = norm(Xend(1:3)-rf);
velError = norm(Xend(4:6)-vf);

%% OUTPUT STRUCTURE

result.method = char(name);

result.initialState = X0;
result.desiredFinalState = Xf;
result.finalState = Xend;

result.DV_mps = DV;

result.T_s = T_s;
result.TOF_s = T_s;
result.TOF_min = T_s/60;

result.t_s = t;
result.trajectory = [r;v];

result.acceleration_mps2 = u;
result.accelerationMagnitude_mps2 = accelMagnitude;

result.peakAcceleration_mps2 = peakAcceleration;
result.maxSpeed_mps = maxSpeed;
result.minRange_m = minRange;

result.posError_m = posError;
result.velError_mps = velError;

result.feasible = posError < 1e-5 && velError < 1e-7;

end
