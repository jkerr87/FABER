
function result = rpo_holdingpt(n, T_holding, mode, param, N)
% RPO_HOLDINGPT
%
% Propagates a stationary holding point or bounded relative ellipse
% using the unperturbed HCW equations.
%
% Coordinate convention:
% x = radial OUTWARD
% y = along-track (V-bar)
% z = cross-track
%
% INPUTS:
% n          Mean orbital motion [rad/s]
% T_holding  Holding duration [s]
% mode       'point' or 'ellipse'
% param      Holding configuration parameters
% N          Number of propagation samples
%
% POINT:
% param.yHold  Along-track holding position [m]
%
% ELLIPSE:
% param.yCenter  Along-track centre [m]
% param.radAmp   Radial amplitude [m]
% param.phase    Initial phase [rad]
%
% OUTPUTS:
% DV_holding_mps   Ideal holding delta-V [m/s]
% TOF_holding_min  Holding duration [min]
% holdTrajectory   Relative state history [6 x N]
% holdFinalState   Final relative state [6 x 1]
%
% Assumptions:
% Circular reference orbit, ideal two-body HCW dynamics.
% No perturbations, stationkeeping corrections or collision avoidance.

%% Defaults and validation

if nargin < 5
    N = 600;
end

assert(n > 0, 'Mean motion must be positive.');
assert(T_holding >= 0, 'Holding time cannot be negative.');
assert(N >= 2, 'N must be at least 2.');

%% Define holding configuration

switch lower(mode)

    case 'point'

        % Stationary holding point on the V-bar
        rHold = [0; param.yHold; 0];
        vHold = zeros(3,1);

    case 'ellipse'

        % Bounded relative ellipse in the orbital plane
        phi = 0;

        if isfield(param,'phase')
            phi = param.phase;
        end

        A = param.radAmp;
        yc = param.yCenter;

        rHold = [ ...
             A*cos(phi);
             yc - 2*A*sin(phi);
             0];

        vHold = [ ...
            -n*A*sin(phi);
            -2*n*A*cos(phi);
             0];

    otherwise

        error('Mode must be point or ellipse.');

end

%% HCW state-space matrix

M = [zeros(3), eye(3);
     3*n^2, 0, 0, 0, 2*n, 0;
     0, 0, 0, -2*n, 0, 0;
     0, 0, -n^2, 0, 0, 0];

%% Propagate holding configuration

tCoast = linspace(0,T_holding,N);

coast = zeros(6,N);

Xhold = [rHold;vHold];

for j = 1:N

    Phi = expm(M*tCoast(j));

    coast(:,j) = Phi*Xhold;

end

%% Holding delta-V and time

% Ideal passive HCW holding
DV_holding = 0;

TOF_holding_s = T_holding;
TOF_holding_min = T_holding/60;

%% Geometry

ranges = vecnorm(coast(1:3,:),2,1);

%% Save outputs

result.mode = mode;

result.holdState = Xhold;
result.holdTrajectory = coast;
result.holdTime_s = tCoast;

result.holdFinalState = coast(:,end);

result.DV_holding_mps = DV_holding;

result.TOF_holding_s = TOF_holding_s;
result.TOF_holding_min = TOF_holding_min;

result.minHoldRange_m = min(ranges);
result.maxHoldRange_m = max(ranges);

end
