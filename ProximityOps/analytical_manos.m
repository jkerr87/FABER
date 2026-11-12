function result = analytical_manos(r0,rf,n,method,N,plotFlag)
% ANALYTICAL_MANOS
%
% Analytical two-impulse closing maneuvers:
%   1. V-bar closing
%   2. R-bar closing
%
% Uses analytical maneuver equations and HCW propagation.
%
% Coordinate convention:
% x = radial OUTWARD
% y = along-track (V-bar)
% z = cross-track
%
% INPUTS:
% r0       Initial relative position [3x1] [m]
% rf       Final relative position [3x1] [m]
% n        Mean motion [rad/s]
% method   'Vbar' or 'Rbar'
% N        Number of trajectory samples
% plotFlag true/false
%
% ASSUMPTIONS:
% - Circular reference orbit
% - Stationary initial and final holding points
% - Initial and final velocities equal zero
% - Both holding points lie on the V-bar
% - Ideal impulsive maneuvers
% - No perturbations or collision avoidance

%% Inputs and defaults

if nargin < 5 || isempty(N)
    N = 600;
end

if nargin < 6 || isempty(plotFlag)
    plotFlag = false;
end

r0 = r0(:);
rf = rf(:);

v0 = zeros(3,1);
vf = zeros(3,1);

assert(n > 0);
assert(N >= 2);
assert(numel(r0)==3 && numel(rf)==3);

% Analytical formulas require V-bar holding points
tol = 1e-9;

assert(all(abs(r0([1 3])) < tol) && ...
       all(abs(rf([1 3])) < tol), ...
       'Analytical maneuver requires V-bar holding points.');

%% Orbital parameters

P = 2*pi/n;

% Signed along-track displacement
deltaY = rf(2)-r0(2);

%% Select analytical maneuver

switch lower(method)

    case 'vbar'

        % Two impulses along V-bar
        % Transfer duration: one orbital period

        T = P;

        burn1 = [0; -n*deltaY/(6*pi); 0];
        burn2 = -burn1;

    case 'rbar'

        % Two impulses along radial axis
        % Transfer duration: half orbital period
        % Note: positive x points OUTWARD

        T = P/2;

        burn1 = [-n*deltaY/4; 0; 0];
        burn2 = burn1;

    otherwise

        error('Method must be Vbar or Rbar.');

end

%% Analytical delta V

DV1 = norm(burn1);
DV2 = norm(burn2);

DV_total = DV1 + DV2;

%% HCW propagation

A = [zeros(3), eye(3);
     3*n^2, 0, 0, 0, 2*n, 0;
     0, 0, 0, -2*n, 0, 0;
     0, 0, -n^2, 0, 0, 0];

% State immediately after first impulse
X0plus = [r0; v0+burn1];

% Propagate coast
t = linspace(0,T,N);

trajectory = zeros(6,N);

for k = 1:N

    Phi = expm(A*t(k));

    trajectory(:,k) = Phi*X0plus;

end

%% Apply second impulse

Xfinal = trajectory(:,end);

Xfinal(4:6) = Xfinal(4:6)+burn2;

%% Verify terminal conditions

posError = norm(Xfinal(1:3)-rf);
velError = norm(Xfinal(4:6)-vf);

feasible = posError < 1e-5 && velError < 1e-7;

%% Additional trajectory information

range = vecnorm(trajectory(1:3,:),2,1);
speed = vecnorm(trajectory(4:6,:),2,1);

%% Save outputs

result.method = method;

result.DV1_mps = DV1;
result.DV2_mps = DV2;
result.DV_mps = DV_total;

result.TOF_s = T;
result.TOF_min = T/60;

result.burn1 = burn1;
result.burn2 = burn2;

result.t_s = t;
result.trajectory = trajectory;

result.initialState = [r0;v0];
result.finalState = Xfinal;

result.posError_m = posError;
result.velError_mps = velError;

result.minRange_m = min(range);
result.maxSpeed_mps = max(speed);

result.feasible = feasible;

%% Optional plotting

if plotFlag

    x = trajectory(1,:);
    y = trajectory(2,:);
    z = trajectory(3,:);

    figure('Color','w', ...
           'Name',['Analytical ' method]);

    tiledlayout(1,2);

    % 2D
    nexttile;
    hold on; grid on; axis equal;

    plot(y,x,'b-','LineWidth',1.8);

    plot(0,0,'kp','MarkerFaceColor','k');
    plot(r0(2),r0(1),'go','MarkerFaceColor','g');
    plot(rf(2),rf(1),'ro','MarkerFaceColor','r');

    xlabel('Along-track y [m]');
    ylabel('Radial outward x [m]');
    title([method ' Closing - 2D']);

    % 3D
    nexttile;
    hold on; grid on; axis equal;

    plot3(y,x,z,'b-','LineWidth',1.8);

    plot3(0,0,0,'kp','MarkerFaceColor','k');
    plot3(r0(2),r0(1),r0(3), ...
          'go','MarkerFaceColor','g');
    plot3(rf(2),rf(1),rf(3), ...
          'ro','MarkerFaceColor','r');

    xlabel('Along-track y [m]');
    ylabel('Radial outward x [m]');
    zlabel('Cross-track z [m]');

    title([method ' Closing - 3D']);
    view(35,25);

end

end
