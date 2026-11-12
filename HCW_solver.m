function result = HCW_solver(r0, v0, rf, vf, T, n, N)
% -------------------------------------------------------------------------
% Function to use the HCW solver for relative motion. This general function
% allows us to use the solver for different mission phases and possible
% degrees of freedom.
%
% INPUTS:
% r0 - Initial relative position [3x1], m
% v0 - Initial relative velocity [3x1], m/s
% rf - Desired final relative position [3x1], m
% vf - Desired final relative velocity [3x1], m/s
% T  - Transfer duration, s
% n  - Chief orbital mean motion, rad/s
%
% OUTPUTS:
% result.dv1        - First burn magnitude [m/s]
% result.dv2        - Second burn magnitude [m/s]
% result.dvTotal    - Total maneuver delta-v [m/s]
% result.burn1      - First burn vector [3x1], m/s
% result.burn2      - Second burn vector [3x1], m/s
% result.t          - Trajectory time vector [s]
% result.trajectory - Coast state history [6xN]
% result.finalState - State after second impulse [6x1]
%--------------------------------------------------------------------------
%% HCW linear relative equations of motion
% xddot = 3*n^2*x + 2*n*ydot + u_x
% yddot = -2*n*xdot + u_y
% zddot = -n^2*z + u_z
% -----------------------------------------------------------------------

%% DEFINE INPUTS
% state vector: X = [x y z xdot ydot zdot]'

% recieve states
r0 = r0(:);
v0 = v0(:);
rf = rf(:);
vf = vf(:);

% Continuous-time HCW state-space model
A = [0 0 0 1 0 0; ...
    0 0 0 0 1 0; ...
    0 0 0 0 0 1; ...
    3*n^2 0 0 0 2*n 0; ...
    0 0 0 -2*n 0 0; ... 
    0 0 -n^2 0 0 0];

% STM - state transition matrix X(T) = Phi(T)*X(0), Phi(T) = e^AT
Phi = expm(A*T);

Prr = Phi(1:3,1:3);             % initial position → final position
Prv = Phi(1:3,4:6);             % initial velocity → final position
Pvr = Phi(4:6,1:3);             % initial position → final velocity
Pvv = Phi(4:6,4:6);             % initial velocity → final velocity

%% SOLVE
% solve backwards from the desired final position (first impulse):
vAfter1 = Prv \ (rf-Prr*r0);
burn1 = vAfter1-v0;
DV1 = norm(burn1);

% velocity right before AZ
vBefore2 = Pvr*r0 + Pvv*vAfter1;

% second impulse
burn2 = vf-vBefore2;
DV2 = norm(burn2);

DVtotal = DV1+DV2;

%% PROPAGATE COAST TRAJECTORY
t = linspace(0,T,N);

Xafter1 = [r0;vAfter1];
trajectory = zeros(6,N);

for k = 1:N
    Phi_t = expm(A*t(k));
    trajectory(:,k) = Phi_t*Xafter1;
end

% Final state after terminal burn
finalState = trajectory(:,end);
finalState(4:6) = finalState(4:6) + burn2;

%% OUTPUT RESULTS

result.dv1 = DV1;
result.dv2 = DV2;
result.dvTotal = DVtotal;

result.burn1 = burn1;
result.burn2 = burn2;

result.t = t;
result.trajectory = trajectory;
result.finalState = finalState;

% Verify final state accuracy
result.positionError = norm(finalState(1:3)-rf);
result.velocityError = norm(finalState(4:6)-vf);

r = result.trajectory(1:3,:);
v = result.trajectory(4:6,:);

range = vecnorm(r,2,1);
relSpeed = vecnorm(v,2,1);

closingSpeed = -sum(r.*v,1)./range;

result.relativeSpeed = relSpeed;
result.closingSpeed = closingSpeed;
result.maxRelativeSpeed = max(relSpeed);
result.maxClosingSpeed = max(closingSpeed);

result.ok=all(isfinite(r(:))) && all(isfinite(v(:))) && isfinite(DVtotal);

end