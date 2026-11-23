
function D = FABER_Taxi_initial_departure(cfg)
% FABER_TAXI_INITIAL_DEPARTURE
%
% Preliminary continuous-thrust deployment from ATLANTE.
%
% Frame: ATLANTE-centred LVLH
% x = radial, y = along-track, z = cross-track
%
% No GNC controller or attitude propagation.
% Required acceleration is calculated from a quintic
% reference and propagated with HCW dynamics.
%
% This is a separate deployment maneuver, NOT
% the Taxi-to-client homing trajectory.

%% 1. Configuration

mu = cfg.mu_km3s2;
Re = cfg.Re_km;
h  = cfg.h_km;

n = sqrt(mu/(Re+h)^3);

g0 = 9.80665;

% Initial and final relative states
X0 = cfg.departure.initialState(:);
Xf = cfg.departure.finalState(:);

% Transfer duration [s]
T = cfg.departure.TOF_min * 60;

% Initial mass
m0 = cfg.massTaxi_kg;

% Selected propulsion group
pcfg = cfg;

switch lower(cfg.departure.propulsion)

    case 'main'
        sys = cfg.propulsion.main;

    case 'rcs'
        sys = cfg.propulsion.rcs;

    otherwise
        error('Unknown departure propulsion system.');
end

pcfg.thrust_N = sys.thrust_N;
pcfg.nActiveThrusters = sys.nActiveThrusters;
pcfg.thrustEfficiency = sys.thrustEfficiency;
pcfg.Isp_s = sys.Isp_s;

% Effective translational thrust.
% Assumes selected engines can deliver force in
% the required direction.
Fmax = pcfg.thrust_N * ...
       pcfg.nActiveThrusters * ...
       pcfg.thrustEfficiency;

ve = g0 * pcfg.Isp_s;

%% 2. Quintic reference coefficients

r0 = X0(1:3);
v0 = X0(4:6);

rf = Xf(1:3);
vf = Xf(4:6);

Dr = rf-r0-v0*T;
Dv = (vf-v0)*T;

A = 10*Dr-4*Dv;
B = -15*Dr+7*Dv;
C = 6*Dr-3*Dv;

%% 3. Numerical propagation

N = 301;
tGrid = linspace(0,T,N);

Y0 = [X0;m0];

opts = odeset( ...
    'RelTol',1e-9, ...
    'AbsTol',1e-10);

[t,Y] = ode45( ...
    @(t,y) dynamics(t,y), ...
    tGrid,Y0,opts);

%% 4. Thrust history and feasibility

U = zeros(numel(t),3);
F = zeros(numel(t),3);

for k = 1:numel(t)

    u = requiredAccel(t(k));

    U(k,:) = u(:)';
    F(k,:) = Y(k,7)*u(:)';

end

Fmag = vecnorm(F,2,2);
Umag = vecnorm(U,2,2);

rho = vecnorm(Y(:,1:3),2,2);

rhoDot = sum( ...
    Y(:,1:3).*Y(:,4:6),2) ./ max(rho,eps);

speed = vecnorm(Y(:,4:6),2,2);

% Total acceleration-equivalent Delta-V
DV = trapz(t,Umag);

% Propellant mass
propellant = m0-Y(end,7);

% Actual thrust-active duration
active = Fmag > 1e-9;

burnTime = trapz(t,double(active));

% Preliminary feasibility checks
thrustOK = all(Fmag <= Fmax+1e-8);

speedOK = all( ...
    speed <= cfg.departure.maxSpeed_mps+1e-8);

clearanceOK = all( ...
    rho >= cfg.departure.minSeparation_m-1e-6);

% Check propagated terminal-state agreement
posError = norm(Y(end,1:3)'-rf);
velError = norm(Y(end,4:6)'-vf);

terminalOK = ...
    posError <= cfg.departure.posTolerance_m && ...
    velError <= cfg.departure.velTolerance_mps;

%% 5. Return results

D = struct();

D.feasible = ...
    thrustOK && speedOK && clearanceOK && terminalOK;

D.thrustFeasible = thrustOK;
D.speedFeasible = speedOK;
D.clearanceFeasible = clearanceOK;
D.terminalFeasible = terminalOK;

D.time_s = t;
D.state = Y(:,1:6);
D.mass_kg = Y(:,7);

D.accel_LVLH_mps2 = U;
D.thrust_LVLH_N = F;

D.rho_m = rho;
D.rhoDot_mps = rhoDot;

D.DV_mps = DV;
D.TOF_min = T/60;

D.propellant_kg = propellant;
D.burnTime_s = burnTime;

D.maxThrust_N = max(Fmag);
D.availableThrust_N = Fmax;

D.maxSpeed_mps = max(speed);
D.minRange_m = min(rho);

D.positionError_m = posError;
D.velocityError_mps = velError;

D.initialState = X0;
D.finalState = Xf;

%% 6. Nested dynamics

    function dy = dynamics(t,y)

        r = y(1:3);
        v = y(4:6);
        m = y(7);

        u = requiredAccel(t);

        % Free HCW acceleration
        aHCW = [ ...
            3*n^2*r(1)+2*n*v(2); ...
            -2*n*v(1); ...
            -n^2*r(3)];

        Freq = m*norm(u);

        mdot = -Freq/ve;

        dy = [v;aHCW+u;mdot];

    end

%% 7. Prescribed continuous acceleration

    function u = requiredAccel(t)

        s = t/T;

        r = r0 + v0*t + ...
            A*s^3+B*s^4+C*s^5;

        v = v0 + ...
            (3*A*s^2+4*B*s^3+5*C*s^4)/T;

        a = ( ...
            6*A*s+12*B*s^2+20*C*s^3)/T^2;

        % Required acceleration after HCW compensation
        u = [ ...
            a(1)-3*n^2*r(1)-2*n*v(2); ...
            a(2)+2*n*v(1); ...
            a(3)+n^2*r(3)];

    end

end
