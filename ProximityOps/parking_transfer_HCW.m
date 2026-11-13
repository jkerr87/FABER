
function result = parking_transfer_HCW(mu, Re, h_km, ...
    xInitial, distance_m, architecture, ...
    radialAmplitude_m, transferTime_s)
% PARKING_TRANSFER_HCW
%
% Two-impulse HCW transfer into a V-bar point or relative orbit.
%
% xInitial: [x;y;z;vx;vy;vz], SI units
% distance_m: negative V-bar parking distance magnitude
% architecture: 'Vbar' or 'orbit'
% radialAmplitude_m: HCW orbit radial amplitude
% transferTime_s: total maneuver duration [s]
%
% Ideal impulsive transfer. No KOZ, finite-thrust, sensor,
% or passive-safety checks are implemented.

n = sqrt(mu/(Re+h_km)^3);

% HCW state matrix (radial, along-track, cross-track)
A = [ ...
    0 0 0 1 0 0;
    0 0 0 0 1 0;
    0 0 0 0 0 1;
    3*n^2 0 0 0 2*n 0;
    0 0 0 -2*n 0 0;
    0 0 -n^2 0 0 0];

% Target state at the end of the transfer
if strcmpi(architecture,'Vbar')

    xTarget = [0;-distance_m;0;0;0;0];

elseif strcmpi(architecture,'orbit')

    amp = radialAmplitude_m;

    theta = n*transferTime_s;

    xTarget = [ ...
        amp*cos(theta);
        -distance_m-2*amp*sin(theta);
        0;
        -amp*n*sin(theta);
        -2*amp*n*cos(theta);
        0];

else
    error('Architecture must be Vbar or orbit.');
end

% State transition matrix
Phi = expm(A*transferTime_s);

Phi_rr = Phi(1:3,1:3);
Phi_rv = Phi(1:3,4:6);

Phi_vr = Phi(4:6,1:3);
Phi_vv = Phi(4:6,4:6);

r0 = xInitial(1:3);
v0 = xInitial(4:6);

rf = xTarget(1:3);
vf = xTarget(4:6);

if rcond(Phi_rv) < 1e-10
    error('Ill-conditioned HCW transfer time.');
end

% Initial velocity required to reach target position
vRequired = Phi_rv\(rf-Phi_rr*r0);

dv1 = vRequired-v0;

% Velocity before final matching impulse
vBeforeFinal = Phi_vr*r0 + Phi_vv*vRequired;

dv2 = vf-vBeforeFinal;

result.DV_initial_mps = norm(dv1);
result.DV_final_mps = norm(dv2);

result.DV_total_mps = norm(dv1)+norm(dv2);
result.TOF_min = transferTime_s/60;

result.xTarget = xTarget;
result.dv1 = dv1;
result.dv2 = dv2;

end
