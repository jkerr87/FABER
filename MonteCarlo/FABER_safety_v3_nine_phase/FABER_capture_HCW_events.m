function events = FABER_capture_HCW_events(x0,xf,T_s,n,phase,vehicle,legIndex,step_s,burn1,burn2)
% FABER_CAPTURE_HCW_EVENTS Capture ACTUAL accepted HCW burn events.
%
% Pass the exact two burn vectors returned by the nominal transfer solver.
% This avoids re-solving a potentially singular HCW boundary value problem,
% particularly for analytical V-bar/R-bar trajectories whose duration is
% chosen by analytical_manos rather than the sampled design duration.
%
% x0,xf: 6x1 relative HCW states [m; m/s]
% T_s:    actual accepted coast duration [s]
% n:      mean motion [rad/s]
% burn1, burn2: 3x1 actual accepted impulsive DV vectors [m/s]
%
% Optionally, the last two inputs may be omitted for older callers. In
% that case a nonsingular two-impulse HCW boundary-value solution is used.

if nargin<8 || isempty(step_s), step_s=600; end
assert(numel(x0)==6 && numel(xf)==6 && isscalar(T_s) && ...
    T_s>0 && step_s>0,'Invalid event capture inputs');
x0=x0(:); xf=xf(:);
assert(all(isfinite([x0;xf;n;T_s;step_s])));

if nargin<10 || isempty(burn1) || isempty(burn2)
    % Legacy fallback ONLY. Prefer passing the accepted solver burns.
    [A,B,C,D]=phi(n,T_s);
    if rcond(B)<=1e-10
        error('FABER:CaptureNeedsBurnVectors', ...
            ['HCW boundary-value matrix ill-conditioned. ', ...
             'Pass actual accepted departure/arrival burns.']);
    end

    v1=B\(xf(1:3)-A*x0(1:3));
    burn1=v1-x0(4:6);
    burn2=xf(4:6)-(C*x0(1:3)+D*v1);
end

burn1=burn1(:);burn2=burn2(:);
assert(numel(burn1)==3 && numel(burn2)==3 && ...
       all(isfinite([burn1;burn2])), 'Invalid accepted burn vectors');

xPost=x0;
xPost(4:6)=xPost(4:6)+burn1;

% Coast at the exact accepted duration. No inversion is performed.
[Ae,Be,Ce,De]=phi(n,T_s);
xBeforeArrival=[Ae Be;Ce De]*xPost;

% Numerical consistency checks; inconsistent solver data must not be
% silently saved as an exact nominal trajectory.
posErr=norm(xBeforeArrival(1:3)-xf(1:3));
velErr=norm(xBeforeArrival(4:6)+burn2-xf(4:6));
if posErr>1e-3 || velErr>1e-6
    error('FABER:CaptureInconsistentNominal', ...
        ['Accepted HCW burns do not reproduce the nominal endpoint ', ...
         '(position error %.6g m, velocity error %.6g m/s).'], ...
        posErr,velErr);
end

events=struct('xPre',{},'dvNom',{},'phase',{},'vehicle',{},'index',{});
events(end+1)=mk(x0,burn1,[char(phase) ' start'],vehicle,legIndex);
for tau=step_s:step_s:(T_s-1e-8)
    [A,B,C,D]=phi(n,tau);
    xCoast=[A B;C D]*xPost;
    events(end+1)=mk(xCoast,nan(3,1), ...
        [char(phase) ' coast'],vehicle,legIndex);
end
events(end+1)=mk(xBeforeArrival,burn2, ...
    [char(phase) ' end'],vehicle,legIndex);
end

function e=mk(x,dv,p,v,i)
e=struct('xPre',x,'dvNom',dv,'phase',p, ...
    'vehicle',v,'index',i);
end

function [A,B,C,D] = phi(n,t)
% HCW state transition matrix partitioned into 3x3 blocks.
%
% x(t) = [A B; C D] * x(0)
% State ordering: [x y z vx vy vz]'

s = sin(n*t);
c = cos(n*t);
nt = n*t;

A = [ ...
    4-3*c,       0, 0;
    6*(s-nt),    1, 0;
    0,           0, c];

B = [ ...
    s/n,           2*(1-c)/n,       0;
    -2*(1-c)/n,   (4*s-3*nt)/n,    0;
    0,             0,              s/n];

C = [ ...
    3*n*s,       0, 0;
    6*n*(c-1),   0, 0;
    0,           0, -n*s];

D = [ ...
    c,       2*s,     0;
    -2*s,    4*c-3,   0;
    0,       0,       c];

end
