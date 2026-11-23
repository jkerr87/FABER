function varargout = FABER_MA_dynamics(command,varargin)
% FABER_MA_DYNAMICS Shared orbital / relative-motion utilities.
% SI units throughout (m, kg, s). ECI inertial, LVLH x radial / y along-track.
% Functions are routed through one file to minimize path dependencies.
switch lower(char(command))
    case 'defaults'
        varargout{1}=defaults();
    case 'chiefinitial'
        varargout{1}=chiefInitial(varargin{:});
    case 'accelerations'
        varargout{1}=accel(varargin{:});
    case 'frame'
        [varargout{1:nargout}]=frame(varargin{:});
    case 'relative'
        varargout{1}=relative(varargin{:});
    case 'deputyeci'
        varargout{1}=deputyECI(varargin{:});
    case 'phi'
        varargout{1}=phi(varargin{:});
    case 'differential'
        varargout{1}=differential(varargin{:});
    otherwise
        error('FABER_MA_dynamics:Command','Unknown command %s',command);
end
end

function e=defaults()
e.mu_m3s2=3.986004418e14;
e.Re_m=6378136.3;
e.J2=1.08262668e-3;
e.useJ2=true;
e.useDrag=true;
e.useSRP=false;
e.rhoRef_kgm3=1e-12; % USER MUST SELECT altitude/solar-activity case
e.rhoRefAlt_m=680120;
e.scaleHeight_m=60000; % simplified exponential density
e.omegaEarth_rads=7.2921150e-5;
e.atmosphereCorotates=true;
e.sunUnitECI=[1;0;0]; % fixed over each analysis (no ephemeris)
e.srpPressure_Nm2=4.56e-6; % 1 AU approximation
% Properties of deputy and chief; can differ significantly after payload capture.
e.deputy=struct('mass_kg',3503.59,'area_m2',8,'Cd',2.2,'Cr',1.5);
e.chief=struct('mass_kg',1000,'area_m2',8,'Cd',2.2,'Cr',1.5);
end

function chief=chiefInitial(h_km,inc_deg,raan_deg,u_deg,e)
r=e.Re_m+1000*h_km; n=sqrt(e.mu_m3s2/r^3);
i=deg2rad(inc_deg); O=deg2rad(raan_deg); u=deg2rad(u_deg);
R3=[cos(O),-sin(O),0;sin(O),cos(O),0;0,0,1];
R1=[1,0,0;0,cos(i),-sin(i);0,sin(i),cos(i)];
Q=R3*R1;
chief=[r*Q*[cos(u);sin(u);0];sqrt(e.mu_m3s2/r)*Q*[-sin(u);cos(u);0]];
end

function a=accel(x,e,which)
% Cartesian central gravity, J2, exponential drag, optional unshadowed SRP.
r=x(1:3); v=x(4:6); rm=norm(r);
a=-e.mu_m3s2*r/rm^3;
if e.useJ2
    z2r2=(r(3)/rm)^2;
    f=1.5*e.J2*e.mu_m3s2*e.Re_m^2/rm^5;
    a=a+f*[r(1)*(5*z2r2-1);r(2)*(5*z2r2-1);r(3)*(5*z2r2-3)];
end
p=e.(char(which));
if e.useDrag
    h=rm-e.Re_m;
    rho=e.rhoRef_kgm3*exp(-(h-e.rhoRefAlt_m)/e.scaleHeight_m);
    if e.atmosphereCorotates
        vrel=v-cross([0;0;e.omegaEarth_rads],r);
    else
        vrel=v;
    end
    a=a-0.5*rho*p.Cd*p.area_m2/p.mass_kg*norm(vrel)*vrel;
end
if e.useSRP
    sun=e.sunUnitECI(:)/norm(e.sunUnitECI);
    a=a-e.srpPressure_Nm2*p.Cr*p.area_m2/p.mass_kg*sun;
end
end

function [Q,omega]=frame(chief)
% Q columns are LVLH radial, along-track, orbit-normal in ECI.
r=chief(1:3); v=chief(4:6);
er=r/norm(r); eh=cross(r,v); eh=eh/norm(eh);
et=cross(eh,er); Q=[er,et,eh];
omega=[0;0;norm(cross(r,v))/norm(r)^2]; % LVLH angular rate
end

function dep=deputyECI(chief,X)
[Q,omega]=frame(chief);
r=X(1:3); v=X(4:6);
dep=[chief(1:3)+Q*r;chief(4:6)+Q*(v+cross(omega,r))];
end

function X=relative(chief,dep)
[Q,omega]=frame(chief);
r=Q'*(dep(1:3)-chief(1:3));
v=Q'*(dep(4:6)-chief(4:6))-cross(omega,r);
X=[r;v];
end

function P=phi(n,t)
s=sin(n*t); c=cos(n*t); nt=n*t;
P=[4-3*c,0,0,s/n,2*(1-c)/n,0; ...
   6*(s-nt),1,0,-2*(1-c)/n,(4*s-3*nt)/n,0; ...
   0,0,c,0,0,s/n; ...
   3*n*s,0,0,c,2*s,0; ...
   6*n*(c-1),0,0,-2*s,4*c-3,0; ...
   0,0,-n*s,0,0,c];
end

function aL=differential(t,Xref,e,h_km,i_deg,raan_deg,u0_deg)
% Short-arc differential NONCENTRAL perturbations, projected in a circular
% chief LVLH reference. Central gravity is ALREADY in the HCW dynamics.
% For long parking intervals use full chief/deputy ECI propagation instead.
chief0=chiefInitial(h_km,i_deg,raan_deg,u0_deg,e);
r0=norm(chief0(1:3)); n=sqrt(e.mu_m3s2/r0^3);
chief=chiefInitial(h_km,i_deg,raan_deg,u0_deg+rad2deg(n*t),e);
dep=deputyECI(chief,Xref);
[Q,~]=frame(chief);
e0=e; e0.useJ2=false; e0.useDrag=false; e0.useSRP=false;
ad=accel(dep,e,'deputy')-accel(dep,e0,'deputy');
ac=accel(chief,e,'chief')-accel(chief,e0,'chief');
aL=Q'*(ad-ac);
end
