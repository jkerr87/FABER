function x = FABER_environment_circular(mu,x0,t)
% Circular chief orbit from an inertial state, evaluated at elapsed seconds.
r = x0(1:3);
v = x0(4:6);
n = sqrt(mu/norm(r)^3);
a = n*t(:)';
x = [r*cos(a)+v/n*sin(a); v*cos(a)-r*n*sin(a)];
end
