function x = FABER_environment_eci(chief,relative)
% Convert a sampled LVLH history into its moving chief's inertial frame.
r = chief(1:3,:);
v = chief(4:6,:);
rmag = vecnorm(r);
R = r./rmag;
H = cross(r,v,1);
Hmag = vecnorm(H);
N = H./Hmag;
T = cross(N,R,1);
motion = Hmag./rmag.^2;
p = relative(1:3,:);
u = relative(4:6,:) + [-motion.*p(2,:); motion.*p(1,:); zeros(1,size(p,2))];
x = [r+R.*p(1,:)+T.*p(2,:)+N.*p(3,:); ...
     v+R.*u(1,:)+T.*u(2,:)+N.*u(3,:)];
end
