function [fraction,aspect,beta] = FABER_environment_light(model,t,x,axisLVLH,Re_m)
% Occulted fraction uses the overlap of the apparent Sun and Earth disks.
s = model.sun(t)-x(1:3,:);
su = s./vecnorm(s);
r = x(1:3,:);
rn = vecnorm(r);
radial = r./rn;
separation = acos(clamp(sum(-radial.*su,1)));
earth = asin(min(1,Re_m./rn));
solar = asin(min(1,695700000./vecnorm(s)));
fraction = ones(size(separation));
covered = separation<=earth-solar;
fraction(covered) = 0;
overlap = ~covered & separation<earth+solar;
d = separation(overlap); a = earth(overlap); b = solar(overlap);
area = a.^2.*acos(clamp((d.^2+a.^2-b.^2)./(2*d.*a))) + ...
    b.^2.*acos(clamp((d.^2+b.^2-a.^2)./(2*d.*b))) - ...
    0.5*sqrt(max(0,(-d+a+b).*(d+a-b).*(d-a+b).*(d+a+b)));
fraction(overlap) = 1-area./(pi*b.^2);
fraction = max(0,min(1,fraction));
fraction(rn<=Re_m) = 0;
if nargout>1
    normal = cross(r,x(4:6,:),1);
    normal = normal./vecnorm(normal);
    transverse = cross(normal,radial,1);
    axis = axisLVLH(1)*radial+axisLVLH(2)*transverse+axisLVLH(3)*normal;
    aspect = acosd(clamp(sum(axis.*su,1)));
    beta = asind(clamp(sum(normal.*su,1)));
end
end

function x = clamp(x)
x = min(1,max(-1,x));
end
