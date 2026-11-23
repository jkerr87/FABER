function state = FABER_relative_state(chief,vehicle,direction,chiefAccel)
% Convert inertial and rotating RTN/LVLH states (m and m/s, same epoch).
% eciToLvlh: vehicle is absolute inertial. lvlhToEci: vehicle is relative.
% Supply chief acceleration [m/s^2] for out-of-plane motion of the RTN frame.
validateattributes(chief,{'numeric'},{'real','finite','numel',6});
validateattributes(vehicle,{'numeric'},{'real','finite','numel',6});
chief = chief(:);
vehicle = vehicle(:);
r = chief(1:3);
v = chief(4:6);
h = cross(r,v);
assert(norm(r)>0 && norm(h)>0,'Chief RTN frame is undefined.');
radial = r/norm(r);
normal = h/norm(h);
along = cross(normal,radial);
basis = [radial,along,normal];
omega = [0;0;norm(h)/norm(r)^2];
if nargin>=4 && ~isempty(chiefAccel)
    validateattributes(chiefAccel,{'numeric'},{'real','finite','numel',3});
    omega(1) = norm(r)*dot(chiefAccel(:),normal)/norm(h);
end
switch validatestring(direction,{'eciToLvlh','lvlhToEci'})
    case 'eciToLvlh'
        rho = basis'*(vehicle(1:3)-r);
        state = [rho;basis'*(vehicle(4:6)-v)-cross(omega,rho)];
    case 'lvlhToEci'
        rho = vehicle(1:3);
        state = [r+basis*rho;v+basis*(vehicle(4:6)+cross(omega,rho))];
end
end
