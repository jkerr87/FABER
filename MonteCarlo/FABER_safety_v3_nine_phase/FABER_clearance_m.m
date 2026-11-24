function clearance = FABER_clearance_m(r,hazards,extraRadius_m)
% Clearance of point r (LVLH meters) from union of spherical and capsule
% obstacles; extraRadius_m = vehicle footprint + additional clearance.
% centers_m [N x 3], radii_m [N x 1] are obstacle sphere radii.
% Optional segmentsA_m / segmentsB_m [K x 3], segmentRadii_m [K x 1]
% describe straight cylindrical/capsule obstacle segments (truss proxies).
% If no per-obstacle radius is supplied, centers are points (zero radius).
assert(numel(r)==3 && extraRadius_m>=0);
r=r(:)';d=inf;
if isfield(hazards,'centers_m') && ~isempty(hazards.centers_m)
    C=hazards.centers_m;
    if isfield(hazards,'radii_m') && ~isempty(hazards.radii_m)
        R=hazards.radii_m(:);
        assert(numel(R)==size(C,1));
    else
        R=zeros(size(C,1),1);
    end

    d=min(d,min(vecnorm(C-r,2,2)-R));
end

if isfield(hazards,'segmentsA_m') && ~isempty(hazards.segmentsA_m)
    A=hazards.segmentsA_m;B=hazards.segmentsB_m;
    assert(all(size(A)==size(B)) && size(A,2)==3);
    if isfield(hazards,'segmentRadii_m')
        R=hazards.segmentRadii_m(:);
    else
        R=zeros(size(A,1),1);
    end

    assert(numel(R)==size(A,1));
    AB=B-A;len2=sum(AB.^2,2);
    t=sum((r-A).*AB,2)./max(len2,eps);
    t=max(0,min(1,t));
    P=A+AB.*t;
    d=min(d,min(vecnorm(P-r,2,2)-R));
end

clearance=d-extraRadius_m;
end
