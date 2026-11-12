function [moduleLarge,moduleSmall] = FABER_make_config_geometry(envelope,N_largeModules,N_smallElements)
% Deterministic PLACEHOLDER attachment points, NOT real client CAD.
% Input envelope = [length; width; height] [m].
% Points are centred on the nominal host origin; clearance not validated.
assert(numel(envelope)==3 && all(envelope>0));
assert(N_largeModules>=0 && N_smallElements>=0);
L=envelope(1); W=envelope(2); H=envelope(3);
% Large modules: roughly uniform slots in the horizontal mid-plane.
if N_largeModules==0
    moduleLarge=zeros(0,3);
else
    nc=ceil(sqrt(N_largeModules*L/W));
    nr=ceil(N_largeModules/nc);
    [XX,YY]=meshgrid(linspace(-.40*L,.40*L,nc), ...
                     linspace(-.40*W,.40*W,nr));
    pts=[XX(:),YY(:),zeros(numel(XX),1)];
    moduleLarge=pts(1:N_largeModules,:);
end
% Small elements: distribute sample sites across faces of the envelope.
% The modular-arithmetic sequence keeps positions repeatable across runs.
moduleSmall=zeros(N_smallElements,3);
for k=1:N_smallElements
    face=mod(k-1,6)+1;
    s=mod(k*0.618033988749895,1); % fraction in [0,1)
    t=mod(k*0.414213562373095,1);
    xx=(s-0.5)*0.90*L;
    yy=(t-0.5)*0.90*W;
    zz=(s-0.5)*0.90*H;
    switch face
      case 1, moduleSmall(k,:)=[xx,yy, H/2];
      case 2, moduleSmall(k,:)=[xx,yy,-H/2];
      case 3, moduleSmall(k,:)=[ L/2,yy,zz];
      case 4, moduleSmall(k,:)=[-L/2,yy,zz];
      case 5, moduleSmall(k,:)=[xx, W/2,zz];
      case 6, moduleSmall(k,:)=[xx,-W/2,zz];
    end
end
end
