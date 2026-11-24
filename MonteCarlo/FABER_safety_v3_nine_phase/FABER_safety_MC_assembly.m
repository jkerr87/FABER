function R = FABER_safety_MC_assembly(A,mu,Re,h,modulePoints,cfg)
% FABER_SAFETY_MC_ASSEMBLY Selected weighted v4 trajectory failure screening.
% Includes two HCW burns and powered-phase cutoff checkpoints.
% Vehicle mass differs for Kiki/Rango and loaded/unloaded legs.
% Geometry: historical module interface spheres + optional user-supplied
% capsule truss-member geometry. NOT a CAD swept-volume collision checker.
if nargin<6 || isempty(cfg),cfg=struct();end
cfg=FABER_safety_defaults(cfg,mu,Re,h,'assembly');
assert(A.feasible,'All selected assembly trajectories must be feasible');
assert(A.Nmodules==size(modulePoints,1));
assert(isfield(cfg,'hostPoint_m') && numel(cfg.hostPoint_m)==3);
assert(isfield(cfg,'massKiki_kg') && isfield(cfg,'massRango_kg'));
assert(isfield(cfg,'moduleMass_kg') && cfg.moduleMass_kg>=0);
proto=struct('xPre',[],'dvNom',[],'phase','','vehicle','', ...
    'index',0,'hazardsCenters_m',[],'mass_kg',nan, ...
    'hazardsSegmentsA_m',[],'hazardsSegmentsB_m',[], ...
    'hazardsSegmentRadii_m',[],'contactUnresolved',false);
events=repmat(proto,0,1);
for i=1:A.Nmodules
    rec=A.bestRecords{i};
    if isempty(rec),continue;end
    if mod(i,2),veh='Kiki';mass=cfg.massKiki_kg;
    else,veh='Rango';mass=cfg.massRango_kg;end
    % Points represent ONLY already-placed interfaces plus the worksite.
    % For structural clearance, supply CAD-derived segment capsules.
    centers=[cfg.hostPoint_m(:)';modulePoints(1:i,:)];
    segA=[];segB=[];segR=[];
    if isfield(cfg.hazards,'segmentsA_m') && ...
            ~isempty(cfg.hazards.segmentsA_m)
        use=true(size(cfg.hazards.segmentsA_m,1),1);
        if isfield(cfg,'trussInstallModule_i') && ...
                ~isempty(cfg.trussInstallModule_i)
            use=cfg.trussInstallModule_i(:)<=i;
        end

        segA=cfg.hazards.segmentsA_m(use,:);
        segB=cfg.hazards.segmentsB_m(use,:);
        segR=cfg.hazards.segmentRadii_m(use);
    end

    legs={rec.outbound,rec.returnTrip};
    for leg=1:2
        a=legs{leg};
        if leg==1
            legName='outbound loaded';m=mass+cfg.moduleMass_kg;
        else
            legName='return unloaded';m=mass;
        end

        Xhcw=a.transferPath;
        if size(Xhcw,1)==6
            events(end+1)=makeEvent([a.exitGate_m;zeros(3,1)], ...
                Xhcw(4:6,1),['HCW departure ' legName],veh,i,centers,m,segA,segB,segR); %#ok<AGROW>
            events(end+1)=makeEvent(Xhcw(:,end), ...
                -Xhcw(4:6,end),['HCW arrival ' legName],veh,i,centers,m,segA,segB,segR); %#ok<AGROW>
        end
        paths={a.exitPath,a.transferPath,a.finalPath};
        durations=[a.exitTOF_min,a.transferTOF_min,a.finalTOF_min]*60;
        names={'controlled exit','HCW coast','controlled arrival'};
        for j=1:3
            X=paths{j};
            if size(X,1)~=6 || size(X,2)<3,continue;end
            N=max(2,ceil(durations(j)/cfg.maxSampleStep_s)+1);
            ix=unique(max(2,min(size(X,2)-1, ...
                round(linspace(2,size(X,2)-1,N)))));
            for k=ix
                x=X(:,k);
                % NEVER silently discard contact-region checkpoints.
                % Flag unresolved contact geometry, keep event in coverage.
                nearContact=FABER_clearance_m(x(1:3), ...
                    struct('centers_m',centers), ...
                    cfg.screen.collisionRadius_m)<=0;
                ev=makeEvent(x,nan(3,1), ...
                    [names{j} ' ' legName],veh,i,centers,m,segA,segB,segR);
                ev.contactUnresolved=nearContact;
                events(end+1)=ev; %#ok<AGROW>
            end
        end
    end
end
assert(~isempty(events),'No assembly events available for failure screening');
R=FABER_safety_MC_core(events,cfg);
R.label='Kiki/Rango selected assembly trajectories: conditional failure screening';
R.warning=['Simplified sphere/capsule obstacles and ideal nominal paths. ', ...
    'Contact-region checkpoints retained but marked UNRESOLVED, ', ...
    'no verified docking authorization, vehicle swept volumes, ', ...
    'or proven 24 h holding state.'];
end

function e=makeEvent(x,dv,p,v,i,centers,m,segA,segB,segR)
e=struct('xPre',x,'dvNom',dv,'phase',p,'vehicle',v, ...
    'index',i,'hazardsCenters_m',centers,'mass_kg',m, ...
    'hazardsSegmentsA_m',segA,'hazardsSegmentsB_m',segB, ...
    'hazardsSegmentRadii_m',segR,'contactUnresolved',false);
end
