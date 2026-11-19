function events = FABER_capture_continuous_events(X,t_s,phase,vehicle,legIndex,step_s)
% Capture thrust-loss checkpoints on a nominal controlled 6-state profile.
% X is 6xN, time t_s [1xN] seconds, in the SAME LVLH frame as the HCW model.
% The off-nominal screen treats thrust as fully removed AFTER each state.
% It does not reconstruct the nominal guidance law or partial-burn errors.
if nargin<6 || isempty(step_s),step_s=60;end
assert(size(X,1)==6 && size(X,2)==numel(t_s));
assert(step_s>0 && all(diff(t_s)>0));
points=unique([t_s(1):step_s:t_s(end),t_s(end)]);
events=struct('xPre',{},'dvNom',{},'phase',{},'vehicle',{},'index',{});
for j=1:numel(points)
    x=interp1(t_s(:),X',points(j),'linear')';
    events(end+1)=struct('xPre',x,'dvNom',nan(3,1), ...
        'phase',[char(phase) ' thrust cutoff'], ...
        'vehicle',char(vehicle),'index',legIndex); %#ok<AGROW>
end
end
