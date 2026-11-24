function cycles=eps_rainflow(x)
% Three-point stack rainflow count of an open history, including half cycles.
% Output: [count, range, mean, start index, end index].
% See MathWorks rainflow documentation for the ASTM E1049 method description.
x=x(:); assert(all(isfinite(x)),'Nonfinite cycle input.');
if numel(x)<2, cycles=zeros(0,5); return; end
% Remove exact plateaus, preserving the first point of each flat segment.
idx=find([true;diff(x)~=0]); y=x(idx);
if numel(y)<2, cycles=zeros(0,5); return; end
turn=[true;diff(y(1:end-1)).*diff(y(2:end))<0;true];
idx=idx(turn); stack=zeros(numel(idx),1); top=0;
cycles=zeros(numel(idx),5); k=0;
for j=1:numel(idx)
    top=top+1; stack(top)=idx(j);
    while top>=3
        a=stack(top-2); b=stack(top-1); z=stack(top);
        old=abs(x(b)-x(a)); recent=abs(x(z)-x(b));
        if recent<old, break; end
        k=k+1;
        if top==3
            cycles(k,:)=[0.5 old (x(a)+x(b))/2 a b];
            stack(1:2)=stack(2:3); top=2;
        else
            cycles(k,:)=[1 old (x(a)+x(b))/2 a b];
            stack(top-2)=stack(top); top=top-2;
        end
    end
end
for j=1:top-1
    a=stack(j); b=stack(j+1); k=k+1;
    cycles(k,:)=[0.5 abs(x(b)-x(a)) (x(a)+x(b))/2 a b];
end
cycles=cycles(1:k,:); cycles=cycles(cycles(:,2)>0,:);
end
