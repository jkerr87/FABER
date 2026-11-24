function run = FABER_environment_eps(c,h,b,e,l,t,sunFraction)
% PM4 constrained dispatch adapted to external geometry and activity intervals.
% Power is in kW, time in minutes, and energy in kWh.
N = numel(t)-1; V=6; dt = diff(t)/60;
keys = c.names;
solarYield = c.solarFlux_W_m2*c.cellEfficiency_BOL*c.packing* ...
    c.solarEOL*c.pmadEfficiency/1000;
EOL = h.Battery_kWh'*c.batteryEOL;
Emin = EOL*(1-c.maxDoD);
energy = zeros(N+1,V); energy(1,:)=EOL.*c.initialSoC;

own = zeros(N,V); peak=own; pv=own; pvBOL=own; rx=own; tx=own;
served=own; unserved=own; missingLink=own;
batOut=own; batIn=own; curtailed=own; soc=zeros(N+1,V);
requestedRx=own; requestedTx=own;
balanceError=own;
mid=(t(1:end-1)+t(2:end))/2;
assert(isequal(size(sunFraction),[N V]) && all(sunFraction(:)>=0 & sunFraction(:)<=1));
inSun = sunFraction>=1-1e-9;
mode=repmat(["T10","A2","KB","KB","R8","R8"],N,1);
for j=1:height(e)
    v=find(keys==e.Vehicle(j),1);
    mode(mid>=e.Start_min(j)&mid<e.End_min(j),v)=e.Mode(j);
end
support=mode(:,2)=="A2" & (any(ismember(mode(:,3:4),["K2","K7","K3"]),2) | ismember(mode(:,1),["T5","T6"]));
mode(support,2)="A3";
for v=1:V
    for m=unique(mode(:,v))'
        row=b.Vehicle==family(keys(v)) & b.Mode==m;
        assert(nnz(row)==1,'Missing or ambiguous mode %s',m);
        mask=mode(:,v)==m;
        own(mask,v)=b.Mean_W(row)/1000; peak(mask,v)=b.Peak_W(row)/1000;
    end
    pv(:,v)=h.PV_m2(v)*solarYield*h.PointingFactor(v)*h.ShadowFactor(v)* ...
        h.PV_Credited(v)*sunFraction(:,v).*pvAvailability(c,v,mode(:,v));
    pvBOL(:,v)=pv(:,v)/c.solarEOL;
end
for j=1:height(l)
    src=find(keys==l.Sender(j),1); dst=find(keys==l.Receiver(j),1);
    assert((src==2 && any(dst==[3,4])) || (any(src==[3,4]) && dst==src+2),'Invalid link.');
    active=mid>=l.Start_min(j)&mid<l.End_min(j);
    eta=c.interfaceEfficiencyKR; if src==2, eta=c.interfaceEfficiencyAK; end
    requestedRx(active,dst)=requestedRx(active,dst)+l.Received_kW(j);
    requestedTx(active,src)=requestedTx(active,src)+l.Received_kW(j)/eta;
end
assert(all(sum(requestedRx(:,3:4)>0,2)<=c.maxActiveKikiPorts),'Too many active Kiki charging ports.');
for k=1:N
    d=dt(k); wants=requestedRx(k,3:4); needs=requestedRx(k,5:6);
    % Atlante feeds both Kikis, sharing constrained power proportionally.
    available = supplyLimit(c,h,energy(k,2),Emin(2),EOL(2),2,d,pv(k,2),0);
    sendWanted = sum(wants)/c.interfaceEfficiencyAK;
    sentAvailable = max(0,available-own(k,2));
    frac = 1;
    if sendWanted>0, frac=min(1,sentAvailable/sendWanted); end
    deliveredK = wants*frac;
    rx(k,3:4)=deliveredK;
    tx(k,2)=sum(deliveredK)/c.interfaceEfficiencyAK;
    missingLink(k,3:4)=wants-deliveredK;

    % Each Kiki can also feed its own Rango. Own load always has priority.
    for pair=1:2
        id=pair+2;
        available = supplyLimit(c,h,energy(k,id),Emin(id),EOL(id),id,d,...
            pv(k,id),rx(k,id));
        sendAvail = max(0,available-own(k,id));
        rx(k,4+pair)=min(needs(pair),sendAvail*c.interfaceEfficiencyKR);
        tx(k,id)=rx(k,4+pair)/c.interfaceEfficiencyKR;
        missingLink(k,4+pair)=needs(pair)-rx(k,4+pair);
    end

    for v=1:V
        net = pv(k,v)+rx(k,v)-own(k,v)-tx(k,v);
        if net<0
            maxOut=maxDischarge(c,h,EOL(v),Emin(v),energy(k,v),v,d);
            batOut(k,v)=min(-net,maxOut);
            unserved(k,v) = max(0,-net-batOut(k,v));
        else
            maxIn = min(h.ChargeC_per_h(v)*EOL(v), ...
                max(0,(EOL(v)-energy(k,v))/(d*c.batteryChargeEfficiency)));
            batIn(k,v)=min(net,maxIn);
            curtailed(k,v)=max(0,net-batIn(k,v));
        end
        served(k,v)=own(k,v)-unserved(k,v);
        energy(k+1,v)=energy(k,v)+d*(c.batteryChargeEfficiency*batIn(k,v) ...
            -batOut(k,v)/c.batteryDischargeEfficiency);
        energy(k+1,v)=max(Emin(v),min(EOL(v),energy(k+1,v)));
        balanceError(k,v)=pv(k,v)+rx(k,v)+batOut(k,v)- ...
            served(k,v)-tx(k,v)-batIn(k,v)-curtailed(k,v);
    end
end
soc=100*energy./EOL;
battLoss=(1-c.batteryChargeEfficiency)*batIn+(1/c.batteryDischargeEfficiency-1)*batOut;
linkLoss=sum(tx,2)-sum(rx,2);
fleetError=sum((sum(pv-served-curtailed-battLoss,2)-linkLoss).*dt) ...
    -sum(energy(end,:)-energy(1,:));
assert(abs(fleetError)<1e-7,'Fleet energy ledger does not close.');
assert(max(abs(balanceError),[],'all')<1e-7,'Power ledger does not close.');
run = struct('t',t,'dt_h',dt,'mode',mode,'sun',inSun,'E',energy,...
    'soc',soc,'EOL',EOL,'Emin',Emin,'own',own,'peak',peak, ...
    'pv',pv,'pvBOL',pvBOL,...
    'rx',rx,'tx',tx,'requestedRx',requestedRx,'requestedTx',requestedTx,...
    'served',served,'unserved',unserved,'missingLink',missingLink,...
    'chargeLimit',repmat(h.ChargeC_per_h'.*EOL,N,1), ...
    'batOut',batOut,'batIn',batIn,'curtailed',curtailed,'balanceError',balanceError, ...
    'fleetEnergyError_kWh',fleetError);
end

function p = supplyLimit(c,h,E,Emin,EOL,v,d,pv,received)
p = pv+received+maxDischarge(c,h,EOL,Emin,E,v,d);
end
function p = maxDischarge(c,h,EOL,Emin,E,v,d)
p = min(h.DischargeC_per_h(v)*EOL*c.batteryDischargeEfficiency,...
    max(0,(E-Emin)*c.batteryDischargeEfficiency/d));
end
function name = family(vehicle)
if startsWith(vehicle,"Kiki"), name="Kiki";
elseif startsWith(vehicle,"Rango"), name="Rango";
else, name=vehicle;
end
end

%% Verification helpers - demand, reserve, power and energy balance
function f=pvAvailability(c,v,mode)
f=ones(size(mode));
if v==1, f(ismember(mode,["TC","T9","T9b"]))=0; end
if v==2
    f(mode=="AC")=0; f(mode=="A7")=c.atlanteRelocationPVFraction;
end
end

