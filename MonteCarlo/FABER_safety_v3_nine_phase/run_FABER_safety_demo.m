function [S,R] = run_FABER_safety_demo()
% PURE DIAGNOSTIC ONLY. Hypothetical omnidirectional 22-N groups and zero
% payload mass; NOT actual Taxi propulsion or physically verified flight data.
warning(['DEMO ONLY: synthetic +/-LVLH force groups, ZERO captured ', ...
    'module mass. Do NOT use output to justify FABER CAM DV budget.']);
F=22*[eye(3),-eye(3)];
[S,R]=run_FABER_safety_taxi(F,0,[]);
end
