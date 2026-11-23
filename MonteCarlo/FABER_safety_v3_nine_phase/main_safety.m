
[S_demo,R_demo] = run_FABER_safety_demo();

disp(S_demo.summary);

% Summary of emergency maneuver results
S_demo.summary

% Individual failure cases
S_demo.trials

% Show cases where an emergency CAM was needed
S_demo.trials(S_demo.trials.CAMRequired,:)

% Show unresolved CAM cases
S_demo.trials( ...
    S_demo.trials.CAMRequired & ...
    ~S_demo.trials.CAMFeasible,:)


%% REAL TAXI SAFETY ANALYSIS - PT2 (once demo runs)

% Each column is one physically achievable
% net force vector in LVLH [N].
% camGroups = YOUR_VERIFIED_3xG_FORCE_MATRIX;
% 
% % Captured module and container mass [kg]
% payloadMass_kg = YOUR_PAYLOAD_MASS;
% 
% % First test CAM without recovery
% [S,R] = run_FABER_safety_taxi( ...
%     camGroups,payloadMass_kg,[]);
% 
% % Display results
% disp(S.summary);
% disp(S.trials);