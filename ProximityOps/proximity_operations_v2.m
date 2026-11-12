
function result = proximity_operations_v2(mu,Re,h,x0,cfg)
% PROXIMITY_OPERATIONS_V2
%
% Configurable FABER close-field rendezvous.
%
% PHASES:
% 1. Homing:     N general HCW transfers
% 2. Holding 1:  stationary point or bounded ellipse
% 3. Closing:    N analytical Vbar/Rbar or general HCW transfers
% 4. Holding 2:  stationary point or bounded ellipse
% 5. Final:      continuously controlled approach
%
% INPUTS:
% mu, Re, h    Orbital parameters [km^3/s^2, km, km]
% x0           Initial RPO state [6x1], [m; m/s]
% cfg          Configuration defined in MA_deltav
%
% OUTPUTS:
% result.homing(i), result.closing(i)
% result.holding1, result.holding2, result.final
% Individual and phase-total DV [m/s] and TOF [min]
%
% Deterministic ideal HCW analysis, NOT Monte Carlo/P95.
% No collision, KOZ, or thruster feasibility certification.

%% CONSTANTS

n = sqrt(mu/(Re+h)^3);
N = cfg.N;

x0 = x0(:);

assert(numel(x0)==6,'Initial state must be 6x1.');

posTol = 1e-4;  % [m] state continuity tolerance
velTol = 1e-7;  % [m/s]

% Current state passed between phases
Xcurrent = x0;

%% ============================================================
% 1. HOMING
% =============================================================

N_homing = size(cfg.homing.targets,2);

assert(size(cfg.homing.targets,1)==6, ...
    'Homing targets must be 6xN.');

assert(numel(cfg.homing.times_s)==N_homing, ...
    'One transfer duration required per homing maneuver.');

homing = struct([]);

DV_homing = 0;
TOF_homing = 0;

for i = 1:N_homing

    Xtarget = cfg.homing.targets(:,i);
    T = cfg.homing.times_s(i);

    assert(T>0,'Homing duration must be positive.');

    sol = HCW_solver( ...
        Xcurrent(1:3),Xcurrent(4:6), ...
        Xtarget(1:3),Xtarget(4:6), ...
        T,n,N);

    homing(i).method = 'HCW';
    homing(i).initialState = Xcurrent;
    homing(i).finalState = Xtarget;

    homing(i).DV_mps = sol.dvTotal;
    homing(i).TOF_min = T/60;
    homing(i).TOF_s = T;
    homing(i).trajectory = sol.trajectory;

    DV_homing = DV_homing + sol.dvTotal;
    TOF_homing = TOF_homing + T/60;

    % State after terminal impulse of HCW transfer
    Xcurrent = Xtarget;

end

%% ============================================================
% 2. HOLDING POINT 1
% =============================================================

h1 = cfg.holding1;

hold1 = rpo_holdingpt( ...
    n,h1.duration_s,h1.mode,h1.param,N);

% Require homing to finish at the holding configuration
checkState(Xcurrent,hold1.holdState, ...
    posTol,velTol,'Holding point 1');

Xcurrent = hold1.holdFinalState;

DV_holding1 = hold1.DV_holding_mps;
TOF_holding1 = hold1.TOF_holding_min;

%% ============================================================
% 3. CLOSING
% =============================================================

N_closing = size(cfg.closing.targets,2);

assert(size(cfg.closing.targets,1)==6, ...
    'Closing targets must be 6xN.');

assert(numel(cfg.closing.methods)==N_closing, ...
    'One closing method required per maneuver.');

closing = struct([]);

DV_closing = 0;
TOF_closing = 0;

for i = 1:N_closing

    Xtarget = cfg.closing.targets(:,i);
    method = char(cfg.closing.methods{i});

    switch lower(method)

        case {'vbar','rbar'}

            % Analytical transfers require:
            % - zero initial and final relative velocities
            % - both positions on V-bar

            assert(norm(Xcurrent(4:6))<velTol && ...
                   norm(Xtarget(4:6))<velTol, ...
                   'Analytical closing requires zero endpoint velocity.');

            sol = analytical_manos( ...
                Xcurrent(1:3), ...
                Xtarget(1:3), ...
                n,method,N,false);

            assert(sol.feasible, ...
                'Analytical closing failed terminal verification.');

            DV = sol.DV_mps;
            TOF = sol.TOF_min;
            T = sol.TOF_s;

            Xend = sol.finalState;

        case 'hcw'

            T = cfg.closing.times_s(i);

            assert(T>0,'Closing duration must be positive.');

            sol = HCW_solver( ...
                Xcurrent(1:3),Xcurrent(4:6), ...
                Xtarget(1:3),Xtarget(4:6), ...
                T,n,N);

            DV = sol.dvTotal;
            TOF = T/60;

            % Desired state following terminal impulse
            Xend = Xtarget;

        otherwise

            error('Unknown closing method: %s',method);

    end

    checkState(Xend,Xtarget, ...
        posTol,velTol,sprintf('Closing %d',i));

    closing(i).method = method;
    closing(i).initialState = Xcurrent;
    closing(i).finalState = Xend;

    closing(i).DV_mps = DV;
    closing(i).TOF_min = TOF;
    closing(i).TOF_s = T;
    closing(i).trajectory = sol.trajectory;

    DV_closing = DV_closing + DV;
    TOF_closing = TOF_closing + TOF;

    Xcurrent = Xend;

end

%% ============================================================
% 4. HOLDING POINT 2
% =============================================================

h2 = cfg.holding2;

hold2 = rpo_holdingpt( ...
    n,h2.duration_s,h2.mode,h2.param,N);

checkState(Xcurrent,hold2.holdState, ...
    posTol,velTol,'Holding point 2');

Xcurrent = hold2.holdFinalState;

DV_holding2 = hold2.DV_holding_mps;
TOF_holding2 = hold2.TOF_holding_min;

%% ============================================================
% 5. FINAL APPROACH
% =============================================================

f = cfg.final;

final = rpo_continuous( ...
    Xcurrent,f.target,n, ...
    f.duration_s,f.method,N);

assert(final.feasible, ...
    'Final approach terminal state not achieved.');

DV_final = final.DV_mps;
TOF_final = final.TOF_min;

Xcurrent = final.finalState;

%% ============================================================
% SAVE INDIVIDUAL STAGES
% =============================================================

result.homing = homing;
result.holding1 = hold1;

result.closing = closing;
result.holding2 = hold2;

result.final = final;

result.N_homing = N_homing;
result.N_closing = N_closing;

result.initialState = x0;
result.finalState = Xcurrent;

%% PHASE DELTA V [m/s]

result.DV_homing = DV_homing;
result.DV_holding1 = DV_holding1;
result.DV_closing = DV_closing;
result.DV_holding2 = DV_holding2;
result.DV_final = DV_final;

result.DV_total = ...
    DV_homing + DV_holding1 + DV_closing + ...
    DV_holding2 + DV_final;

%% PHASE TOF [min]

result.TOF_homing = TOF_homing;
result.TOF_holding1 = TOF_holding1;
result.TOF_closing = TOF_closing;
result.TOF_holding2 = TOF_holding2;
result.TOF_final = TOF_final;

result.TOF_total = ...
    TOF_homing + TOF_holding1 + TOF_closing + ...
    TOF_holding2 + TOF_final;

%% ============================================================
% PHASE SUMMARY TABLE
% =============================================================

Phase = [
    "Homing"
    "Holding 1"
    "Closing"
    "Holding 2"
    "Final approach"
];

Count = [
    N_homing
    1
    N_closing
    1
    1
];

DV_mps = [
    DV_homing
    DV_holding1
    DV_closing
    DV_holding2
    DV_final
];

TOF_min = [
    TOF_homing
    TOF_holding1
    TOF_closing
    TOF_holding2
    TOF_final
];

result.table = table(Phase,Count,DV_mps,TOF_min);

fprintf('\n========================================\n');
fprintf('FABER PROXIMITY OPERATIONS\n');
fprintf('========================================\n');

disp(result.table);

fprintf('Total RPO DV: %.5f m/s\n',result.DV_total);
fprintf('Total RPO TOF: %.2f min\n',result.TOF_total);

end

%% ============================================================
% LOCAL HELPER: CHECK STATE CONTINUITY
% =============================================================

function checkState(Xactual,Xexpected,posTol,velTol,label)

    posError = norm(Xactual(1:3)-Xexpected(1:3));
    velError = norm(Xactual(4:6)-Xexpected(4:6));

    assert(posError<=posTol && velError<=velTol, ...
        ['State mismatch at %s: ', ...
         'position error %.4g m, ', ...
         'velocity error %.4g m/s.'], ...
         label,posError,velError);

end
