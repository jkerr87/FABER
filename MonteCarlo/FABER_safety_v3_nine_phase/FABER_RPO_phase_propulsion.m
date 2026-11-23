
function pcfg = FABER_RPO_phase_propulsion(cfg,phaseNumber)
% FABER_RPO_PHASE_PROPULSION
%
% Return propulsion properties for a specified RPO phase.
%
% INPUTS:
% cfg         - Full Monte Carlo configuration
% phaseNumber - Phase number from 1 to 9
%
% OUTPUT:
% pcfg        - Configuration with active propulsion properties

pcfg = cfg;

% Legacy mode: use original common propulsion settings
if ~isfield(cfg,'propulsionMode') || ...
        strcmpi(cfg.propulsionMode,'legacy')
    return;
end

assert(strcmpi(cfg.propulsionMode,'phase'), ...
    'Propulsion mode must be phase or legacy.');

assert(phaseNumber >= 1 && phaseNumber <= 9 && ...
    phaseNumber == round(phaseNumber), ...
    'Phase number must be an integer from 1 to 9.');

% Identify propulsion group
name = cfg.phasePropulsion{phaseNumber};

assert(isfield(cfg.propulsion,name), ...
    'Unknown propulsion group: %s',name);

sys = cfg.propulsion.(name);

% Replace common properties with phase-specific ones
pcfg.thrust_N = sys.thrust_N;
pcfg.nActiveThrusters = sys.nActiveThrusters;
pcfg.Isp_s = sys.Isp_s;
pcfg.thrustEfficiency = sys.thrustEfficiency;

% Optional directional thrust layout
pcfg.thrusterDirections_LVLH = ...
    sys.thrusterDirections_LVLH;

pcfg.thrusterForces_N = ...
    sys.thrusterForces_N;

end
