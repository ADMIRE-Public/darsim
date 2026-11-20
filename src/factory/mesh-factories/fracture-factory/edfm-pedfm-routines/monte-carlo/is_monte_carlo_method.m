%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel
%TU Delft
%
% True when the fracture input selects EDFM_METHOD monte_carlo. The
% analytical method ("analytical", or legacy "geometric") is the default.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Flag = is_monte_carlo_method(Settings)

Flag = false;
if ~isfield(Settings, "edfm_method") || isempty(Settings.edfm_method)
    return
end

Method = Settings.edfm_method;
if iscell(Method)
    Method = Method{1};
end
Method = lower(strtrim(string(Method)));

switch Method
    case "monte_carlo"
        Flag = true;
    case {"analytical", "geometric"}
        Flag = false;
    otherwise
        error("DARSim Error: EDFM_METHOD '%s' not recognized. " + ...
              "Available options: analytical, monte_carlo", Method);
end

end
