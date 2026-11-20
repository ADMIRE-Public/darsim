%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel
%TU Delft
%
% Normalize the Monte Carlo (p)EDFM settings coming from the fracture
% input file (values may arrive as text or cells) and fill defaults.
%
%   num_area_points   (FRACTURE_NUM_POINTS)  n -> n^2 points per frac cell
%   num_frac_boundbox (FRACTURE_BOUNDING_BOX) n -> n^2 sub-boxes per
%                     fracture, processed one at a time to bound memory
%   d_avg_refinement  (FRACTURE_AVERAGE_DISTANCE_REFINEMENT) n -> n^3
%                     points per rock cell for the average distance
%   useHalf           (HALF_PRECISION) half-precision test points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Settings = monte_carlo_settings(Settings)

Defaults = struct("num_area_points", 100, ...
                  "num_frac_boundbox", 1, ...
                  "d_avg_refinement", 35, ...
                  "useHalf", false);

for Name = string(fieldnames(Defaults))'
    if ~isfield(Settings, Name) || isempty(Settings.(Name))
        Settings.(Name) = Defaults.(Name);
    end
    Settings.(Name) = to_scalar(Settings.(Name));
end
Settings.useHalf = logical(Settings.useHalf);

end

function Value = to_scalar(Value)
if iscell(Value)
    Value = Value{1};
end
if ischar(Value) || isstring(Value)
    Text = lower(strtrim(string(Value)));
    if Text == "true"
        Value = true;
    elseif Text == "false"
        Value = false;
    else
        Value = str2double(Text);
    end
end
end
