%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel
%TU Delft
%
% Copy the EDFM method keywords of the fracture input file into the
% rock-fracture intersection settings. Missing keywords are left out, so
% the analytical method (default) never depends on them.
%
%   EDFM_METHOD                           -> edfm_method
%   FRACTURE_NUM_POINTS                   -> num_area_points
%   FRACTURE_BOUNDING_BOX                 -> num_frac_boundbox
%   FRACTURE_AVERAGE_DISTANCE_REFINEMENT  -> d_avg_refinement
%   FRACTURE_RANDOM_SEED                  -> seed
%   HALF_PRECISION (simulator settings)   -> useHalf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Settings = add_monte_carlo_settings(Settings, FracInput, ...
                                             SimulatorSettings)

Keywords = ["EDFM_METHOD", "edfm_method";
            "FRACTURE_NUM_POINTS", "num_area_points";
            "FRACTURE_BOUNDING_BOX", "num_frac_boundbox";
            "FRACTURE_AVERAGE_DISTANCE_REFINEMENT", "d_avg_refinement";
            "FRACTURE_RANDOM_SEED", "seed"];

for iKey = 1:size(Keywords, 1)
    if isfield(FracInput, Keywords(iKey, 1))
        Settings.(Keywords(iKey, 2)) = FracInput.(Keywords(iKey, 1));
    end
end

if isfield(SimulatorSettings, "half_precision")
    Settings.useHalf = SimulatorSettings.half_precision;
end

end
