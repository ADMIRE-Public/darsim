%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel
%TU Delft
%
% Monte Carlo alternative to find_distance_frac_rock (EDFM_METHOD
% monte_carlo): average distance between each connected rock cell and the
% fracture plane, estimated from random points sampled inside the cell.
%
% Input / Output: same contract as find_distance_frac_rock.
%   DistOutput{f,1}: connected rock-cell IDs of fracture f
%   DistOutput{f,2}: corresponding average distances
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DistOutput] = find_distance_frac_rock_MC(AFracRockEDFM, ...
                                                   AFracRockPEDFM, ...
                                                   RockMesh, ...
                                                   FracGeometry, ...
                                                   Settings)

Settings = monte_carlo_settings(Settings);

edfm = struct('vol_ind', {AFracRockEDFM.RockID});
pedfm = struct();
if Settings.pedfm
    pedfm.vol_ind = AFracRockPEDFM.RockID;
end

DistOutput = calculate_fracture_average_distance_MC(edfm, pedfm, ...
                                RockMesh, FracGeometry(:), Settings);

end
