%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Ryan Haagenson
%TU Delft
%
% Finding average distance from fracture cells to fracture intersections. 
% Distances areas are calculated numerically with grid sampling.
%
% Input:
%   fracmesh:       struct
%                   Contains relevant frac mesh information
%   frac_int:       struct
%                   Contains information about the fracture intersections
%   plane_geo:      struct
%                   Contains information about the fracture planes
%   plane_int:      struct
%                   Contains information about the fracture plan
%                   intersections.
%   settings:       struct
%                   Contains information about the simulation settings
%
% Output:
%   frac_int_avg_d: cell
%                   Contains the distance from each frac cell to each frac
%                   cell intersection line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [frac_int_avg_d] = ...
          find_distance_frac_frac_int(fracmesh, ...
                                         frac_int, ...
                                         plane_geo, ...
                                         plane_int, ...
                                         settings)
% Number of fractures
nFrac = max(size(fracmesh));

% Number of points for grid sampling
% Note: In CPG, found with reg grid sampling using 36 (6x6) points
grid_size = 6;

% Initialize output
[frac_avg_d_ind] = create_avg_d_frac_structure(frac_int);
frac_int_avg_d = cell(size(frac_int));

for iFrac = 1:nFrac
    
    [test_points, frac_el] =...
        generate_testpoint_parameters_grid(fracmesh(iFrac), ...
        grid_size);

    if ~all(all(cellfun(@isempty,frac_int)))
        frac_int_avg_d = calculate_avg_frac_d(frac_int_avg_d, test_points,...
        frac_el, frac_avg_d_ind, frac_int, plane_int, iFrac);
    end

end

end