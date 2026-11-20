%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel;
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xyz_testpoints, frac_ind] = ...
           generate_testpoint_parameters(fracmesh,plane_normal, n, useHalf)
%% creates the elements described by the adj and coords and populates it 
% with n x n equally spaced points around the bounding box of each element.
% frac_ind indicates to which frac element the points in xyz_testpoints
% belong.  bound_box_area contains the area of the bound_box of each
% element.
adj = fracmesh.elements.vertices;
coords = fracmesh.vertex.coords;
[xy_coords] = proj_coplanar_coords_on_xy(coords(adj(1,1),:), ...
                                                plane_normal,coords,1);    
zmean = mean(xy_coords(:,3));
[min_vec, max_vec] = generate_boundbox(adj,xy_coords, useHalf);

frac_ind = uint32(repelem([1:size(adj,1)]',n*n,1));

xy_points(:,1) = min_vec(frac_ind,1) + ...
                    (max_vec(frac_ind,1) - min_vec(frac_ind,1)).* ...
                    rand(size(frac_ind,1),1);
xy_points(:,2) = min_vec(frac_ind,2) + ...
                    (max_vec(frac_ind,2) - min_vec(frac_ind,2)).* ...
                    rand(size(frac_ind,1),1);
%% if available, use half precision
if useHalf
    xy_points = half(xy_points);
end

%% check points inside
[ref] = vectorize_inpolygon(xy_points, frac_ind,xy_coords,...
  fracmesh.elements.faces, adj, [fracmesh.faces.inner.vertices; ...
                                        fracmesh.faces.bound.vertices]);
%% removing points not contained on the fracture mesh element
xy_points = xy_points(ref,:);
frac_ind = frac_ind(ref,:);
xy_testpoints = [xy_points zmean*ones(size(xy_points,1),1)];
[xyz_testpoints] = proj_coplanar_coords_on_xy(coords(adj(1,1),:), ...
                                             plane_normal,xy_testpoints,0);


end

