%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Ryan Haagenson
%TU Delft
%
% Finding matrix cells intersected by a fracture plane.
%
% Input:
%   CellVertIDs:    2D array (N_cells x N_verts_per_cell)
%                   Each row contains the vertex IDs of all vertices for
%                   for the cell with ID equal to the row index of
%                   el_verts.
%   VertCoords:     2D array (N_verts x 3)
%                   The (x,y,z) coordinates for all the vertices in the
%                   mesh. 
%   FracPlane:      struct
%                   4 fields containing (A,B,C,D) for fracture plane eq.
%
% Output:
%   PIRockCellFlags:    1D array (N_cells x 1)
%                       Flags that represent which matrix cells are
%                       intersected by the fracture plane. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PIRockCellFlags = identify_plane_crossed_volumes(CellVertIDs, ...
                                                          VertCoords, ...
                                                          FracPlane)

tol = 1e-3;

% Evaluate plane equation for every vertex
R = FracPlane.A*VertCoords(:,1) + FracPlane.B*VertCoords(:,2) + FracPlane.C*VertCoords(:,3) + FracPlane.D;

% Gather vertex values per cell
Rcell = R(CellVertIDs);

% Mark padding vertices (0) as NaN
Rcell(CellVertIDs == 0) = NaN;

% Compute extrema per cell
maxR = max(Rcell,[],2,'omitnan');
minR = min(Rcell,[],2,'omitnan');

% Plane intersects if values span both sides or touch plane
PIRockCellFlags = (maxR > tol & minR < -tol) | ...
               (abs(maxR) <= tol) | ...
               (abs(minR) <= tol);

end
