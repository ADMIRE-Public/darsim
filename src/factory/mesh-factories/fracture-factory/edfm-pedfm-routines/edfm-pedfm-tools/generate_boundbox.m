%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Ryan Haagenson; Artur Castiel
%TU Delft
%
% Generating axis-aligned bounding boxes for polyhedra defined by a set
% of vertices. The bounding boxes are determined by the minimum and
% maximum Cartesian coordinates of the vertices of each polyhedron.
%
% Input:
%   adj:            array (N x M)
%                   Vertex indices defining each polyhedron. Each row
%                   corresponds to a polyhedron, with zero entries
%                   indicating unused vertex positions
%   coords:         array (N_vert x 3)
%                   Cartesian coordinates of all vertices in the mesh
%
% Output:
%   minVec:         array (N x 3)
%                   Minimum X, Y, and Z coordinates of each bounding box
%   maxVec:         array (N x 3)
%                   Maximum X, Y, and Z coordinates of each bounding box
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [minVec, maxVec] = generate_boundbox(adj, coords)

% Initialize bounding box coordinates using the first vertex
minVec = coords(adj(:,1),:);
maxVec = coords(adj(:,1),:);

% Number of vertices defining each polyhedron
nVert = size(adj,2);

% Loop through each vertex position
for iVert = 1:nVert

    % Identify polyhedra with a valid vertex at this position
    ref = adj(:,iVert) ~= 0;

    % Current minimum and maximum bounding box coordinates
    LocalMin = minVec(ref,:);
    LocalMax = maxVec(ref,:);

    % Coordinates of current vertices
    CurrentVertCoords = coords(adj(ref,iVert),:);

    % Identify coordinates smaller than the current minimum
    MinRef = CurrentVertCoords - LocalMin <= 0;
    LocalMin(MinRef) = CurrentVertCoords(MinRef);

    % Identify coordinates larger than the current maximum
    MaxRef = CurrentVertCoords - LocalMax >= 0;
    LocalMax(MaxRef) = CurrentVertCoords(MaxRef);

    % Update bounding box coordinates
    minVec(ref,:) = LocalMin;
    maxVec(ref,:) = LocalMax;

end

end
