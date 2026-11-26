%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel
%TU Delft
%
% Flags the rock cells having a face that lies on the fracture plane. For
% a convex cell, the fracture part inside such a cell is that face, i.e.
% the EDFM connection is coplanar.
%
% Input:
%   RockID:     rock cell IDs
%   RockMesh:   rock mesh
%   Eq:         fracture plane equation (fields A, B, C, D)
% Output:
%   Coplanar:   logical, one entry per RockID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Coplanar = coplanar_rock_cells(RockID, RockMesh, Eq)

Normals = [RockMesh.face.inner.normal; RockMesh.face.bound.normal];
Normals = Normals ./ vecnorm(Normals, 2, 2);
FaceCentroids = [RockMesh.face.inner.centroid; RockMesh.face.bound.centroid];

Faces = RockMesh.element.faces(RockID, :);
Valid = Faces ~= 0;
Faces(~Valid) = 1;

PlaneNormal = [Eq.A Eq.B Eq.C];
Offset = Eq.D / norm(PlaneNormal);
PlaneNormal = PlaneNormal / norm(PlaneNormal);

Parallel = reshape(abs(Normals(Faces, :) * PlaneNormal'), size(Faces)) ...
           > 1 - 1e-8;
Distance = reshape(abs(FaceCentroids(Faces, :) * PlaneNormal' + Offset), ...
                   size(Faces));
Tol = 1e-6 * RockMesh.element.volume(RockID(:)).^(1/3);

Coplanar = any(Valid & Parallel & Distance <= Tol, 2);

end
