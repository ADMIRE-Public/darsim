%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Ryan Haagenson
%TU Delft
%
% Calculating the average distance between the matrix cell to the fracture
% cell using regular grid sampling. 
%
% Input:
%   edfm:           struct
%                   Fields include: vol_ind, frac_ind, area
%                   Describes the EDFM connections
%   pedfm:          struct
%                   Fields include: vol_ind, frac_ind, area, face_ind,
%                   alpha_ind, alpha_area
%                   Describes the pEDFM connections and alpha-adjusted
%                   projection faces
%   rockmesh:       struct
%                   Contains relevant rock mesh information
%   frac_geo:       struct
%                   Information on the fracture plane
%   settings:       struct
%                   General simulation settings. 
%
% Output:
%   dist:           cell (N_frac x 2)
%                   Field 1: Matrix cell ID (N_connected_matrix_cells x 1)
%                   Field 2: Avg norm dist from matrix cell to frac plane
%                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [DistOutput] = find_distance_frac_rock_int(edfm, ...
                                                       pedfm, ...
                                                       rockmesh, ...
                                                       frac_geo, ...
                                                       settings)

% Number of fractures
NumFrac = size(frac_geo,1);
% Initialize output
DistOutput = cell(NumFrac,2);

% Parameters
n = 7;  % Grid size in x, y, and z directions

% Find unique list of matrix cells connected to each fracture
[ConnectedCells, ConnectedCellsByFrac] = ...
                ConnectedCellFinder(edfm, pedfm, settings.pedfm);
NumConnectedCells = size(ConnectedCells,1);

% Find rock mesh properties
AllCellFaceIDs = rockmesh.element.faces;
AllCellFaceNormals =  [rockmesh.face.inner.normal;...
    rockmesh.face.bound.normal];
AllCellFaceCentroids =  [rockmesh.face.inner.centroid;...
    rockmesh.face.bound.centroid];
AllCellCentroids = rockmesh.element.centroid;

% Generate bounding box for each connected matrix cell
[MinPoints, MaxPoints ] = ...
    generate_boundbox(rockmesh.element.vertices(ConnectedCells,:), ...
    rockmesh.vertex.coord);

% Generate a n x n x n grid of test points within each bounding box
% Filter out points inside bounding boxes, but outside matrix cell polygon
PointsPerBox = n^3;
TestPoints = zeros(PointsPerBox*NumConnectedCells, 3);
TestPointCellIDs = zeros(PointsPerBox*NumConnectedCells, 1);
Filter = zeros(PointsPerBox*NumConnectedCells, 1);

idx = 1;

for i = 1:NumConnectedCells

    CellID = ConnectedCells(i);

    % Get bounds
    minP = MinPoints(i, :);
    maxP = MaxPoints(i, :);
    
    % Generate sampling grid
    x = linspace(minP(1), maxP(1), 7);
    y = linspace(minP(2), maxP(2), 7);
    z = linspace(minP(3), maxP(3), 7);
    
    [X, Y, Z] = ndgrid(x, y, z);
    Points = [X(:), Y(:), Z(:)];
    
    % Grab cell faces and centroid
    CellFaceIDs = AllCellFaceIDs(CellID,:);
    CellFaceNormals = AllCellFaceNormals(CellFaceIDs,:);
    CellFaceCentroids = AllCellFaceCentroids(CellFaceIDs,:);
    CellCentroid = AllCellCentroids(CellID,:);

    % Orient face normals outward
    InnerVec = CellCentroid - CellFaceCentroids;
    dot = sum(CellFaceNormals .* InnerVec, 2);
    flipMask = dot > 0;
    CellFaceNormals(flipMask, :) = -CellFaceNormals(flipMask, :);

    % Filter points outside cell polygon
    InsideFilter = PointsInsidePolyhedron(Points, ...
                                          CellFaceCentroids, ...
                                          CellFaceNormals);

    % Store points
    TestPoints(idx:idx+PointsPerBox-1, :) = Points;    
    TestPointCellIDs(idx:idx+PointsPerBox-1) = repmat(CellID,size(Points,1),1);
    Filter(idx:idx+PointsPerBox-1) = InsideFilter;
    
    idx = idx + PointsPerBox;

end

% Loop over fracs. Find the average distance to each connected matrix cell
for f = 1:NumFrac

    % Fracture info
    FracNorm = frac_geo(f).normal;
    FracPoint = frac_geo(f).center;
    
    % Find connected cells for this fracture
    ConnectedCellsThisFrac = cell2mat(ConnectedCellsByFrac(f));

    TestPointsThisFrac = TestPoints(ismember( ...
        TestPointCellIDs, ConnectedCellsThisFrac),:);
    TestPointsCellIDsThisFrac = TestPointCellIDs(ismember( ...
        TestPointCellIDs, ConnectedCellsThisFrac));
    FilterThisFrac = Filter(ismember( ...
        TestPointCellIDs, ConnectedCellsThisFrac));
    
    % Calculate the distance to each point for each connected cell
    DistancePointToPlane = distance_points_plane( ...
        TestPointsThisFrac,FracNorm,FracPoint);

    % Find the average for each connected cell
    NumPointsByCell = accumarray(TestPointsCellIDsThisFrac, ...
        FilterThisFrac);
    NumPointsByCell = NumPointsByCell(ConnectedCellsThisFrac);
    SumDistanceByCell = accumarray(TestPointsCellIDsThisFrac, ...
        DistancePointToPlane);
    SumDistanceByCell = SumDistanceByCell(ConnectedCellsThisFrac);
    AvgDistByCell = SumDistanceByCell./NumPointsByCell;
    
    % Store results
    DistOutput{f,1} = ConnectedCellsThisFrac;
    DistOutput{f,2} = AvgDistByCell;

end

end

function [vol_list, vol_list_by_frac] = ConnectedCellFinder(edfm, ...
                                                        pedfm, flag)
nfrac = max(size(edfm.vol_ind));
vol_list_by_frac = cell(nfrac,1);
for index = 1:nfrac
    lvol = unique(edfm.vol_ind{index,1});
    if flag
        lvol = union(lvol,unique(pedfm.vol_ind{index,1}));
    end
    vol_list_by_frac{index,1} = lvol;
end
vol_list = unique(cell2mat(vol_list_by_frac));
end

function mask = PointsInsidePolyhedron(points,fc,fn)

d = (points - fc(1,:)) * fn(1,:)';

mask = true(size(points,1),1);

for f = 1:size(fc,1)
    mask = mask & ((points - fc(f,:)) * fn(f,:)' <= 0);
end

end

function dist = distance_points_plane(points,normal,plane_point)
A = normal(:,1);
B = normal(:,2);
C = normal(:,3);
D = -dot(normal,plane_point,2);
num = (A.*points(:,1) + B.*points(:,2) + C.*points(:,3)) + D;
num = abs(num);
den = (A^2 + B^2 + C^2)^.5;
dist = num./den;
end