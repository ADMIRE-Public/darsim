%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Ryan Haagenson
%TU Delft
%
% Calculate the average perpendicular distance from each connected rock
% matrix cell to its associated fracture plane using voxel integration.
% Each convex polyhedral cell is sampled at the midpoints of a regular
% grid in its bounding box; points outside the cell are excluded.
%
% Input:
%   AFracRockEDFM:  struct
%                   EDFM connection data (including RockID).
%   AFracRockPEDFM: struct
%                   pEDFM connection data (including RockID).
%   RockMesh:       struct
%                   Contains relevant rock mesh information
%   FracGeometry:    struct
%                   Fracture plane geometry (normal and center).
%   Settings:       struct
%                   Contains relevant simulation settings
%
% Output:
%   DistOutput:  Cell (NumFrac-by-2) For frac f:
%                {f,1}: List of connected rock-cell IDs
%                {f,2}: List of corresponding distances to frac f
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DistOutput] = find_distance_frac_rock(AFracRockEDFM, ...
                                                    AFracRockPEDFM, ...
                                                    RockMesh, ...
                                                    FracGeometry, ...
                                                    Settings)

% Parameters and output
n = 7;  % Number of voxels along each bounding-box direction
NumFrac = numel(FracGeometry);
DistOutput = cell(NumFrac, 2);

% Identify connected rock cells (including pEDFM connections if enabled)
[ConnectedCells, ConnectedCellsByFrac] = ...
    FindConnectedCells(AFracRockEDFM, ...
                       AFracRockPEDFM, ...
                       Settings.pedfm, ...
                       NumFrac);

if isempty(ConnectedCells)
    for iFrac = 1:NumFrac
        DistOutput{iFrac,1} = ConnectedCellsByFrac{iFrac};
        DistOutput{iFrac,2} = zeros(0,1);
    end
    return
end

% Generate midpoint samples and retain only those inside each polyhedron
[VoxelCentroids, VoxelCellIDs] = ...
    GenerateRockVoxelCentroids(RockMesh, ConnectedCells, n);

% Average point-to-plane distances for each fracture
for iFrac = 1:NumFrac
    ConnectedCellsThisFrac = ConnectedCellsByFrac{iFrac};
    DistOutput{iFrac,1} = ConnectedCellsThisFrac;

    if isempty(ConnectedCellsThisFrac)
        DistOutput{iFrac,2} = zeros(0,1);
        continue
    end

    FracNormal = FracGeometry(iFrac).normal;
    FracPoint = FracGeometry(iFrac).center;
    VoxelMask = ismember(VoxelCellIDs, ConnectedCellsThisFrac);

    PointDistances = DistancePointsToPlane( ...
        VoxelCentroids(VoxelMask,:), FracNormal, FracPoint);
    AvgDistByCell = AverageDistancesByCell( ...
        VoxelCellIDs(VoxelMask), PointDistances, ConnectedCellsThisFrac);

    DistOutput{iFrac,2} = AvgDistByCell;
end

end

function [ConnectedCells, ConnectedCellsByFrac] = ...
    FindConnectedCells(edfm, pedfm, UsePEDFM, NumFrac)
    
    % Collect unique rock-cell IDs associated with each fracture.
    ConnectedCellsByFrac = cell(NumFrac,1);
    
    for iFrac = 1:NumFrac
        CellIDs = unique(edfm.RockID{iFrac,1});
        if UsePEDFM
            CellIDs = union(CellIDs, pedfm.RockID{iFrac,1});
        end
        ConnectedCellsByFrac{iFrac} = CellIDs(:);
    end

    ConnectedCells = unique(vertcat(ConnectedCellsByFrac{:}));
    
end

function [VoxelCentroids, VoxelCellIDs] = ...
    GenerateRockVoxelCentroids(rockmesh, ConnectedCells, n)
    
    % Sample bounding boxes and keep only points inside convex rock cells
    NumConnectedCells = numel(ConnectedCells);
    PointsPerBox = n^3;
    
    [MinPoints, MaxPoints] = generate_boundbox( ...
        rockmesh.element.vertices(ConnectedCells,:), ...
        rockmesh.vertex.coord);
    
    AllCellFaceIDs = rockmesh.element.faces;
    AllCellFaceNormals = ...
        [rockmesh.face.inner.normal; rockmesh.face.bound.normal];
    AllCellFaceCentroids = ...
        [rockmesh.face.inner.centroid; rockmesh.face.bound.centroid];
    AllCellCentroids = rockmesh.element.centroid;
    
    % Preallocate maximum capacity; trim rejected points at the end
    VoxelCentroids = zeros(NumConnectedCells*PointsPerBox,3);
    VoxelCellIDs = zeros(NumConnectedCells*PointsPerBox,1);
    NumAccepted = 0;
    
    for iCell = 1:NumConnectedCells
        CellID = ConnectedCells(iCell);
        Points = VoxelGridCentroids(MinPoints(iCell,:), ...
                                    MaxPoints(iCell,:), ...
                                    n);
    
        CellFaceIDs = AllCellFaceIDs(CellID,:);
        CellFaceNormals = AllCellFaceNormals(CellFaceIDs,:);
        CellFaceCentroids = AllCellFaceCentroids(CellFaceIDs,:);
        CellCentroid = AllCellCentroids(CellID,:);
    
        % Orient face normals outward from the convex polyhedron
        InnerVectors = CellCentroid - CellFaceCentroids;
        FlipMask = sum(CellFaceNormals .* InnerVectors,2) > 0;
        CellFaceNormals(FlipMask,:) = -CellFaceNormals(FlipMask,:);
    
        InsideMask = PointsInsidePolyhedron(Points, ...
                                            CellFaceCentroids, ...
                                            CellFaceNormals);
        Points = Points(InsideMask,:);
        NumPoints = size(Points,1);
    
        if NumPoints == 0
            error('find_distance_frac_rock_int:NoInteriorSamples', ...
                'No voxel centroids inside rock cell %d; increase n.', ...
                CellID);
        end
    
        Range = NumAccepted + (1:NumPoints);
        VoxelCentroids(Range,:) = Points;
        VoxelCellIDs(Range) = CellID;
        NumAccepted = NumAccepted + NumPoints;
    end
    
    VoxelCentroids = VoxelCentroids(1:NumAccepted,:);
    VoxelCellIDs = VoxelCellIDs(1:NumAccepted);
end

function [Points] = VoxelGridCentroids(MinPoint, MaxPoint, n)
    % Return centroids of n-by-n-by-n equal-volume bounding-box voxels
    x = MinPoint(1) + ((1:n)-0.5)*(MaxPoint(1)-MinPoint(1))/n;
    y = MinPoint(2) + ((1:n)-0.5)*(MaxPoint(2)-MinPoint(2))/n;
    z = MinPoint(3) + ((1:n)-0.5)*(MaxPoint(3)-MinPoint(3))/n;
    [X,Y,Z] = ndgrid(x,y,z);
    Points = [X(:),Y(:),Z(:)];
end

function [InsideMask] = ...
    PointsInsidePolyhedron(Points, FaceCentroids, FaceNormals)

    % Point is inside convex polyhedron if it is behind every outward face
    InsideMask = true(size(Points,1),1);
    for iFace = 1:size(FaceCentroids,1)
        InsideMask = InsideMask & ...
            ((Points - FaceCentroids(iFace,:))*FaceNormals(iFace,:)' <= 0);
    end
end

function [Distances] = DistancePointsToPlane(Points, Normal, PlanePoint)
    % Perpendicular distance to the infinite plane.
    Normal = Normal(:);
    PlanePoint = PlanePoint(:);
    Distances = abs((Points - PlanePoint.')*Normal)/norm(Normal);
end

function [AvgDistByCell] = ...
    AverageDistancesByCell(VoxelCellIDs, PointDistances, ConnectedCells)
    
    % Arithmetic mean of equally weighted interior voxel samples per cell
    [~, CellPositions] = ismember(VoxelCellIDs, ConnectedCells);
    NumCells = numel(ConnectedCells);
    SumDistances = accumarray(CellPositions, PointDistances, [NumCells,1]);
    NumPointsByCell = accumarray(CellPositions, 1, [NumCells,1]);
    
    if any(NumPointsByCell == 0)
        error('find_distance_frac_rock_int:NoSamplesForConnection', ...
            'A connected rock cell has no interior voxel samples.');
    
    end
    
    AvgDistByCell = SumDistances ./ NumPointsByCell;

end
