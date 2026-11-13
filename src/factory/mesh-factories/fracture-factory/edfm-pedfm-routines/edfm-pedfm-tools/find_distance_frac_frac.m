%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Ryan Haagenson
%TU Delft
%
% Calculate the average perpendicular distance from each connected fracture
% cell to its associated fracture-fracture intersection line using voxel
% integration. Each rectangular fracture cell is sampled at the midpoints
% of a regular grid in its plane.
%
% Input:
%   LFracFracEDFM:      truct
%                       Contains information about the frac-frac EDFM
%                       intersections
%   FracMesh:           struct
%                       Contains relevant frac mesh information
%   GeoIntersection:    struct
%                       Contains information about the fracture plan
%                       intersections.
%   Settings:           struct
%                       Contains relevant simulation settings
%
% Output:
%   DistOutput: Cell (NumFrac-by-1); Each row represents a frac f1
%               Each entry is Cell (NumFrac-by-2); Each row is for frac f2
%               For f1 and f2:
%               {f1,1}{f2,1}: List of frac-cell IDs from f1 connected to f2
%               {f1,1}{f2,2}: List of distances from f1-cells to f2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DistOutput] = find_distance_frac_frac(LFracFracEDFM, ...
                                                LFracFracPEDFM, ...
                                                FracMesh, ...
                                                GeoIntersection, ...
                                                Settings)

% Parameters and output
n = 7;  % Number of midpoint samples along each fracture-cell edge
NumFrac = max(size(FracMesh));
DistOutput = cell(NumFrac,1);

% Generate midpoint samples for each fracture and evaluate its connections
for iFrac1 = 1:NumFrac

    % Initialize the cell for this fracture
    DistOutput{iFrac1,1} = cell(NumFrac,2);

    % Find which fractures this fracture intersects with
    IntersectsWithThisFrac = find(GeoIntersection.intersectWith(:,iFrac1));
    nIntersectingFracs = numel(IntersectsWithThisFrac);
    
    % Check if there are no intersecting fractures
    if nIntersectingFracs > 0

        % Loop through intersecting fractures
        for iFrac2 = IntersectsWithThisFrac'

            % Find all iFrac1 cells connected to iFrac2
            if iFrac1 < iFrac2
                FracConnectionInfo = LFracFracEDFM.Frac1ID{iFrac1,iFrac2};
                ConnectedCellsFrac1 = unique(FracConnectionInfo(:));
                if Settings.pedfm
                    FracConnectionInfo = ...
                        LFracFracPEDFM.Frac1ID{iFrac1,iFrac2};
                    if ~isempty(FracConnectionInfo)
                        ConnectedCellsFrac1 = ...
                        unique([ConnectedCellsFrac1; ...
                                FracConnectionInfo(:)]);
                    end
                end
            else
                FracConnectionInfo = LFracFracEDFM.Frac2ID{iFrac2,iFrac1};
                ConnectedCellsFrac1 = unique(FracConnectionInfo(:));
                if Settings.pedfm
                    FracConnectionInfo = ...
                        LFracFracPEDFM.Frac2ID{iFrac2,iFrac1};
                    if ~isempty(FracConnectionInfo)
                        ConnectedCellsFrac1 = ...
                        unique([ConnectedCellsFrac1; ...
                                FracConnectionInfo(:)]);
                    end
                end
            end
            
            
            % Find the intersection line
            IntLineEndPoints = cell2mat(...
                GeoIntersection.get_boundary_line_intersection(iFrac1, ...
                                                               iFrac2));
            
            % Generate voxels for all cells from Frac1
            [VoxelCentroids, VoxelCellIDs] = ...
                GenerateFractureVoxelCentroids(FracMesh(iFrac1), n);

            % Remove voxels for unconnected cells from Frac1
            VoxelCentroids = ...
                VoxelCentroids(...
                ismember(VoxelCellIDs,ConnectedCellsFrac1),:);
            VoxelCellIDs = ...
                VoxelCellIDs(...
                ismember(VoxelCellIDs,ConnectedCellsFrac1));

            % Calculate distance from each voxel centroid
            VoxelDistance = ...
                DistanceVoxelToLine(VoxelCentroids, IntLineEndPoints);

            % Calculate the average distance by frac cell
            CellDistancesFrac1 = ...
                accumarray(VoxelCellIDs,VoxelDistance,[],@mean);
            CellDistancesFrac1 = CellDistancesFrac1(ConnectedCellsFrac1);

            DistOutput{iFrac1,1}{iFrac2,1} = ConnectedCellsFrac1;
            DistOutput{iFrac1,1}{iFrac2,2} = CellDistancesFrac1;

        end
    end

end

end

function [VoxelCentroids, VoxelCellIDs] = ...
    GenerateFractureVoxelCentroids(Fracture, n)
    
    % Generate n-by-n midpoint samples on each rectangular fracture cell
    NumCells = size(Fracture.elements.area,1);
    PointsPerCell = n^2;
    VoxelCentroids = zeros(NumCells*PointsPerCell,3);
    VoxelCellIDs = zeros(NumCells*PointsPerCell,1);
    
    for iCell = 1:NumCells
        Vertices = ...
            Fracture.vertex.coords(Fracture.elements.vertices(iCell,:),:);
        Range = (iCell-1)*PointsPerCell + (1:PointsPerCell);
        VoxelCentroids(Range,:) = VoxelGridCentroids(Vertices,n);
        VoxelCellIDs(Range) = iCell;
    end
end

function [Centroids] = VoxelGridCentroids(Vertices,n)
    % Affine midpoint grid. Vertices 2 and 4 are adjacent to vertex 1
    PointA = Vertices(1,:);
    EdgeU = Vertices(2,:) - PointA;
    EdgeV = Vertices(4,:) - PointA;
    [S,T] = ndgrid(((1:n)-0.5)/n, ((1:n)-0.5)/n);
    Centroids = PointA + S(:)*EdgeU + T(:)*EdgeV;
end

function [Distances] = DistanceVoxelToLine(Points, Line)
    % Perpendicular distance to the infinite 3D line through two endpoints
    LineDirection = Line(2,:) - Line(1,:);
    LineLength = norm(LineDirection);
    
    if LineLength == 0
        error('find_distance_frac_frac_int:DegenerateLine', ...
            'Fracture intersection line has zero length.');
    end

    Distances = vecnorm(cross(Points - Line(1,:), ...
        repmat(LineDirection,size(Points,1),1),2),2,2)/LineLength;
end
