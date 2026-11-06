%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Ryan Haagenson
%TU Delft
%
% Find pairs of fracture mesh cells connected along the intersection line 
% between two rectangular fractures. Also return the length of the line 
% segments with intersection pairs. 
%
% Input:
%   RockMesh:           struct
%                       Contains relevant rock mesh information
%   FracMesh:           struct
%                       Contains relevant frac mesh information
%   iFrac:              int
%                       Index of fracture
%   FracGeometry:       struct
%                       Contains information about the fracture plane
%
% Output:
%   IntersectOutput:    struct; with 4 fields
%                       Information regarding EDFM connections frac-frac
%                       FracCellID
%                       RockCellID
%                       AreaEDFM
%                       Coplanar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function IntersectOutput = find_edfm_frac_rock(RockMesh, ...
                                               FracMesh, ...
                                               iFrac, ...
                                               FracGeometry)

% Initialize output arrays
FracCellID = [];
RockCellID = [];
AreaEDFM = [];
Coplanar = [];

% Fracture Normal
FracNormal = FracGeometry.normal;

% Identify the rock cells intersected by this fracture plane ("PI")
[PIRockCellFlags] = ...
    find_plane_intersected_rock_cells(RockMesh.element.vertices,...
                                      RockMesh.vertex.coord, ...
                                      FracGeometry.eq);

% Plane intersected rock cell indices
PIRockCellIDs = find_uint(PIRockCellFlags);

% Plane intersected rock cell centroids
PIRockCellCentroids = ...
    RockMesh.element.centroid(PIRockCellFlags,:);

% Number of plane intersected rock cells
nPIRockCell = size(PIRockCellIDs,1);

% Fracture cell centroids
FracCellCentroids = FracMesh(iFrac).elements.centroid;

% Number of fracture cells
nFracCell = size(FracCellCentroids,1);

% Loop through fracture cells
for iFracCell = 1:nFracCell

    % Current fracture cell ID
    CurrentFracCellID = iFracCell;

    % Find the fracture cell area
    FracCellArea = FracMesh(iFrac).elements.area(iFracCell);

    % Find indices of the vertices of this frac cell
    FracCellVertIndex = FracMesh(iFrac).elements.vertices(iFracCell,:);
    
    % Global coordinates of fracture cell vertices
    GlobalFracCellVert = ...
        FracMesh(iFrac).vertex.coords(FracCellVertIndex,:);

    % Find local coordinates of fracture cell vertices (rectangular)
    % Local coordinate system origin at fracture cell centroid
    % X and Y axes parallel to fracture cell edges
    
    % Find local X, Y and Z directions
    P0 = GlobalFracCellVert(1,:);            % Reference vert
    PVec = GlobalFracCellVert(2:end,:) - P0; % Vector to other verts
    PDist = vecnorm(PVec,2,2);               % Distance to other verts
    [~,id] = sort(PDist);                    % Sort by distance
    V1 = PVec(id(1),:);                      % First closet vert (x vec)
    V2 = PVec(id(2),:);                      % Second closet vert - Y
    ex = V1 / norm(V1);                      % Local x-axis
    ey = V2 / norm(V2);                      % Local y-axis
    if ~dot(ey,ex) == 0                      % Ensure orthogonality
        ey = ey - dot(ey,ex)*ex;
        ey = ey / norm(ey);
    end
    ez = cross(ex,ey) / norm(cross(ex,ey));  % Local z-axis

    % Find rotation matrix -- rotates vertices from global to local
    R = [ex(:) ey(:) ez(:)];

    % Local coordinates of fracture cell vertices
    LocalFracCellVert = (R.' * ...
                        (GlobalFracCellVert - ...
                         FracCellCentroids(iFracCell,:)).').';

    % Polygon object of the fracture cell in local coordinates
    PgonFrac = polyshape(LocalFracCellVert(:,1), ...
                          LocalFracCellVert(:,2));
    
    % Organize intersected rock cells by distance to frac cell
    PIRockCellFracCellVectors = PIRockCellCentroids - ...
                                FracCellCentroids(iFracCell,:);
    PIRockCellFracCellDistance =...
        vecnorm(PIRockCellFracCellVectors,2,2);
    [~, SortedIndex] = sort(PIRockCellFracCellDistance);
    SortedPIRockCellIDs = PIRockCellIDs(SortedIndex);

    % Loop through sorted rock cell list (prioritize close rock cells).
    % Exit loop once cummulative intersection area = fracture cell area
    CummArea = 0;
    for iPIRockCell = 1:nPIRockCell
        
        % Create convex hull of rock cell iPIRockCell
        CurrentPIRockCellID = SortedPIRockCellIDs(iPIRockCell);
        CurrentRockCellVertIDs = ...
            RockMesh.element.vertices(CurrentPIRockCellID,:);
        CurrentHullVerts = ...
            RockMesh.vertex.coord(CurrentRockCellVertIDs,:);
        CurrentHull = convhulln(CurrentHullVerts);
        
        % Find the edges of CurrentHull
        HullEdges = unique(sort([
            CurrentHull(:,[1 2])
            CurrentHull(:,[2 3])
            CurrentHull(:,[3 1])
        ],2),'rows');
        
        % Initialize global coordinates of intersection polygon verts
        GlobalPlaneIntPolyVert = [];

        % Tolerance for geometric calculation
        tol = 1e-10;

        % Loop through each edge of the current hull
        for iEdge = 1:size(HullEdges,1)

            % End points of edge
            p1 = CurrentHullVerts(HullEdges(iEdge,1),:);
            p2 = CurrentHullVerts(HullEdges(iEdge,2),:);
            
            % Half space test of each end point with fracture plane
            phi1 = dot(FracNormal,p1) + FracGeometry.eq.D;
            phi2 = dot(FracNormal,p2) + FracGeometry.eq.D;
        
            % Identify Intersections
            % Proper intersection
            if phi1 * phi2 < -tol
                t = phi1 / (phi1 - phi2);
                GlobalPlaneIntPolyVert(end+1,:) = p1 + t * (p2 - p1); 
        
            % Endpoint lies on plane
            elseif abs(phi1) < tol && abs(phi2) >= tol
                GlobalPlaneIntPolyVert(end+1,:) = p1; 
            elseif abs(phi2) < tol && abs(phi1) >= tol
                GlobalPlaneIntPolyVert(end+1,:) = p2;
        
            % Entire edge lies on plane
            elseif abs(phi1) < tol && abs(phi2) < tol
                GlobalPlaneIntPolyVert(end+1,:) = p1; 
                GlobalPlaneIntPolyVert(end+1,:) = p2;
            end

        end

        % Remove any duplicates points in GlobalIntersectPolyVert
        GlobalPlaneIntPolyVert = ...
            uniquetol(GlobalPlaneIntPolyVert,tol,'ByRows',true);
        
        if size(GlobalPlaneIntPolyVert,1) < 3
            AreaIntersect = 0;
        else
            % Project intersection polygon verts from global to local
            LocalPlaneIntPolyVert = (R.' * ...
                                 (GlobalPlaneIntPolyVert - ...
                                  FracCellCentroids(iFracCell,:)).').';

            K = convhull(LocalPlaneIntPolyVert(:,1), ...
                         LocalPlaneIntPolyVert(:,2));
            LocalPlaneIntPolyVert = ...
                LocalPlaneIntPolyVert(K(1:end-1),:);

            % Polygon object of the fracture cell in local coordinates
            warning('off', 'MATLAB:polyshape:repairedBySimplify')
            PgonPlaneInt = polyshape(LocalPlaneIntPolyVert(:,1), ...
                                     LocalPlaneIntPolyVert(:,2));
            
            % Find the polygon describing the intersection between the 
            % fracture polygon and the polygon that describes the
            % intersection of the rock cell and the fracture plane in local
            % coordinates. This has equivalent area to the polygon
            % describing the intersection of the fracture and rock cell. 
            PgonIntersect = intersect(PgonPlaneInt, PgonFrac);
            AreaIntersect = area(PgonIntersect);
        end

        % Check if fracture plane is coplanar with any faces
        Coplanar_flag = 0;
        for iFace = 1:size(CurrentHull,1)

            % Face vertices
            p1 = CurrentHullVerts(CurrentHull(iFace,1),:);
            p2 = CurrentHullVerts(CurrentHull(iFace,2),:);
            p3 = CurrentHullVerts(CurrentHull(iFace,3),:);

            % Check if all vertices are on fracture plane
            phi1 = dot(FracNormal,p1) + FracGeometry.eq.D;
            phi2 = dot(FracNormal,p2) + FracGeometry.eq.D;
            phi3 = dot(FracNormal,p3) + FracGeometry.eq.D;

            if abs(phi1) < tol && abs(phi2) < tol && abs(phi3) < tol
                Coplanar_flag = 1;
            end

        end

        if AreaIntersect > 0

            % Track this EDFM connection
            RockCellID = [RockCellID; CurrentPIRockCellID];
            FracCellID = [FracCellID; CurrentFracCellID];
            % Track if any cell faces are coplanar with fracture plane
            Coplanar = [Coplanar; Coplanar_flag];
            % Is the fracture coplanar with a matrix cell face?
            if Coplanar_flag == 1
                % Coplanar -- only count one side of fracture for area
                AreaEDFM = [AreaEDFM; AreaIntersect];
            else
                % Inside cell -- count both sides of fracture for area
                AreaEDFM = [AreaEDFM; 2*AreaIntersect];
            end
            
            % Check if the entire fracture cell is now accounted for
            % Need to account for both sides of fracture
            CummArea = CummArea + AreaIntersect;
            if CummArea/(2*FracCellArea) > 0.9999
                % Skip remaining rock cells in sorted list
                break;
            end

        end

    end

end

% Construct output
IntersectOutput.FracCellID = FracCellID;
IntersectOutput.RockCellID = RockCellID;
IntersectOutput.AreaEDFM = AreaEDFM;
IntersectOutput.Coplanar = Coplanar;

end