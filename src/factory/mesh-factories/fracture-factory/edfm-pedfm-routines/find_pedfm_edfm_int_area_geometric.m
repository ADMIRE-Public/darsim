%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Ryan Haagenson
%TU Delft
%
% Finding EDFM and pEDFM connections and intersection areas. 
% Intersection areas are calculated geometrically.
%
% Input:
%   fracmesh:       struct
%                   Contains relevant frac mesh information
%   rockmesh:       struct
%                   Contains relevant rock mesh information
%   plane_geo:      struct
%                   Contains relevant information about fracture properties
%   setting:        struct
%                   Contains relevant information about simulation settings
%
% Output:
%   edfm:           struct
%                   Contains three cell fields:
%                   vol_ind - rock cell index of EDFM connection
%                   frac_ind - fracture cell index of EDFM connection
%                   area - area of EDFM connection
%   pedfm:          struct
%                   Contains six cell field: 
%                   vol_ind - rock cell index of pEDFM connection
%                   frac_ind - fracture cell index of pEDFM connection
%                   area - area of pEDFM connection
%                   face_ind - face index with altered alpha values
%                   organized by each fracture
%                   alpha_ind - face index with altered alpha values in
%                   single list
%                   alpha_area - area value to be used for alpha alteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [edfm, ...
          pedfm ] = ...
          find_pedfm_edfm_int_area_geometric(fracmesh, ...
                                             rockmesh, ...
                                             plane_geo, ...
                                             settings)

% Number of fractures in simulation
nFrac = max(size(fracmesh));

% Initialize output variables
edfm = struct();
RockCellEDFMOutput = cell(nFrac,1);
FracCellEDFMOutput = cell(nFrac,1);
AreaEDFMOutput = cell(nFrac,1);
pedfm = struct();
RockCellPEDFMOutput = cell(nFrac,1);
FracCellPEDFMOutput = cell(nFrac,1);
AreaPEDFMOutput = cell(nFrac,1);
FacePEDFMOutput = cell(nFrac,1);
AlphaFaceOutput = cell(nFrac,1);
AlphaAreaOutput = cell(nFrac,1);

% Loop over each fracture individually
for iFrac = 1:nFrac

    % Fracture Normal
    FracNormal = plane_geo(iFrac).normal;

    % Identify the rock cells intersected by this fracture plane ("PI")
    [PIRockCellFlags] = ...
        identify_plane_crossed_volumes(rockmesh.element.vertices,...
                                       rockmesh.vertex.coord, ...
                                       plane_geo(iFrac).eq);
    
    % Plane intersected rock cell indices
    PIRockCellIDs = find_uint(PIRockCellFlags);

    % Plane intersected rock cell centroids
    PIRockCellCentroids = ...
        rockmesh.element.centroid(PIRockCellFlags,:);

    % Number of plane intersected rock cells
    nPIRockCell = size(PIRockCellIDs,1);

    % Fracture cell centroids
    FracCellCentroids = fracmesh(iFrac).elements.centroid;

    % Number of fracture cells
    nFracCell = size(FracCellCentroids,1);

    % Arrays to track EDFM connections
    RockCellEDFM = [];
    FracCellEDFM = [];
    AreaEDFM = [];
    Coplanar = [];

    % Loop through fracture cells
    for iFracCell = 1:nFracCell

        % Current fracture cell ID
        CurrentFracCellID = iFracCell;

        % Find the fracture cell area
        FracCellArea = fracmesh(iFrac).elements.area(iFracCell);

        % Find indices of the vertices of this frac cell
        FracCellVertIndex = fracmesh(iFrac).elements.vertices(iFracCell,:);
        
        % Global coordinates of fracture cell vertices
        GlobalFracCellVert = ...
            fracmesh(iFrac).vertex.coords(FracCellVertIndex,:);

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
                rockmesh.element.vertices(CurrentPIRockCellID,:);
            CurrentHullVerts = ...
                rockmesh.vertex.coord(CurrentRockCellVertIDs,:);
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
                phi1 = dot(FracNormal,p1) + plane_geo(iFrac).eq.D;
                phi2 = dot(FracNormal,p2) + plane_geo(iFrac).eq.D;
            
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

            % Project intersection polygon vertices from global to local
            LocalPlaneIntPolyVert = (R.' * ...
                                     (GlobalPlaneIntPolyVert - ...
                                      FracCellCentroids(iFracCell,:)).').';

            % Order polygon vertices CCW
            if size(LocalPlaneIntPolyVert,1) < 3
                AreaIntersect = 0;
            else
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
                phi1 = dot(FracNormal,p1) + plane_geo(iFrac).eq.D;
                phi2 = dot(FracNormal,p2) + plane_geo(iFrac).eq.D;
                phi3 = dot(FracNormal,p3) + plane_geo(iFrac).eq.D;

                if abs(phi1) < tol && abs(phi2) < tol && abs(phi3) < tol
                    Coplanar_flag = 1;
                end

            end

            if AreaIntersect > 0

                % Track this EDFM connection
                RockCellEDFM = [RockCellEDFM; CurrentPIRockCellID];
                FracCellEDFM = [FracCellEDFM; CurrentFracCellID];
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

                % Place in output Cell objects
                RockCellEDFMOutput{iFrac,1} = RockCellEDFM;
                FracCellEDFMOutput{iFrac,1} = FracCellEDFM;
                AreaEDFMOutput{iFrac,1} = AreaEDFM;
                
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

    if settings.pedfm
        EDFMconnections = [RockCellEDFM, FracCellEDFM];
        [rock_frac_ind_pdf,int_area_pedfm, rock_face_ind, rock_face_area]= ...
            find_pdefm_projection_faces_centroid_intersect(EDFMconnections, ...
            AreaEDFM, Coplanar, rockmesh, plane_geo(iFrac));
        RockCellPEDFMOutput{iFrac,1} = rock_frac_ind_pdf(:,1); 
        FracCellPEDFMOutput{iFrac,1} = rock_frac_ind_pdf(:,2);
        AreaPEDFMOutput{iFrac,1} = int_area_pedfm;  
        FacePEDFMOutput{iFrac} = rock_face_ind;
        AlphaFaceOutput{iFrac,1} = rock_face_ind;
        AlphaAreaOutput{iFrac,1} = rock_face_area;
    end 

end

% Group alpha values
AlphaAreaOutput = cell2mat(AlphaAreaOutput);
AlphaFaceOutput = cell2mat(AlphaFaceOutput);
[AlphaFaceOutput,~,AlphaGroup] = unique(AlphaFaceOutput);
AlphaAreaOutput = accumarray(AlphaGroup, AlphaAreaOutput);

edfm.vol_ind = RockCellEDFMOutput;
edfm.frac_ind = FracCellEDFMOutput;
edfm.area = AreaEDFMOutput;

if settings.pedfm
    pedfm.vol_ind = RockCellPEDFMOutput;
    pedfm.frac_ind = FracCellPEDFMOutput;
    pedfm.area = AreaPEDFMOutput;
    pedfm.face_ind = FacePEDFMOutput;
    pedfm.alpha_ind = AlphaFaceOutput;
    pedfm.alpha_area = AlphaAreaOutput;
end
