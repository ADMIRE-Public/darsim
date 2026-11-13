%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Ryan Haagenson
%TU Delft
%
% Finding pEDFM connections and rockmesh faces for alpha adjustments. 
% Method checks if the fracture plane intersects the line connecting the 
% rock cell centroid with all cell centroids of neighboring rock cells. 
%
% Input:
%   RockIntersection:   struct; with 4 fields
%                       Information regarding EDFM connections frac-frac
%                       FracCellID
%                       RockCellID
%                       AreaEDFM
%                       Coplanar
%   RockMesh:           struct
%                       Contains relevant rock mesh information
%   FracGeometry:       struct
%                       Contains information about the fracture plane
%
% Output:
%   RockFracCellIDs:    array (N_pEDFM x 2)
%                       Column 1 lists rock cell ID of pEDFM connections
%                       Column 2 lists frac cell ID of pEDFM connections
%   AreaPEDFM:          array (N_pEDFM x 1)
%                       List of area to be used in pEDFM CI values
%   AlphaFaceID:        array (N_alpha_face x 1)
%                       List of faces in rockmesh selected for projection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RockFracCellIDs, ...
          AreaPEDFM, ...
          AlphaFaceID] = ...
          find_pdefm_frac_rock(RockIntersections, ...
                               RockMesh, ...
                               FracGeometry)
% Initialize output variables
RockFracCellIDs = [ ];
AreaPEDFM = [ ];
AlphaFaceID = [ ];

% Create array of Rock and Frac Cell ID pairings
RockFracID = [RockIntersections.RockCellID, RockIntersections.FracCellID];

% List of neighboring rock cells
% Row index of list = Face ID
CellNeighborPairs = [RockMesh.face.inner.upstream, ...
                     RockMesh.face.inner.downstream; ...
                     RockMesh.face.bound.upstream, ...
                     zeros(size(RockMesh.face.bound.upstream))];

% List of rockmesh face areas
FaceArea = [RockMesh.face.inner.area;...
            RockMesh.face.bound.area];

% List of rockmesh face normal vectors
FaceNormals = [RockMesh.face.inner.normal; RockMesh.face.bound.normal];
FaceNormals = FaceNormals./vecnorm(FaceNormals, 2, 2);

% Normal vector of fracture plane
eq = FracGeometry.eq;
FractureNormal = [eq.A, eq.B, eq.C];

% Find the intersection area of EDFM connections, counting only one side
AreaIntersectEDFM = zeros(size(RockIntersections.AreaEDFM));
AreaIntersectEDFM(RockIntersections.Coplanar == 1) = ...
    RockIntersections.AreaEDFM(RockIntersections.Coplanar == 1);
AreaIntersectEDFM(RockIntersections.Coplanar == 0) = ...
    RockIntersections.AreaEDFM(RockIntersections.Coplanar == 0)/2;

% Loop through the list of intersected rock cells
for i = 1:length(RockFracID(:,1))
    CellID = RockFracID(i,1);
    FracCellID = RockFracID(i,2);
    CellCentroid = RockMesh.element.centroid(CellID,:);
   
    % Find all nieghboring cells
    Mask1 = CellNeighborPairs(:,1) == CellID;
    Mask2 = CellNeighborPairs(:,2) == CellID;
    Neighbors = nonzeros([CellNeighborPairs(Mask1, 2); ...
                          CellNeighborPairs(Mask2, 1)]);

    % Test if rock cell centroid is on the fracture plane
    OnFracture = OnFractureTest(eq, CellCentroid);
    while OnFracture
        % If the cell centroid is on the fracture then it will intersect
        % with each centroid connection line drawn with neighboring cells. 
        % In this case, adjust the cell centroid slightly for testing
        % the intersections to avoid all cell faces being selected for
        % pEDFM projection.
        LengthScale = (mean(RockMesh.element.volume(CellID)))^(1/3);
        CellCentroid = CellCentroid + (ones(1,3)*1e-3)*LengthScale;
        OnFracture = OnFractureTest(eq, CellCentroid);
    end

    % Loop through neighboring cells
    for NeighborID = Neighbors'
        % Finding if the fracture plate crosses between the matrix cell 
        % and any of its neighbors. If it crosses, the connectivity 
        % between these two matrix cells is affected. 
        
        NeighborCentroid = RockMesh.element.centroid(NeighborID,:);
        DoesIntersect = IntersectionTest(FracGeometry, CellCentroid,...
                                         NeighborCentroid);
         
        if DoesIntersect
            % The row index in CellNeighborPairs will provide the ID of the
            % face which separates these two cells.
            FaceID = find((CellNeighborPairs(:,1)==CellID & ...
                           CellNeighborPairs(:,2)==NeighborID) | ...
                          (CellNeighborPairs(:,1)==NeighborID & ...
                           CellNeighborPairs(:,2)==CellID));

            % Create pEDFM connection
            new_pedfm_connection = [NeighborID, FracCellID];

            % Calculate the area for pEDFM connection
            % For simplicity, ignore cases where projection of frac-rock 
            % intersection area doesn't completely overlap with selected 
            % interface. This approach below is taken for CPG. 
            % Note: For coplanar fractures, no pEDFM connection should
            % exist. Thus, pEDFM area is set to zero. However, we still
            % identify the associated projection face for alpha adjustment.
            FaceNormal = FaceNormals(FaceID,:);
            ProjectedArea = (AreaIntersectEDFM(i)) * ...
                            abs(dot(FractureNormal, FaceNormal)) * ...
                            (1.0 - RockIntersections.Coplanar(i));

            % Is this the first pEDFM connection?
            if isempty(RockFracCellIDs)
                RockFracCellIDs = new_pedfm_connection;
                AreaPEDFM = ProjectedArea;
            else
                % Check if this pEDFM connection is already listed
                pEDFM_row = find(RockFracCellIDs(:,1) == ...
                                 new_pedfm_connection(1) & ...
                                 RockFracCellIDs(:,2) == ...
                                 new_pedfm_connection(2));
                if isempty(pEDFM_row)
                    RockFracCellIDs = ...
                        [RockFracCellIDs; new_pedfm_connection];
                    AreaPEDFM = [AreaPEDFM; ProjectedArea];
                else
                    AreaPEDFM(pEDFM_row) = AreaPEDFM(pEDFM_row) + ...
                                                ProjectedArea;
                end
            end

            % Check if FaceID is already listed
            if ~any(AlphaFaceID(:) == FaceID)
                AlphaFaceID(end+1) = FaceID;
            end
        end
    end
end

AlphaFaceID = AlphaFaceID';

end

function [IntersectBool] = IntersectionTest(plane_geo, CellCentroid, ...
                                            NeighborCentroid)
    eq = plane_geo.eq;
    Val1 = eq.A*CellCentroid(1) + eq.B*CellCentroid(2) + ...
           eq.C*CellCentroid(3) + eq.D;
    Val2 = eq.A*NeighborCentroid(1) + eq.B*NeighborCentroid(2) + ...
           eq.C*NeighborCentroid(3) + eq.D;

    IntersectBool = Val1*Val2 < 0;
    
    % Check that intersection point lies inside the fracture polygon
    if IntersectBool

        % Intersection Point
        d = NeighborCentroid - CellCentroid;
        t = dot(plane_geo.coords(1,:) - ...
            CellCentroid, plane_geo.normal) / dot(d,plane_geo.normal);
        IntPoint = CellCentroid + t*d;

        % Test if inside
        w = IntPoint - plane_geo.coords(1,:);
        u = plane_geo.coords(2,:) - plane_geo.coords(1,:);
        v = plane_geo.coords(3,:) - plane_geo.coords(1,:);
        uu = dot(w, u) / dot(u, u);
        vv = dot(w, v) / dot(v, v);
        IntersectBool = (uu >= 0 && uu <= 1 && vv >= 0 && vv <= 1);

    end

end

function [OnFractureBool] = OnFractureTest(eq, CellCentroid)
    val = eq.A*CellCentroid(1) + eq.B*CellCentroid(2) + ...
          eq.C*CellCentroid(3) + eq.D;

    OnFractureBool = val == 0;
end
