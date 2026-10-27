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
%   RockFracInd:    array (N_EDFM x 2)
%                   Column 1 lists rock cell ID of EDFM connections
%                   Column 2 lists frac cell ID of EDFM connections
%   AreaEDFM:       array (N_EDFM x 1)
%                   List of intersection area of EDFM connections
%   Coplanar:       array (N_EDFM x 1)
%                   Flag to denote EDFM connections that are coplanar with
%                   matrix cell face
%   rockmesh:       struct
%                   Contains relevant rock mesh information
%   plane_geo:      struct
%                   Contains information about the fracture plane
%
% Output:
%   int_rock_frac:  array (N_pEDFM x 2)
%                   Column 1 lists rock cell ID of pEDFM connections
%                   Column 2 lists frac cell ID of pEDFM connections
%   pedfm_area:     array (N_pEDFM x 1)
%                   List of area to be used in pEDFM CI values
%   alpha_ind:      array (N_alpha_face x 1)
%                   List of faces in rockmesh selected for alpha adjustment
%   alpha_area:     array (N_alpha_face x 1)
%                   List of area values to be used in alpha calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [int_rock_frac, ...
          pedfm_area, ...
          alpha_ind, ...
          alpha_area] = ...
          find_pdefm_projection_faces(RockFracInd, ...
                                                    AreaEDFM, ...
                                                    Coplanar, ...
                                                    rockmesh, ...
                                                    plane_geo)
% Initialize output variables
int_rock_frac = [ ];
pedfm_area = [ ];
alpha_ind = [ ];
alpha_area = [ ];

% List of neighboring rock cells
% Row index of list = Face ID
CellNeighborPairs = [rockmesh.face.inner.upstream, ...
                     rockmesh.face.inner.downstream; ...
                     rockmesh.face.bound.upstream, ...
                     zeros(size(rockmesh.face.bound.upstream))];

% List of rockmesh face areas
FaceArea = [rockmesh.face.inner.area;...
            rockmesh.face.bound.area];

% List of rockmesh face normal vectors
FaceNormals = [rockmesh.face.inner.normal; rockmesh.face.bound.normal];
FaceNormals = FaceNormals./vecnorm(FaceNormals, 2, 2);

% Normal vector of fracture plane
eq = plane_geo.eq;
FractureNormal = [eq.A, eq.B, eq.C];

% Find the intersection area of EDFM connections, counting only one side
AreaIntersectEDFM = zeros(size(AreaEDFM));
AreaIntersectEDFM(Coplanar == 1) = AreaEDFM(Coplanar == 1);
AreaIntersectEDFM(Coplanar == 0) = AreaEDFM(Coplanar == 0)/2;

% Loop through the list of intersected rock cells
for i = 1:length(RockFracInd(:,1))
    CellId = RockFracInd(i,1);
    FracCellId = RockFracInd(i,2);
    CellCentroid = rockmesh.element.centroid(CellId,:);
   
    % Find all nieghboring cells
    Mask1 = CellNeighborPairs(:,1) == CellId;
    Mask2 = CellNeighborPairs(:,2) == CellId;
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
        LengthScale = (mean(rockmesh.element.volume(CellId)))^(1/3);
        CellCentroid = CellCentroid + (ones(1,3)*1e-3)*LengthScale;
        OnFracture = OnFractureTest(eq, CellCentroid);
    end

    % Loop through neighboring cells
    for NeighborId = Neighbors'
        % Finding if the fracture plate crosses between the matrix cell 
        % and any of its neighbors. If it crosses, the connectivity 
        % between these two matrix cells is affected. 
        
        NeighborCentroid = rockmesh.element.centroid(NeighborId,:);
        DoesIntersect = IntersectionTest(plane_geo, CellCentroid,...
                                         NeighborCentroid);
         
        if DoesIntersect
            % The row index in CellNeighborPairs will provide the ID of the
            % face which separates these two cells.
            FaceId = find((CellNeighborPairs(:,1)==CellId & ...
                           CellNeighborPairs(:,2)==NeighborId) | ...
                          (CellNeighborPairs(:,1)==NeighborId & ...
                           CellNeighborPairs(:,2)==CellId));

            % Create pEDFM connection
            new_pedfm_connection = [NeighborId, FracCellId];

            % Calculate the area for pEDFM connection
            % For simplicity, ignore cases where projection of frac-rock 
            % intersection area doesn't completely overlap with selected 
            % interface. The approach below is taken for CPG. 
            % Note: For coplanar fractures, no pEDFM connection should
            % exist. Thus, pEDFM area is set to zero. However, we still
            % identify the associated projection face for alpha adjustment.
            FaceNormal = FaceNormals(FaceId,:);
            ProjectedArea = (AreaIntersectEDFM(i)) * ...
                            abs(dot(FractureNormal, FaceNormal)) * ...
                            (1.0 - Coplanar(i));

            % Is this the first pEDFM connection?
            if isempty(int_rock_frac)
                int_rock_frac = new_pedfm_connection;
                pedfm_area = ProjectedArea;
            else
                % Check if this pEDFM connection is already listed
                pEDFM_row = find(int_rock_frac(:,1) == ...
                                 new_pedfm_connection(1) & ...
                                 int_rock_frac(:,2) == ...
                                 new_pedfm_connection(2));
                if isempty(pEDFM_row)
                    int_rock_frac = [int_rock_frac; new_pedfm_connection];
                    pedfm_area = [pedfm_area; ProjectedArea];
                else
                    pedfm_area(pEDFM_row) = pedfm_area(pEDFM_row) + ...
                                                ProjectedArea;
                end
            end

            % Check if FaceId is already listed
            if ~any(alpha_ind(:) == FaceId)
                alpha_ind(end+1) = FaceId;
    
                % In this implementation, we would like alpha=1 everywhere. 
                % Therefore, set alpha_area = FaceId
                alpha_area(end+1) = FaceArea(FaceId);
            end
        end
    end
end

alpha_ind = alpha_ind';
alpha_area = alpha_area';

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
        t = dot(plane_geo.coords(1,:) - CellCentroid, plane_geo.normal) / ...
            dot(d,plane_geo.normal);
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
