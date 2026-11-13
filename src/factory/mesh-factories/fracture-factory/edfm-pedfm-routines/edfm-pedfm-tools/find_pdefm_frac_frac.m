%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Ryan Haagenson
%TU Delft
%
% Identifies the fracture cells in iFrac that connect to cells in another
% fracture through a pEDFM connection, along with the associated projected 
% lengths of those pEDFM connections. The pEDFM connections are found based 
% on the EDFM connections between iFrac and another fracture, provided 
% through FractureIntersections. The internal edges/faces in the iFrac mesh
% selected for alpha projection are also found. 
%
% Input:
%   FractureIntersections:  struct; with 7 fields
%                           Info regarding EDFM connections frac-frac
%                           FracID1
%                           FracID2
%                           CellID1
%                           CellID2
%                           Length
%                           OnEdge1
%                           OnEdge2
%   FracMesh:               struct
%                           Contains relevant frac mesh information
%   iFrac:                  int
%                           Index of fracture of interest
%   IntLineEndPoints:       array (2-by-3)
%                           Each row gives (x,y,z) coords of end point of
%                           the intersection line between iFrac and iFrac2
%
% Output:
%   FracFracCellIDs:    array (N_pEDFM x 2)
%                       Column 1 lists iFrac cell ID of pEDFM connections
%                       Column 2 lists iFrac2 cell ID of pEDFM connections
%   LengthPEDFM:        array (N_pEDFM x 1)
%                       List of length to be used in pEDFM CI values
%   AlphaFaceID:        array (N_alpha_face x 1)
%                       List of faces in iFrac selected for projection    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [FracFracCellIDs, LengthPEDFM, AlphaFaceID] = ...
    find_pdefm_frac_frac(FractureIntersections, ...
                         FracMesh, ...
                         iFrac, ...
                         IntLineEndPoints)

% Initialize outputs
FracFracCellIDs = [];
LengthPEDFM = [];
AlphaFaceID = [];

% Find intersection line direction
P0 = IntLineEndPoints(1,:);
P1 = IntLineEndPoints(2,:);
D = P1 - P0;
L = norm(D);

if L == 0
    error('Intersection line has zero length.');
end

LineDir = D/L;

% Grab cell IDs & length for EDFM connections of each fracture
EDFMCellID1 = FractureIntersections.CellID1(:);
EDFMCellID2 = FractureIntersections.CellID2(:);
EDFMLength = FractureIntersections.LengthEDFM(:);

% Find the mesh information for iFrac
VertCoords = FracMesh(iFrac).vertex.coords;
ElemVertIDs = FracMesh(iFrac).elements.vertices;
ElemFaceIDs = FracMesh(iFrac).elements.faces;
ElemCentroids = FracMesh(iFrac).elements.centroid;

% Build neighboring-cell connectivity
[Neighbors, NeighborFaces] = ...
    build_cell_neighbors(ElemFaceIDs);

% Process unique EDFM-connected cells for iFrac
uniqueEDFMCells = unique(EDFMCellID1,'stable');

% Loop through each EDFM connected cell within iFrac
for iEDFMCell = 1:numel(uniqueEDFMCells)
    
    % ID and centroid of this cell
    CellID = uniqueEDFMCells(iEDFMCell);
    CellCentroid = ElemCentroids(CellID,:);

    % Examine each neighboring cell
    for NeighborID = 1:numel(Neighbors{CellID})
        % Finding if the fracture intersection line crosses between the
        % frac cell  and any of its neighbors. If it crosses, the 
        % connectivity between these two frac cells is affected. 
        
        % ID and centroid of neighbor cell
        NeighborCellID = Neighbors{CellID}(NeighborID);
        NeighborCentroid = ElemCentroids(NeighborCellID,:);

        % ID of face between two frac cells
        FaceID = NeighborFaces{CellID}(NeighborID);

        % Does fracture intersection line intersect centroid connector?
        [intersects] = intersect_segments_3D(P0, ...
                                             P1, ...
                                             CellCentroid, ...
                                             NeighborCentroid);
        
        % Check if centroid connector is not intersected
        % If not, the face between them is unaffected
        if ~intersects
            continue
        end

        % Record this projection face only once
        AlreadyRecorded = any(AlphaFaceID == FaceID);
        if ~AlreadyRecorded
            AlphaFaceID = [AlphaFaceID; FaceID];
        end

        % Identify shared face vertices
        CommonVertices = intersect(ElemVertIDs(CellID,:), ...
                                   ElemVertIDs(NeighborCellID,:));
        VA = VertCoords(CommonVertices(1),:);
        VB = VertCoords(CommonVertices(2),:);
        
        % Find direction of the face
        FaceVector = VB-VA;
        FaceLength = norm(FaceVector);
        FaceDir = FaceVector/FaceLength;

        % Determine face coincidence
        % A coincident face overlaps the fracture intersection line 
        % If coincident, we record this face as a projection face
        % But new pEDFM connections will have a zero length, as all
        % connections are already provided by the EDFM method alone
        IsCoincident = is_coincident_face(P0, P1, VA, VB);

        % Add pEDFM connections to all cells in iFrac2 that have an
        % EDFM connection with this cell
        EDFMCellID2ThisCell = find(EDFMCellID1 == CellID);

        % Project EDFM lengths onto projection face
        % This will be the length associated with the pEDFM connection
        % Note: Cells that are coincident with internal fracture mesh
        % edge/face will have a 0 projected length, effectively nullifying
        % the influence of pEDFM connection. This method is reflected in
        % the rock-frac pEDFM connections as well. 
        ProjectedLength = EDFMLength(EDFMCellID2ThisCell) ...
                          * abs(dot(LineDir,FaceDir)) ...
                          * (1.0 - IsCoincident);

        % Generate pEDFM connections
        for iEDFMConn = 1:numel(EDFMCellID2ThisCell)
            
            % Find cell IDs that neighbor cell should connect with
            % inside iFrac2.
            % This is the cell from iFrac2 that this cell already has
            % EDFM connection with. There may be multiple.
            PEDFMCellID = EDFMCellID2(EDFMCellID2ThisCell(iEDFMConn));

            % Cell ID output
            FracFracCellIDs(end+1,:) = [NeighborCellID, PEDFMCellID];

            % Save length of conneciton
            LengthPEDFM = [LengthPEDFM; ProjectedLength(iEDFMConn)];

        end

    end

end

end

function [Neighbors, NeighborFaces] = build_cell_neighbors(F)
    % Build neighboring-cell connectivity using shared face IDs.

    % Every internal face belongs to exactly two cells.
    % Boundary faces belong to one cell.

    nCells = size(F,1);

    Neighbors = cell(nCells,1);
    NeighborFaces = cell(nCells,1);

    faceIDs = F(:);

    cellIDs = repmat((1:nCells)',size(F,2),1);

    valid = faceIDs > 0;

    faceIDs = faceIDs(valid);
    cellIDs = cellIDs(valid);

    [uniqueFaces,~,group] = unique(faceIDs);

    for k = 1:numel(uniqueFaces)

        cells = cellIDs(group == k);

        if numel(cells) == 1
            % External boundary face
            continue
        elseif numel(cells) == 2
            i = cells(1);
            j = cells(2);

            faceID = uniqueFaces(k);

            Neighbors{i}(end+1) = j;
            Neighbors{j}(end+1) = i;

            NeighborFaces{i}(end+1) = faceID;
            NeighborFaces{j}(end+1) = faceID;
        else
            error('Face %d belongs to more than two elements.', ...
                  uniqueFaces(k));
        end

    end

end

function [intersects] = intersect_segments_3D( ...
    P0,P1,Q0,Q1)

    % Test intersection between two 3D line segments.
    
    % Handles ordinary intersections, endpoint intersections, and collinear
    % overlap. For collinear overlap, the midpoint of the overlap is 
    % returned as the intersection point.

    intersects = false;

    D = P1-P0;
    E = Q1-Q0;

    lenD = norm(D);
    lenE = norm(E);

    if lenD == 0 || lenE == 0
        return
    end

    scale = max([lenD,lenE,1]);
    tol = 1e-10*scale;

    % Parallel or collinear case
    crossDE = cross(D,E);

    if norm(crossDE) <= 1e-12*lenD*lenE

        % Check collinearity
        if norm(cross(Q0-P0,D))/lenD > tol
            return
        end

        % Project Q endpoints onto P
        t0 = dot(Q0-P0,D)/dot(D,D);
        t1 = dot(Q1-P0,D)/dot(D,D);

        tLow = max(0,min(t0,t1));
        tHigh = min(1,max(t0,t1));

        if tLow <= tHigh + tol/lenD
            intersects = true;
            tm = 0.5*(tLow+tHigh);
            Pint = P0+tm*D;
        end

        return

    end

    % Nonparallel case
    M = [D(:),-E(:)];

    x = M\(Q0(:)-P0(:));

    t = x(1);
    s = x(2);

    PA = P0+t*D;
    PB = Q0+s*E;

    residual = norm(PA-PB);

    if residual <= tol && ...
       t >= -tol/lenD && t <= 1+tol/lenD && ...
       s >= -tol/lenE && s <= 1+tol/lenE

        intersects = true;

    end

end

function coincident = is_coincident_face(P0,P1,VA,VB)

    % Determine whether the fracture intersection line overlaps a mesh face
    % over a finite physical length. Isolated point intersections are NOT 
    % coincidence.

    coincident = false;

    D = P1-P0;
    E = VB-VA;

    LD = norm(D);
    LE = norm(E);

    if LD == 0 || LE == 0
        return
    end

    scale = max([LD,LE,1]);
    tol = 1e-10*scale;

    % Check parallelism
    if norm(cross(D,E)) > tol*LE
        return
    end

    % Check collinearity
    if norm(cross(VA-P0,D))/LD > tol
        return
    end

    % Calculate overlap along intersection line
    tA = dot(VA-P0,D)/dot(D,D);
    tB = dot(VB-P0,D)/dot(D,D);

    tLow = max(0,min(tA,tB));
    tHigh = min(1,max(tA,tB));

    overlapLength = (tHigh-tLow)*LD;

    % Require positive-length overlap
    coincident = overlapLength > tol;

    if coincident
        coincident = 1;
    else
        coincident = 0;
    end

end

function pEDFM = empty_pedfm()

    pEDFM.FracID1 = zeros(0,1);
    pEDFM.FracID2 = zeros(0,1);
    pEDFM.CellID1 = zeros(0,1);
    pEDFM.CellID2 = zeros(0,1);
    pEDFM.Length = zeros(0,1);

end

function ProjectionFaces = empty_projection_faces()

    ProjectionFaces.FracID = zeros(0,1);
    ProjectionFaces.FaceID = zeros(0,1);
    ProjectionFaces.IsCoincident = false(0,1);

end
