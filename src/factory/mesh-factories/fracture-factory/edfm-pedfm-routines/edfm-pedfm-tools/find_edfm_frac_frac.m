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
%   FracMesh:           struct
%                       Contains relevant frac mesh information
%   iFrac1:             int
%                       Index of fracture 1
%   iFrac2:             int
%                       Index of fracture 2
%   IntLineEndPoints:   array (2-by-3)
%                       Each row gives the (x,y,z) coords of end point of
%                       the intersection line between iFrac1 and iFrac2
%
% Output:
%   IntersectOutput:    struct; with 8 fields
%                       Information regarding EDFM connections frac-frac
%                       SegmentID
%                       FracID1
%                       FracID2
%                       CellID1
%                       CellID2
%                       LengthEDFM
%                       OnEdge1
%                       OnEdge2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function IntersectOutput = ...
    find_edfm_frac_frac(FracMesh, ...
                        iFrac1, ...
                        iFrac2, ...
                        IntLineEndPoints)

% Initialize output arrays
SegmentID = []; % Segment of intersec line w/ unique f1-f2 cell pairing
FracID1   = []; % Frac ID of f1 attached to segment
FracID2   = []; % Frac ID of f2 attached to segment
CellID1   = []; % Cell ID of f1 cell attached to segment
CellID2   = []; % Cell ID of f2 cell attached to segment

LengthEDFM = []; % Length of segment 

OnEdge1 = []; % Does the segment coincide with edge in f1
OnEdge2 = []; % Does the segment coincide with edge in f2

tStart = [];
tEnd   = [];

PStart = zeros(0,3);
PEnd   = zeros(0,3);

% Intersection geometry
P0 = IntLineEndPoints(1,:);
P1 = IntLineEndPoints(2,:);

D = P1 - P0;
L = norm(D);

if L == 0
    error('Intersection line has zero length.');
end

% Build local grids
Grid1 = build_local_grid(FracMesh(iFrac1));
Grid2 = build_local_grid(FracMesh(iFrac2));

% Project intersection endpoints
[u10,v10] = project_point(Grid1,P0);
[u11,v11] = project_point(Grid1,P1);

[u20,v20] = project_point(Grid2,P0);
[u21,v21] = project_point(Grid2,P1);

% Find crossings with internal grid lines
t1 = get_grid_crossings(Grid1,u10,v10,u11,v11);
t2 = get_grid_crossings(Grid2,u20,v20,u21,v21);

% Merge crossing parameters
t = sort([0; 1; t1(:); t2(:)]);

% Parameter tolerance corresponding to a smallfraction of the smallest
% cell dimension
hMin = min([diff(Grid1.uEdges), ...
            diff(Grid1.vEdges), ...
            diff(Grid2.uEdges), ...
            diff(Grid2.vEdges)]);

tolT = min(1e-10, 1e-8*hMin/L);

tUnique = t(1);

for k = 2:numel(t)

    if t(k)-tUnique(end) > tolT
        tUnique(end+1,1) = t(k);
    end

end

tUnique(1) = 0;
tUnique(end) = 1;

% Process geometric segments
for s = 1:numel(tUnique)-1

    ta = tUnique(s);
    tb = tUnique(s+1);

    tm = 0.5*(ta+tb);

    % Segment endpoints and midpoint
    PA = P0 + ta*D;
    PB = P0 + tb*D;
    PM = P0 + tm*D;

    segLength = (tb-ta)*L;

    % Find all cells touching this segment
    [cells1,edge1] = locate_cells_multi( ...
        Grid1,PM,PA,PB);

    [cells2,edge2] = locate_cells_multi( ...
        Grid2,PM,PA,PB);

    if isempty(cells1) || isempty(cells2)
        error(['Unable to identify cells for ' ...
               'intersection segment %d.'],s);
    end

    % Construct all cell pair combinations
    [A,B] = ndgrid(cells1,cells2);
    pairs = [A(:),B(:)];
    nPairs = size(pairs,1);

    % Append results
    SegmentID = [SegmentID; repmat(s,nPairs,1)];

    FracID1 = [FracID1; repmat(iFrac1,nPairs,1)];
    FracID2 = [FracID2; repmat(iFrac2,nPairs,1)];

    CellID1 = [CellID1; pairs(:,1)];
    CellID2 = [CellID2; pairs(:,2)];

    LengthEDFM = [LengthEDFM; repmat(segLength,nPairs,1)];

    OnEdge1 = [OnEdge1; repmat(edge1,nPairs,1)];
    OnEdge2 = [OnEdge2; repmat(edge2,nPairs,1)];

    tStart = [tStart; repmat(ta,nPairs,1)];
    tEnd   = [tEnd; repmat(tb,nPairs,1)];

    PStart = [PStart; repmat(PA,nPairs,1)];
    PEnd   = [PEnd; repmat(PB,nPairs,1)];

end

% Construct output
IntersectOutput = struct();
IntersectOutput.SegmentID = SegmentID;
IntersectOutput.FracID1 = FracID1;
IntersectOutput.FracID2 = FracID2;
IntersectOutput.CellID1 = CellID1;
IntersectOutput.CellID2 = CellID2;
IntersectOutput.LengthEDFM = LengthEDFM;
IntersectOutput.OnEdge1 = OnEdge1;
IntersectOutput.OnEdge2 = OnEdge2;

end


function Grid = build_local_grid(Frac)

V = Frac.vertex.coords;
E = Frac.elements.vertices;
C = Frac.elements.centroid;

% Determine local axes from first quadrilateral
verts = V(E(1,:),:);

origin = verts(1,:);

e1 = verts(2,:) - verts(1,:);
e2 = verts(4,:) - verts(1,:);

U = e1/norm(e1);

e2 = e2-dot(e2,U)*U;
W = e2/norm(e2);

% Project vertices
dV = V-origin;

uV = dV*U';
vV = dV*W';

% Project centroids
dC = C-origin;

uC = dC*U';
vC = dC*W';

% Estimate tolerance from mesh dimensions
scale = max([max(uV)-min(uV), max(vV)-min(vV), 1]);

tol = 1e-10*scale;

% Unique centroid positions
uCenters = unique_tol(uC,tol);
vCenters = unique_tol(vC,tol);

Nu = numel(uCenters);
Nv = numel(vCenters);

if Nu*Nv ~= size(E,1)
    error('Mesh is not a complete rectangular grid.');
end

% Grid boundaries
uEdges = linspace(min(uV),max(uV),Nu+1);
vEdges = linspace(min(vV),max(vV),Nv+1);

% Build lookup matrix
CellMap = zeros(Nu,Nv);

for k = 1:size(E,1)

    [~,i] = min(abs(uCenters-uC(k)));
    [~,j] = min(abs(vCenters-vC(k)));

    if CellMap(i,j) ~= 0
        error('Duplicate cell location in grid map.');
    end

    CellMap(i,j) = k;

end

% Store grid
Grid.origin = origin;

Grid.U = U;
Grid.W = W;

Grid.uEdges = uEdges;
Grid.vEdges = vEdges;

Grid.Nu = Nu;
Grid.Nv = Nv;

Grid.CellMap = CellMap;

% Coordinate tolerance, kept well below cell spacing
hMin = min([diff(uEdges),diff(vEdges)]);

Grid.tol = min(1e-8*hMin, ...
               1e-10*scale);

end

function values = unique_tol(x,tol)

    x = sort(x(:));
    values = x(1);

    for k = 2:numel(x)

        if abs(x(k)-values(end)) > tol
            values(end+1,1) = x(k);
        end

    end

end

function [u,v] = project_point(Grid,P)

    dP = P-Grid.origin;

    u = dP*Grid.U';
    v = dP*Grid.W';

end

function tCross = get_grid_crossings( ...
    Grid,u0,v0,u1,v1)

    du = u1-u0;
    dv = v1-v0;

    tCross = [];

    % Crossing parameters for constant-u grid lines
    if abs(du) > Grid.tol
        uLines = Grid.uEdges(2:end-1);
        tu = (uLines-u0)/du;
        tu = tu(tu > 0 & tu < 1);
        tCross = [tCross;tu(:)];
    end

    % Crossing parameters for constant-v grid lines
    if abs(dv) > Grid.tol
        vLines = Grid.vEdges(2:end-1);
        tv = (vLines-v0)/dv;
        tv = tv(tv > 0 & tv < 1);
        tCross = [tCross;tv(:)];
    end

end

function [CellIDs,OnEdge] = locate_cells_multi( ...
    Grid,PM,PA,PB)

    % Returns all cells touching the segment over a finite
    % length, including both sides of coincident internal edges.
    %
    % PM = segment midpoint
    % PA = segment start
    % PB = segment end

    % Project segment points
    [u,v] = project_point(Grid,PM);

    [ua,va] = project_point(Grid,PA);
    [ub,vb] = project_point(Grid,PB);

    tol = Grid.tol;

    % Find primary cell indices
    i = discretize(u,Grid.uEdges);
    j = discretize(v,Grid.vEdges);

    if isnan(i) || isnan(j)
        CellIDs = [];
        OnEdge = false;
        return
    end

    % Determine whether segment follows internal grid lines
    onUEdge = false;
    onVEdge = false;

    iu = [];
    jv = [];

    % Check internal constant-u grid lines
    for k = 2:numel(Grid.uEdges)-1

        uk = Grid.uEdges(k);

        if abs(u-uk) <= tol && ...
           abs(ua-uk) <= tol && ...
           abs(ub-uk) <= tol

            onUEdge = true;
            iu = k;
            break

        end

    end

    % Check internal constant-v grid lines
    for k = 2:numel(Grid.vEdges)-1

        vk = Grid.vEdges(k);

        if abs(v-vk) <= tol && ...
           abs(va-vk) <= tol && ...
           abs(vb-vk) <= tol

            onVEdge = true;
            jv = k;
            break

        end

    end

    % Determine touching grid indices
    if onUEdge
        % Cells immediately left and right of grid line
        iList = [iu-1,iu];
    else
        iList = i;
    end

    if onVEdge
        % Cells immediately below and above grid line
        jList = [jv-1,jv];
    else
        jList = j;
    end

    % Construct all touching cells
    [I,J] = ndgrid(iList,jList);
    ind = sub2ind( ...
        size(Grid.CellMap),I(:),J(:));
    CellIDs = Grid.CellMap(ind);
    CellIDs = unique(CellIDs(:),'stable');

    % Edge coincidence flag
    OnEdge = onUEdge || onVEdge;

end