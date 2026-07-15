classdef translateMeshCore
    %translateMeshCore Creates a compatibility layers for the meshcore/
    %  This layer allows MeshCore to be used with the pEDFM and EDFM
    %  routines.

    properties
        vertex
        face
        element
        props
        core
    end

    methods
        function obj = translateMeshCore(core, ex)
            obj.core = core;
            obj.vertex = transalte_vertex(core);
            obj.face = transalte_face(core);
            obj.element = translate_element(core,obj.face.num);
            obj.props = translate_props(core);
        end        

    end
end

function props = translate_props(core)
    if isprop(core,"props")        
        props = core.props;
    else
        props  = struct();
    end
    % if  isempty(core.props)
    %     props.perm = ones(size(core.cells.centroid,1),1);
    % end
end


function element = translate_element(core,face_num)
element = struct();
%% translating vertices
element.vertices = ...
    sparse_to_full_zeromat(core.cells.element_nodes_bool);
for c = 1:size(element.vertices,1)
    v = element.vertices(c,:);
    v = v(v>0);
    if numel(v)==8
        cellFaceIds = find(core.cells.element_faces_bool(c,:));
        v = reorder_hex8_for_vtk( ...
            v, ...
            cellFaceIds, ...
            core.faces.face_nodes, ...
            core.nodes.coords);
        element.vertices(c,1:8) = v;
    end
end
%% defining the number of elements
element.num = size(element.vertices,1);
%% creating face translation vectors
trans_face = zeros(size(core.faces.bound_faces));
num_ifaces = sum(~core.faces.bound_faces);
trans_face(~core.faces.bound_faces) = 1:num_ifaces;
trans_face(core.faces.bound_faces) = ...
    num_ifaces+1:size(trans_face,1);
%% translating element_face
auxmat = sparse_to_full_zeromat(core.cells.element_faces_bool);
auxmat(auxmat ~= 0) = trans_face(auxmat(auxmat ~= 0));
element.faces = auxmat;
%% defining element faces bool
element.faces_bool = ...
    create_element_face_bool(element.faces,...
    element.num, face_num);
%% translating volume
element.volume = core.cells.volume;
%% translating centroid
element.centroid = core.cells.centroid;
end

function [faces_bool] = create_element_face_bool(faces,num_cells,num_faces)
%% lacks this
num_rep = size(faces,2);
ind_mat = [1:size(faces,1)]';
ii = repelem(ind_mat,num_rep,1);
jj = reshape(faces',[],1);
ref = (ii~=0) & (jj~=0);
faces_bool = sparse(ii(ref),jj(ref),true, num_cells,num_faces);
end

function face = transalte_face(core)
face = struct();
inner = struct();
bound = struct();
%% Defining face vertices
bound.vertices = ...
    core.faces.face_nodes(core.faces.bound_faces,:);
inner.vertices = ...
    core.faces.face_nodes(~core.faces.bound_faces,:);
%% Defining num of elements
bound.num = size(bound.vertices,1);
inner.num = size(inner.vertices,1);
face.num = bound.num + inner.num;
%% Defining area
bound.area = ...
    core.faces.areas(core.faces.bound_faces,:);
inner.area = ...
    core.faces.areas(~core.faces.bound_faces,:);
%% Defining normals (pEDFM requires area vector normal)
% bound.normal = bound.area .* ...
%     core.faces.normals(core.faces.bound_faces,:);
% inner.normal = inner.area .* ...
%     core.faces.normals(~core.faces.bound_faces,:);

bound.normal = core.faces.normals(core.faces.bound_faces,:);

inner.normal = core.faces.normals(~core.faces.bound_faces,:);
%% Defining centroid
bound.centroid = ...
    core.faces.centroid(core.faces.bound_faces,:);
inner.centroid = ...
    core.faces.centroid(~core.faces.bound_faces,:);
%% Defining neighbors
bound.upstream = ...
    core.faces.face_neigh(core.faces.bound_faces,1);
inner.upstream = ...
    core.faces.face_neigh(~core.faces.bound_faces,1);
inner.downstream = ...
    core.faces.face_neigh(~core.faces.bound_faces,2);
%% defining face
face.bound = bound;
face.inner = inner;
end

function vertex = transalte_vertex(core)
%METHOD1 Summary of this method goes here
%   Detailed explanation goes here
vertex = struct();
vertex.coord = core.nodes.coords;
vertex.num = core.nodes.num;
end

function v_out = reorder_hex8_for_vtk( ...
        cellNodeIds, cellFaceIds, allFaceNodes, nodeCoords)
    %REORDER_HEX8_FOR_VTK Reconstruct general HEX8 connectivity for VTK.
    %
    % This method is independent of:
    %   - global coordinate axes,
    %   - extrusion direction,
    %   - element rotation,
    %   - element skewness,
    %   - node ID ordering.
    %
    % It reconstructs the hexahedron's edge graph from its six quadrilateral
    % faces, then selects a right-handed VTK_HEXAHEDRON ordering.
    %
    % VTK ordering:
    %
    %          8----------7
    %         /|         /|
    %        5----------6 |
    %        | |        | |
    %        | 4--------|-3
    %        |/         |/
    %        1----------2
    %
    % Faces:
    %   [1 2 3 4]
    %   [5 6 7 8]
    %   [1 2 6 5]
    %   [2 3 7 6]
    %   [3 4 8 7]
    %   [4 1 5 8]

    cellNodeIds = cellNodeIds(:).';
    cellFaceIds = cellFaceIds(:).';

    if numel(cellNodeIds) ~= 8
        error('HEX8 reorder requires exactly eight cell nodes.');
    end

    if numel(unique(cellNodeIds)) ~= 8
        error('HEX8 cell contains repeated node IDs.');
    end

    if numel(cellFaceIds) ~= 6
        error( ...
            'HEX8 reorder requires six faces, but cell has %d faces.', ...
            numel(cellFaceIds));
    end

    % ---------------------------------------------------------------
    % Obtain the six face-node sets in cell-local node numbering.
    % ---------------------------------------------------------------
    faceLocal = zeros(6,4);

    for f = 1:6
        globalFaceNodes = allFaceNodes(cellFaceIds(f), :);
        globalFaceNodes = globalFaceNodes(globalFaceNodes > 0);

        if numel(globalFaceNodes) ~= 4
            error( ...
                ['HEX8 face %d has %d nodes. A first-order HEX8 cell ' ...
                 'must have six quadrilateral faces.'], ...
                cellFaceIds(f), numel(globalFaceNodes));
        end

        [isInCell, localIds] = ismember(globalFaceNodes, cellNodeIds);

        if ~all(isInCell)
            error( ...
                'HEX8 face %d contains a node not belonging to the cell.', ...
                cellFaceIds(f));
        end

        if numel(unique(localIds)) ~= 4
            error('HEX8 face contains repeated nodes.');
        end

        faceLocal(f,:) = localIds;
    end

    % Every HEX8 node must belong to exactly three cell faces.
    nodeFaceCount = accumarray(faceLocal(:), 1, [8,1]);

    if any(nodeFaceCount ~= 3)
        error( ...
            ['Invalid HEX8 topology: every corner must belong to exactly ' ...
             'three faces. Face counts are [%s].'], ...
            sprintf('%d ', nodeFaceCount));
    end

    % ---------------------------------------------------------------
    % Build the HEX8 edge graph.
    %
    % Adjacent HEX8 vertices share exactly two cell faces.
    % Face diagonals share one; opposite vertices share zero.
    % ---------------------------------------------------------------
    sharedFaceCount = zeros(8,8);

    for f = 1:6
        nodes = faceLocal(f,:);

        for i = 1:3
            for j = i+1:4
                a = nodes(i);
                b = nodes(j);

                sharedFaceCount(a,b) = sharedFaceCount(a,b) + 1;
                sharedFaceCount(b,a) = sharedFaceCount(b,a) + 1;
            end
        end
    end

    edgeGraph = sharedFaceCount == 2;
    edgeGraph(1:9:end) = false;

    nodeDegree = sum(edgeGraph, 2);

    if any(nodeDegree ~= 3)
        error( ...
            ['Invalid HEX8 topology: each vertex must have exactly three ' ...
             'edge-connected neighbours. Degrees are [%s].'], ...
            sprintf('%d ', nodeDegree));
    end

    % A valid HEX8 has exactly 12 edges.
    if nnz(triu(edgeGraph,1)) ~= 12
        error( ...
            'Invalid HEX8 topology: expected 12 edges, found %d.', ...
            nnz(triu(edgeGraph,1)));
    end

    % ---------------------------------------------------------------
    % Generate all orderings originating at each possible VTK node 1.
    %
    % At VTK node 1, the three graph neighbours become nodes 2, 4 and 5.
    % Their permutation controls the local coordinate directions.
    % ---------------------------------------------------------------
    candidates = zeros(0,8);
    quality    = zeros(0,1);

    cellPts = nodeCoords(cellNodeIds,:);

    cellScale = max(vecnorm(cellPts - mean(cellPts,1), 2, 2));

    if cellScale <= eps
        error('HEX8 cell has zero geometric extent.');
    end

    determinantTolerance = 1e-12 * cellScale^3;

    for vtk1 = 1:8
        neighbours = find(edgeGraph(vtk1,:));

        neighbourPermutations = perms(neighbours);

        for k = 1:size(neighbourPermutations,1)
            vtk2 = neighbourPermutations(k,1);
            vtk4 = neighbourPermutations(k,2);
            vtk5 = neighbourPermutations(k,3);

            vtk3 = uniqueCommonNeighbour( ...
                vtk2, vtk4, vtk1, edgeGraph);

            vtk6 = uniqueCommonNeighbour( ...
                vtk2, vtk5, vtk1, edgeGraph);

            vtk8 = uniqueCommonNeighbour( ...
                vtk4, vtk5, vtk1, edgeGraph);

            if any([vtk3 vtk6 vtk8] == 0)
                continue
            end

            used = [vtk1 vtk2 vtk3 vtk4 vtk5 vtk6 vtk8];

            vtk7Candidates = setdiff(1:8, used);

            if numel(vtk7Candidates) ~= 1
                continue
            end

            vtk7 = vtk7Candidates;

            ordering = ...
                [vtk1 vtk2 vtk3 vtk4 vtk5 vtk6 vtk7 vtk8];

            if ~isValidVtkHexGraph(ordering, edgeGraph)
                continue
            end

            if ~facesMatchVtkHex(ordering, faceLocal)
                continue
            end

            % Local right-handedness at VTK corner 1.
            p1 = cellPts(vtk1,:);
            p2 = cellPts(vtk2,:);
            p4 = cellPts(vtk4,:);
            p5 = cellPts(vtk5,:);

            detJ = dot(cross(p2-p1, p4-p1), p5-p1);

            if detJ <= determinantTolerance
                continue
            end

            candidates(end+1,:) = ordering; %#ok<AGROW>
            quality(end+1,1) = detJ; %#ok<AGROW>
        end
    end

    if isempty(candidates)
        error( ...
            ['Could not construct a right-handed VTK ordering for this ' ...
             'HEX8 cell. The cell may be inverted, self-intersecting, ' ...
             'topologically invalid, or nearly degenerate.']);
    end

    % Remove duplicate local orderings.
    [candidates, uniqueRows] = unique(candidates, 'rows', 'stable');
    quality = quality(uniqueRows);

    % ---------------------------------------------------------------
    % Choose deterministically.
    %
    % Multiple VTK orderings can describe the same physical hexahedron.
    % First prefer the candidate with the strongest local determinant.
    % For essentially equal candidates, use global node IDs as a stable
    % tie-break.
    % ---------------------------------------------------------------
    maxQuality = max(quality);
    qualityTolerance = 1e-10 * maxQuality;

    bestCandidates = find(quality >= maxQuality - qualityTolerance);

    candidateGlobalIds = zeros(numel(bestCandidates),8);

    for i = 1:numel(bestCandidates)
        candidateGlobalIds(i,:) = ...
            cellNodeIds(candidates(bestCandidates(i),:));
    end

    [~, lexOrder] = sortrows(candidateGlobalIds);
    best = bestCandidates(lexOrder(1));

    v_out = cellNodeIds(candidates(best,:));

    % ---------------------------------------------------------------
    % Geometric validity check at all eight VTK corners.
    % ---------------------------------------------------------------
    validateHexJacobians(v_out, nodeCoords, cellScale);
end


function common = uniqueCommonNeighbour(a, b, excluded, edgeGraph)
% Find the cube vertex connected to both a and b, excluding one known node.

    candidates = find(edgeGraph(a,:) & edgeGraph(b,:));
    candidates(candidates == excluded) = [];

    if numel(candidates) == 1
        common = candidates;
    else
        common = 0;
    end
end


function valid = isValidVtkHexGraph(ordering, edgeGraph)
% Verify that the candidate has exactly the twelve VTK HEX8 edges.

    vtkEdges = [
        1 2
        2 3
        3 4
        4 1
        5 6
        6 7
        7 8
        8 5
        1 5
        2 6
        3 7
        4 8
    ];

    expected = false(8,8);

    for e = 1:size(vtkEdges,1)
        a = ordering(vtkEdges(e,1));
        b = ordering(vtkEdges(e,2));

        expected(a,b) = true;
        expected(b,a) = true;
    end

    valid = isequal(expected, edgeGraph);
end


function valid = facesMatchVtkHex(ordering, actualFaces)
% Verify that the candidate VTK face sets equal the cell's actual faces.

    vtkFacePositions = [
        1 2 3 4
        5 6 7 8
        1 2 6 5
        2 3 7 6
        3 4 8 7
        4 1 5 8
    ];

    expectedFaces = zeros(6,4);

    for f = 1:6
        expectedFaces(f,:) = ...
            sort(ordering(vtkFacePositions(f,:)));
    end

    actualFaces = sort(actualFaces, 2);

    expectedFaces = sortrows(expectedFaces);
    actualFaces   = sortrows(actualFaces);

    valid = isequal(expectedFaces, actualFaces);
end


function validateHexJacobians(v, nodeCoords, cellScale)
% Check the trilinear HEX8 Jacobian at all natural-coordinate corners.
%
% A positive determinant only at node 1 is insufficient for a heavily
% warped HEX8. This catches inverted or folded cells at the other corners.

    p = nodeCoords(v,:);

    naturalCorners = [
        -1 -1 -1
         1 -1 -1
         1  1 -1
        -1  1 -1
        -1 -1  1
         1 -1  1
         1  1  1
        -1  1  1
    ];

    tolerance = 1e-12 * cellScale^3;

    determinantValues = zeros(8,1);

    for sample = 1:8
        xi   = naturalCorners(sample,1);
        eta  = naturalCorners(sample,2);
        zeta = naturalCorners(sample,3);

        dN = hex8ShapeDerivatives(xi, eta, zeta);
        J  = p.' * dN;

        determinantValues(sample) = det(J);
    end

    if any(determinantValues <= tolerance)
        error( ...
            ['HEX8 has a non-positive or nearly zero Jacobian at one or ' ...
             'more corners. Determinants: [%s].'], ...
            sprintf('%.6e ', determinantValues));
    end
end


function dN = hex8ShapeDerivatives(xi, eta, zeta)
% Shape-function derivatives for the standard VTK HEX8 ordering.
%
% Rows correspond to nodes 1...8.
% Columns are derivatives with respect to xi, eta and zeta.

    signs = [
        -1 -1 -1
         1 -1 -1
         1  1 -1
        -1  1 -1
        -1 -1  1
         1 -1  1
         1  1  1
        -1  1  1
    ];

    dN = zeros(8,3);

    for i = 1:8
        sx = signs(i,1);
        sy = signs(i,2);
        sz = signs(i,3);

        dN(i,1) = ...
            0.125 * sx * (1 + sy*eta) * (1 + sz*zeta);

        dN(i,2) = ...
            0.125 * sy * (1 + sx*xi) * (1 + sz*zeta);

        dN(i,3) = ...
            0.125 * sz * (1 + sx*xi) * (1 + sy*eta);
    end
end

% function v_out = reorder_hex8_for_vtk(v, P)
%     % v: 1x8 node ids (1-based), P: POINTS array (N x 3)
%     pts = P(v, :);
%     c = mean(pts, 1);
% 
%     bx = pts(:,1) > c(1);
%     by = pts(:,2) > c(2);
%     bz = pts(:,3) > c(3);
% 
%     mask = bx + 2*by + 4*bz;  % 0..7
% 
%     % mask -> VTK_HEX position (1..8)
%     % 000->1, 100->2, 110->3, 010->4, 001->5, 101->6, 111->7, 011->8
%     lut = [1 2 4 3 5 6 8 7];
% 
%     v_out = zeros(1,8);
%     for i = 1:8
%         pos = lut(mask(i)+1);
%         if v_out(pos) ~= 0
%             error('HEX reorder failed: duplicate corner classification (degenerate cell?)');
%         end
%         v_out(pos) = v(i);
%     end
% end
