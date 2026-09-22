%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Ryan Haagenson; Artur Castiel
%TU Delft
%
% Finding the intersection between a fracture-fracture intersection line
% and the edges of a fracture mesh. The function also identifies the
% fracture cells containing the endpoints of the intersection line.
%
% Input:
%   p1:             array (N x 3)
%                   First point defining the fracture intersection line
%   p2:             array (N x 3)
%                   Second point defining the fracture intersection line
%   normal:         array (N x 3)
%                   Normal vector of the fracture plane
%   fracmesh:       struct
%                   Contains fracture mesh information
%
% Output:
%   numIntersection: Sparse double column vector (N x 1)
%                    Number of intersections between the fracture
%                    intersection line and fracture mesh edges
%   p1_3d:           Sparse double column vector (N x 3)
%                    First intersection point in 3D coordinates
%   p2_3d:           Sparse double column vector (N x 3)
%                    Second intersection point in 3D coordinates
%   v1:              array
%                    Fracture cells containing first intersection point
%   v2:              array
%                    Fracture cells containing second intersection point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [numIntersection, p1_3d, p2_3d, v1, v2] = ...
         calculate_vol_edge_intersection(p1, p2, normal, fracmesh)

% Fracture mesh coordinates
coords = fracmesh.vertex.coords;

% Fracture mesh edge connectivity
adjList = [fracmesh.faces.inner.vertices; ...
           fracmesh.faces.bound.vertices];

% Center of fracture mesh used as the reference point for the projection
center = mean(coords);

% Project the fracture mesh and fracture intersection line to 2D.
% This allows all edge intersection calculations to be performed in the
% plane of the fracture rather than directly in 3D.
xyCoords = proj_coplanar_coords_on_xy(center, normal, coords, 1);
xyP1 = proj_coplanar_coords_on_xy(center, normal, p1, 1);
xyP2 = proj_coplanar_coords_on_xy(center, normal, p2, 1);

% Find where the fracture intersection line crosses the fracture mesh
% edges in the projected 2D coordinate system.
[node1, node2, numIntersection] = ...
    edge_analysis_in_2d(adjList, xyCoords, xyP1, xyP2);

% Identify the fracture cells containing the endpoints of the
% intersection line. Bounding boxes are used first to reduce the number
% of cells requiring a polygon containment test.
[v1, v2] = find_in_volume_frac_line(xyP1, xyP2, xyCoords, ...
                                    fracmesh.elements.vertices);

% Add the third coordinate required to transform the projected points
% back to the original fracture plane.
node1(numIntersection > 0,3) = center(:,3);
node2(numIntersection > 1,3) = center(:,3);

% Initialize sparse arrays for the intersection points in 3D.
p1_3d = sparse(size(node1,1),3);
p2_3d = sparse(size(node2,1),3);

% Transform the first intersection points from the projected coordinate
% system back to the original 3D fracture plane.
ref = numIntersection > 0;
p1_3d(ref,:) = ...
    proj_coplanar_coords_on_xy(center, normal, node1(ref,:), 0);

% Transform the second intersection points back to 3D.
ref = numIntersection > 0;
p2_3d(ref,:) = ...
    proj_coplanar_coords_on_xy(center, normal, node2(ref,:), 0);

end

function [v1, v2] = find_in_volume_frac_line(p1, p2, coords, volNodes)

    tolerance = 10^-12;

    % Generate bounding boxes for all fracture cells. The bounding boxes
    % are used to exclude cells that cannot contain the intersection point.
    [minVec, maxVec] = generate_boundbox(volNodes, coords, false);

    % Identify cells whose bounding boxes contain the first intersection
    % point in the projected 2D coordinate system.
    refP1X = (minVec(:,1) - tolerance <= p1(:,1)) & ...
             (maxVec(:,1) + tolerance >= p1(:,1));
    refP1Y = (minVec(:,2) - tolerance <= p1(:,2)) & ...
             (maxVec(:,2) + tolerance >= p1(:,2));
    refP1 = refP1X & refP1Y;

    % Identify cells whose bounding boxes contain the second intersection
    % point in the projected 2D coordinate system.
    refP2X = (minVec(:,1) - tolerance <= p2(:,1)) & ...
             (maxVec(:,1) + tolerance >= p2(:,1));
    refP2Y = (minVec(:,2) - tolerance <= p2(:,2)) & ...
             (maxVec(:,2) + tolerance >= p2(:,2));
    refP2 = refP2X & refP2Y;

    % Test the first intersection point against candidate fracture cells.
    % A point is considered part of a cell if it lies inside the polygon
    % or directly on one of its edges.
    v1 = [];

    for localVol = find(refP1)'
        localNodes = volNodes(localVol,:);
        polygon = coords(localNodes,:);

        [flagIn, flagOn] = ...
            inpolygon(p1(:,1), p1(:,2), polygon(:,1), polygon(:,2));

        if flagIn || flagOn
            v1 = [v1; localVol];
        end
    end

    % Test the second intersection point against candidate fracture cells.
    v2 = [];

    for localVol = find(refP2)'
        localNodes = volNodes(localVol,:);
        polygon = coords(localNodes,:);

        [flagIn, flagOn] = ...
            inpolygon(p2(:,1), p2(:,2), polygon(:,1), polygon(:,2));

        if flagIn || flagOn
            v2 = [v2; localVol];
        end
    end

end

function [intP1, intP2, numIntersection] = ...
         edge_analysis_in_2d(adj, coords, fracP1, fracP2)

    tolerance = 10^-13;

    % Initialize intersection points and number of intersections for
    % every fracture mesh edge.
    intP1 = sparse(size(adj,1),2);
    intP2 = sparse(size(adj,1),2);
    numIntersection = sparse(size(adj,1),1);

    % Coordinates of the endpoints of each fracture mesh edge
    p1 = coords(adj(:,1),1:2);
    p2 = coords(adj(:,2),1:2);

    % Coordinates of the endpoints of the fracture-fracture
    % intersection line
    fracP1 = fracP1(:,1:2);
    fracP2 = fracP2(:,1:2);

    % Repeat the fracture intersection line for every fracture mesh edge
    % so that all edge-line intersections can be calculated simultaneously.
    refFracP1 = repmat(fracP1, size(p1,1), 1);
    refFracP2 = repmat(fracP2, size(p2,1), 1);

    % Define parametric equations for the mesh edges and fracture line.
    % Edge:      p = edgeOrigin + t1 * edgeVector
    % Fracture:  p = fractureOrigin + t2 * fractureVector
    edgeOrigin = p1;
    edgeVector = p2 - p1;

    fractureOrigin = refFracP1;
    fractureVector = refFracP2 - refFracP1;

    % Identify mesh edges that are parallel to the fracture intersection
    % line. Parallel edges require separate treatment because the standard
    % two-line intersection calculation has a zero denominator.
    edgeVectorNorm = abs(edgeVector ./ vecnorm(edgeVector,2,2));
    fractureVectorNorm = ...
        abs(fractureVector ./ vecnorm(fractureVector,2,2));

    refParallel = ...
        all(abs(edgeVectorNorm - fractureVectorNorm) < tolerance, 2);

    refCross = ~refParallel;

    % Find intersections between non-parallel mesh edges and the fracture
    % intersection line by solving the two parametric equations.
    [t1, t2] = ...
        solve_system_parametric(edgeVector(refCross,:), ...
                                fractureVector(refCross,:), ...
                                fractureOrigin(refCross,:) - ...
                                edgeOrigin(refCross,:));

    % Retain only intersections that occur within both line segments.
    refT1 = (t1 >= (0 - tolerance)) & ...
            (t1 <= (1 + tolerance));
    refT2 = (t2 >= (0 - tolerance)) & ...
            (t2 <= (1 + tolerance));
    ref = refT1 & refT2;

    % Each valid non-parallel edge has a single intersection point.
    refCross(refCross) = ref;
    numIntersection(refCross) = 1;

    % Calculate the intersection point using the mesh-edge parameter.
    intP1(refCross,:) = ...
        edgeOrigin(refCross,:) + ...
        t1(ref,:) .* (p2(refCross,:) - p1(refCross,:));

    % For parallel edges, determine whether the mesh edge lies on the
    % fracture intersection line. If so, the intersection is an interval
    % rather than a single point and both endpoints are retained.
    tP1 = interpolate_t(p1(refParallel,:), ...
                         fractureOrigin(refParallel,:), ...
                         fractureVector(refParallel,:));

    % A parallel mesh edge lies on the fracture line when its endpoint
    % parameters calculated from the fracture line are equal.
    refTP1 = abs(tP1(:,1) - tP1(:,2)) < tolerance;
    tP1 = mean(tP1,2);

    % Repeat the test using the second mesh-edge endpoint. This ensures
    % that both endpoints of the mesh edge lie on the fracture line.
    tP2 = interpolate_t(p2(refParallel,:), ...
                         fractureOrigin(refParallel,:), ...
                         edgeVector(refParallel,:));

    refTP2 = abs(tP2(:,1) - tP2(:,2)) < tolerance;
    tP2 = mean(tP2,2);

    % Retain parallel edges that are coincident with the fracture
    % intersection line.
    refFracEdge = refTP1 & refTP2;

    % Bound the intersection interval to the fracture intersection line
    % and ensure that the two endpoints are ordered consistently.
    [tP1(refFracEdge), tP2(refFracEdge)] = ...
        sort_and_bound_interval(tP1(refFracEdge), ...
                                tP2(refFracEdge));

    refParallel(refParallel) = refFracEdge;

    % A coincident edge produces two intersection points corresponding to
    % the endpoints of the overlapping interval.
    numIntersection(refParallel) = 2;

    intP1(refParallel,:) = ...
        fractureOrigin(refParallel,:) + ...
        tP1(refFracEdge) .* ...
        (fractureVector(refParallel,:));

    intP2(refParallel,:) = ...
        edgeOrigin(refParallel,:) + ...
        tP2(refFracEdge) .* ...
        (fractureVector(refParallel,:));

end

function [newT1, newT2] = sort_and_bound_interval(t1, t2)

    % Bound intersection parameters to the fracture intersection line.
    auxMat = [t1 t2];
    auxMat(auxMat > 1) = 1;
    auxMat(auxMat < 0) = 0;

    % Sort the parameters so that newT1 always represents the beginning
    % of the intersection interval and newT2 the end.
    auxMat = sort(auxMat,2);

    newT1 = auxMat(:,1);
    newT2 = auxMat(:,2);

end

function [t1, t2] = solve_system_parametric(m, n, q)

    % Split the two-dimensional parametric system into scalar components.
    m1 = m(:,1);
    m2 = m(:,2);
    n1 = n(:,1);
    n2 = n(:,2);
    q1 = q(:,1);
    q2 = q(:,2);

    % Determinant of the coefficient matrix. This is non-zero for the
    % non-parallel edges handled by this function.
    denominator = m1 .* n2 - m2 .* n1;

    % Solve for the parameters describing the intersection point on
    % each of the two parametric lines.
    t1 = (n2 .* q1) ./ denominator - ...
         (n1 .* q2) ./ denominator;

    t2 = (m2 .* q1) ./ denominator - ...
         (m1 .* q2) ./ denominator;

end

function t = interpolate_t(x, a, b)

    % Calculate the parametric coordinate of x along the line
    % x = a + t*b.
    t = (x - a) ./ b;

end