%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Ryan Haagenson; Artur Castiel
%TU Delft
%
% Calculating the average distance between fracture cells and the line of
% intersection between two fracture planes.
%
% Input:
%   frac_int_avg_d:     cell
%                       Contains average fracture cell distances to
%                       fracture-fracture intersection lines
%   test_points:        array
%                       Coordinates of points used to calculate distances
%   frac_el:            array
%                       Fracture cell indices associated with test points
%   frac_avg_d_ind:     cell
%                       Contains fracture cell indices associated with
%                       each fracture-fracture intersection
%   frac_int:           cell
%                       Contains fracture-fracture intersection data
%   plane_int:          object
%                       Contains fracture plane intersection information
%   index:              scalar
%                       Index of the current fracture
%
% Output:
%   frac_int_avg_d:     cell
%                       Updated average distances between fracture cells
%                       and fracture-fracture intersection lines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [frac_int_avg_d] = calculate_avg_frac_d(frac_int_avg_d, ...
                                                 test_points, ...
                                                 frac_el, ...
                                                 frac_avg_d_ind, ...
                                                 frac_int, ...
                                                 plane_int, ...
                                                 index)

% Identify fracture-fracture intersections
MatRef = ~cellfun(@isempty,frac_int);
[ii,jj] = find(MatRef);

% Find intersections involving the current fracture
Ref1 = ii == index;
Ref2 = jj == index;
Ref = Ref1 | Ref2;
ii = ii(Ref);
jj = jj(Ref);

% Loop through each fracture-fracture intersection
for Counter = 1:size(ii,1)

    LocalII = ii(Counter);
    LocalJJ = jj(Counter);

    % Find fracture cells associated with the current intersection
    IndVol = ...
        frac_avg_d_ind{LocalII,LocalJJ}{~(LocalII == index)+1};

    Ref = ismember(frac_el,IndVol);
    LocalFracEl = frac_el(Ref);

    % Find the line of intersection between the two fracture planes
    FracLine = ...
        cell2mat(plane_int.get_boundary_line_intersection(LocalII, ...
                                                          LocalJJ));

    % Calculate distance from each test point to the intersection line
    LAvgDist = dist_point_to_line(test_points(Ref,:), FracLine);

    % Sum distances for each fracture cell
    SumDist = accumarray(LocalFracEl,LAvgDist);

    % Count number of points for each fracture cell
    NDist = accumarray(LocalFracEl,1);

    % Calculate average distance for each fracture cell
    AvgDist = SumDist(IndVol)./NDist(IndVol);

    % Initialize output cell if necessary
    if isempty(frac_int_avg_d{LocalII,LocalJJ})
        frac_int_avg_d{LocalII,LocalJJ} = cell(2,1);
    end

    % Store average distance for current fracture
    AuxCell = frac_int_avg_d{LocalII,LocalJJ};
    AuxCell{~(LocalII == index)+1} = AvgDist;
    frac_int_avg_d{LocalII,LocalJJ} = AuxCell;

end

end

function [d] = dist_point_to_line(p,line)

    % Calculate distance from points to line in three dimensions
    % https://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
    X0 = p;
    X1 = repmat(line(1,:),size(X0,1),1);
    X2 = repmat(line(2,:),size(X0,1),1);
    
    % Area of parallelogram formed by point and line vectors
    A = vecnorm(cross((X0-X1),(X0-X2)),2,2);
    
    % Length of line
    B = norm(line(1,:) - line(2,:),2);
    
    % Distance from points to line
    d = A./B;

end
