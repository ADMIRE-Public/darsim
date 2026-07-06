%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Ryan Haagenson
%TU Delft
%
% Seeding a regularized grid of test points for fracture-fracture distance.
%
% Input:
%   fracture:       struct
%                   Contains relevant frac mesh information for given
%                   fracture
%   grid_size:      double
%                   Size of test point grid
%
% Output:
%   test_points:    array (N_test_point x 3)
%                   (x,y,z) coordinates of the gridded test points
%   frac_el:        array (N_test_point x 1)
%                   Fracture cell containing each test point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [test_points, ...
          frac_el ] = ...
          generate_testpoint_parameters_grid(fracture, ...
                                             grid_size)

% Number of fracutre cells
nFracCell = size(fracture.elements.area,1);

% Initialize output
test_points = zeros(nFracCell*grid_size^2,3);
frac_el = zeros(nFracCell*grid_size^2,1);

% Loop over fracture cells
for iCell = 1:nFracCell
    
    % Grab the four vertices of the cell
    verts = fracture.vertex.coords(fracture.elements.vertices(iCell,:),:);
    
    % Break fracture cell into grid and return centroids of sub-cells
    subGridCentroids = rectangleGridCentroids(verts,grid_size);

    % Place in output array
    test_points(1+(iCell-1)*grid_size^2:iCell*grid_size^2,:) = subGridCentroids;
    frac_el(1+(iCell-1)*grid_size^2:iCell*grid_size^2,1) = iCell*ones(grid_size^2,1);

end

end

function centroids = rectangleGridCentroids(verts,N)
% A,B,C,D are 1x3 vectors (rectangle vertices)
% N is grid resolution
% centroids is NxNx3 array

A = verts(1,:);
B = verts(2,:);
C = verts(3,:);
D = verts(4,:);

A = A(:)'; B = B(:)'; D = D(:)';

u = B - A;
v = D - A;

centroids = zeros(N*N,3);

for i = 0:N-1
    for j = 0:N-1
        s = (i + 0.5)/N;
        t = (j + 0.5)/N;

        P = A + s*u + t*v;
        centroids((i+1)+j*N,:) = P;
    end
end
end
