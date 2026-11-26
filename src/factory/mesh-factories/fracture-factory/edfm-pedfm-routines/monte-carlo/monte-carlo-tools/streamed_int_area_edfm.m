%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel
%TU Delft
%
% EDFM intersection areas of one fracture with streamed test points: the
% fracture cells are processed in batches; for each batch the test points
% are generated, located in the rock volumes and reduced to weighted hit
% counts per (rock cell, fracture cell), then discarded. Peak memory is
% set by the batch size, not by the fracture or mesh size.
%
% A point found in k rock cells (on a shared face) weighs 1/k; the area of
% a fracture cell is normalized by all points generated in it, so points
% outside the rock domain carry their share of area out of the domain.
%
% Output:
%   RockFracID: (N x 2) [rock cell ID, fracture cell ID]
%   IntArea:    (N x 1) one-sided intersection area
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RockFracID, IntArea] = streamed_int_area_edfm(FracMesh, ...
            PlaneGeo, TestVolumes, MinVec, MaxVec, Tol, RockMesh, ...
            RockNormals, RockFaceCent, Settings)

MaxPoints = Settings.batch_size;
n = Settings.num_area_points;
nEl = size(FracMesh.elements.vertices, 1);
BatchSize = max(1, floor(MaxPoints / n^2));
VolumeIDs = find_uint(TestVolumes);

PointCount = zeros(nEl, 1);
Keys = cell(0, 1);
Weights = cell(0, 1);
for First = 1:BatchSize:nEl
    Elements = (First:min(nEl, First + BatchSize - 1))';
    [Points, PointEl] = generate_testpoint_parameters(FracMesh, ...
        PlaneGeo.normal, n, Settings.useHalf, Elements);
    PointCount = PointCount + accumarray(double(PointEl), 1, [nEl 1]);
    IsIn = find_points_in_volumes(Points, uint32(1:size(Points, 1))', ...
        PointEl, VolumeIDs, MinVec, MaxVec, Tol, RockMesh, ...
        RockNormals, RockFaceCent, Settings.batch_size);
    if isempty(IsIn)
        continue
    end
    Mult = accumarray(double(IsIn(:, 2)), 1, [size(Points, 1) 1]);
    [Pairs, ~, iPair] = unique([double(IsIn(:, 1)) ...
                                double(PointEl(IsIn(:, 2)))], "rows");
    Keys{end+1, 1} = Pairs; %#ok<AGROW>
    Weights{end+1, 1} = accumarray(iPair, 1 ./ Mult(IsIn(:, 2))); %#ok<AGROW>
end

if isempty(Keys)
    RockFracID = zeros(0, 2);
    IntArea = zeros(0, 1);
    return
end
% batches hold disjoint fracture cells, so pairs never repeat
RockFracID = vertcat(Keys{:});
[RockFracID, Order] = sortrows(RockFracID);
W = vertcat(Weights{:});
W = W(Order);
IntArea = W ./ PointCount(RockFracID(:, 2)) .* ...
          FracMesh.elements.area(RockFracID(:, 2));

end
