%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel
%TU Delft
%
% Find which rock volumes contain which test points, with bounded memory.
% Points are grouped by the fracture cell they were sampled in, and each
% group is only paired with the rock volumes whose bounding box overlaps
% the bounding box of the group. The (volume, point) pairs are then
% processed in batches of at most MaxPairs, instead of building the full
% (volumes x points) arrays.
%
% Input:
%   TestPoints:     (N_points x 3) coordinates of all test points
%   PointIDs:       indices (into TestPoints) of the points to test
%   FracEl:         (N_points x 1) fracture cell of every test point
%   VolumeIDs:      indices of the candidate rock volumes
%   MinVec, MaxVec: (N_rock x 3) rock volume bounding boxes
%   Tol:            bounding box tolerance
%   RockMesh, RockNormals, RockFaceCent: rock faces for is_in_vol3d
%
% Output:
%   IsIn:           (N_hits x 2) [rock volume ID, point ID]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function IsIn = find_points_in_volumes(TestPoints, PointIDs, FracEl, ...
                                       VolumeIDs, MinVec, MaxVec, Tol, ...
                                       RockMesh, RockNormals, RockFaceCent)

MaxPairs = 4e6;
MaxOverlap = 5e7;

PointIDs = PointIDs(:);
VolumeIDs = VolumeIDs(:);
IsIn = zeros(0, 2, 'like', PointIDs);
if isempty(PointIDs) || isempty(VolumeIDs)
    return
end

% Group the points by fracture cell
[Group, Order] = sort(double(FracEl(PointIDs)));
PointIDs = PointIDs(Order);
[~, GroupFirst] = unique(Group, 'first');
GroupCount = diff([GroupFirst; numel(Group) + 1]);
GroupID = repelem((1:numel(GroupFirst))', GroupCount);

% Bounding box of the points of each group
Coords = double(TestPoints(PointIDs, :));
GroupMin = zeros(numel(GroupFirst), 3);
GroupMax = zeros(numel(GroupFirst), 3);
for iDim = 1:3
    GroupMin(:, iDim) = accumarray(GroupID, Coords(:, iDim), [], @min);
    GroupMax(:, iDim) = accumarray(GroupID, Coords(:, iDim), [], @max);
end
clear Coords

% Candidate (volume, group) pairs: overlapping bounding boxes, tested
% for chunks of groups to bound the size of the overlap matrix
VolMin = double(MinVec(VolumeIDs, :)) - Tol;
VolMax = double(MaxVec(VolumeIDs, :)) + Tol;
nGroup = numel(GroupFirst);
GroupChunk = max(1, floor(MaxOverlap / numel(VolumeIDs)));
CandVol = cell(0, 1);
CandGroup = cell(0, 1);
for iFirst = 1:GroupChunk:nGroup
    Groups = iFirst:min(nGroup, iFirst + GroupChunk - 1);
    Overlap = true(numel(VolumeIDs), numel(Groups));
    for iDim = 1:3
        Overlap = Overlap & ...
            (VolMin(:, iDim) <= GroupMax(Groups, iDim)') & ...
            (VolMax(:, iDim) >= GroupMin(Groups, iDim)');
    end
    [LocalVol, LocalGroup] = find(Overlap);
    CandVol{end+1, 1} = LocalVol; %#ok<AGROW>
    CandGroup{end+1, 1} = reshape(Groups(LocalGroup), [], 1); %#ok<AGROW>
end
CandVol = vertcat(CandVol{:});
CandGroup = vertcat(CandGroup{:});
if isempty(CandVol)
    return
end

% Process the candidate pairs in batches of about MaxPairs point tests
PairCount = GroupCount(CandGroup);
Batch = floor(cumsum(PairCount) / MaxPairs);
Batch = Batch - Batch(1) + 1;
IsIn = repmat({zeros(0, 2, 'like', PointIDs)}, max(Batch), 1);
for iBatch = 1:max(Batch)
    InBatch = Batch == iBatch;
    if ~any(InBatch)
        continue
    end
    BVol = VolumeIDs(CandVol(InBatch));
    BGroup = CandGroup(InBatch);
    BCount = GroupCount(BGroup);

    % Expand every (volume, group) into (volume, point) pairs
    ii = repelem(BVol, BCount);
    Offset = (1:sum(BCount))' - repelem(cumsum([0; BCount(1:end-1)]), BCount);
    jj = PointIDs(repelem(GroupFirst(BGroup), BCount) + Offset - 1);

    % Point inside the volume bounding box
    Points = TestPoints(jj, :);
    Ref = all(Points - (MinVec(ii, :) - Tol) >= 0, 2) & ...
          all((MaxVec(ii, :) + Tol) - Points >= 0, 2);
    ii = ii(Ref);
    jj = jj(Ref);

    % Point inside the volume
    InRef = is_in_vol3d(TestPoints(jj, :), RockMesh.element.faces(ii, :), ...
                        RockNormals, RockFaceCent, ...
                        RockMesh.element.centroid(ii, :));
    IsIn{iBatch} = [ii(InRef) jj(InRef)];
end
IsIn = cell2mat(IsIn);

end
