%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel;
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [min_vec, max_vec] = generate_boundbox(adj,coords, useHalf)
min_vec = coords(adj(:,1),:);
max_vec = coords(adj(:,1),:);
if useHalf
    min_vec = half(min_vec);
    max_vec = half(max_vec);
end
size_c = size(adj,2);
for index = 1:size_c
    ref = adj(:,index)  ~= 0;
    local_min = min_vec(ref,:);
    local_max = max_vec(ref,:);
    index_adj_coord = coords(adj(ref,index),:);
    lmin_ref = index_adj_coord  - local_min  <= 0;
    local_min(lmin_ref) = index_adj_coord(lmin_ref);
    lmax_ref = index_adj_coord  - local_max  >= 0;
    local_min(lmin_ref) = index_adj_coord(lmin_ref);
    local_max(lmax_ref) = index_adj_coord(lmax_ref);
    min_vec(ref,:) = local_min;
    max_vec(ref,:) = local_max;
end
end
