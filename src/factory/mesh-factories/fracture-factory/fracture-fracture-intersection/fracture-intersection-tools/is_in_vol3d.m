%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel;
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function is_in = is_in_vol3d(points, el_faces_adj, normals, ...
    normal_base_point, el_center)
tol = 10^-5;
is_in_bool = false(size(el_faces_adj));
for index = 1:size(el_faces_adj,2)
    ref = el_faces_adj(:,index) ~= 0;
    is_in_bool(~ref,index) = true;
    fadj = el_faces_adj(ref,index);
    value = semi_plane_value(points(ref,:), normals(fadj,:), ...
        normal_base_point(fadj,:) );
    value_center = semi_plane_value(el_center(ref,:), normals(fadj,:), ...
        normal_base_point(fadj,:) );    
    %value_total = value .*value_center;
    %lref_pos = value_total > 0;
    lref_pos = compare_signs(value,value_center,tol);
    ref(ref) = lref_pos;
    is_in_bool(ref,index) = true;
end
is_in = all(is_in_bool,2);
end
