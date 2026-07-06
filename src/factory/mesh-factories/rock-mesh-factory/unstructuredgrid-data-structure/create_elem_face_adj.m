function [element_face_adj,element_face_adj_bool] = ...
             create_elem_face_adj(inner_up, inner_down, bound_up, el_size)
%create_elem_face_adj create the face_element adjacency matrix
%% slighly modifed than the tetrahedron version 
num_inner = size(inner_up,1);
num_bound = size(bound_up,1);
ind_inner_vec = [1:num_inner]';
ind_bound_vec = [1:num_bound]';
ind_bound_vec = ind_bound_vec + num_inner;
[face_elem_adj_bool] = sparse([ind_inner_vec; ind_inner_vec; ...
                     ind_bound_vec],[inner_up; inner_down;bound_up], true);
element_face_adj_bool = face_elem_adj_bool';
[element_face_adj,~] = find(element_face_adj_bool');
[element_face_adj] = reshape(element_face_adj,el_size,[])';
end