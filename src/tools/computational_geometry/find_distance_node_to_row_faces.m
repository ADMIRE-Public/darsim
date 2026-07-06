function [distance_matrix] = find_distance_node_to_row_faces(... 
                rock_frac_center, possible_projected_faces, all_fcenter)
distance_matrix = ones(size(possible_projected_faces)) *10^16;
for index = 1:size(possible_projected_faces,2)
    ref = possible_projected_faces(:,index) ~= 0;
    distance_matrix(ref,index) = vecnorm(rock_frac_center - ...
        all_fcenter(possible_projected_faces(ref,index),:),2,2);  
end
end