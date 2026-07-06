function [dot_product] = find_dotproduct_frac_normal_faces_normal( ...
                       frac_normal, possible_projected_faces, all_normals)
frac_normal = repmat(frac_normal, size(possible_projected_faces,1),1);
dot_product = ones(size(possible_projected_faces))*10^16;
for index = 1:size(possible_projected_faces,2)
    ref = possible_projected_faces(:,index) ~= 0;
    dot_product(ref,index) = dot(frac_normal(ref,:), ...
                all_normals(possible_projected_faces(ref,index),:),2);
end
end