function [pedfm_face, pedfm_volume, vol_ind] = ...
        find_proj_face(vol_faces,face_ref, dist_point, vol_ind,rockmesh)
face_ref = face_ref ~=0;
possible_faces = any(face_ref,2);
vol_faces = vol_faces(possible_faces,:);
face_ref = face_ref(possible_faces,:);
dist_point = dist_point(possible_faces,:);
%area = area(possible_faces,:);
vol_ind = vol_ind(possible_faces,:);
min_ref = zeros(size(vol_faces,1),1);
non_zeros = sum(face_ref,2);
for index = unique(non_zeros)'
    ref = non_zeros == index;
    local_faces_ref = face_ref(ref,:);
    local_dist_values = reshape(dist_point(local_faces_ref),[],index);
    min_ref(ref) = min(local_dist_values,[],2);
end
face_ref = dist_point == min_ref;
pedfm_face = sum(vol_faces .* face_ref,2);
%pedfm_area = sum(area .* face_ref,2);
pedfm_volume = find_opposite_volume(vol_ind(:,1), pedfm_face, rockmesh);
end

