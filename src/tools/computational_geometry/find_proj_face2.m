function [pedfm_face, pedfm_volume, vol_ind] = ...
    find_proj_face2(vol_faces,face_ref, dist_point, vol_ind,rockmesh, rule)
face_ref = face_ref ~=0;
face_ref = face_ref & rule;
possible_faces = any(face_ref,2);
vol_faces = vol_faces(possible_faces,:);
face_ref = face_ref(possible_faces,:);
dist_point = dist_point(possible_faces,:);
%area = area(possible_faces,:);
vol_ind = vol_ind(possible_faces,:);
min_ref = Inf*ones(size(vol_faces,1),1);
%non_zeros = sum(face_ref,2);
for index = 1:size(dist_point,2)
    lref = face_ref(:,index); 
    min_ref(lref) = min(min_ref(lref), dist_point(lref,index),2);
end
face_ref = dist_point == min_ref;
pedfm_face = sum(vol_faces .* face_ref,2);
%pedfm_area = sum(area .* face_ref,2);
pedfm_volume = find_opposite_volume(vol_ind(:,1), pedfm_face, rockmesh);
end

