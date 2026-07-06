function [face_value] = face_average(neighbours,values)
all_neigh_not_zero = all(neighbours ~=0,2);
num_lfaces = size(neighbours,1);
face_value = zeros(num_lfaces,size(values,2));
ref = neighbours(:,1) ~= 0;
face_value(ref,:) = face_value(ref,:) + values(neighbours(ref,1));
ref = neighbours(:,2) ~= 0;
face_value(ref,:) =  face_value(ref,:) + values(neighbours(ref,2));
face_value(all_neigh_not_zero) = face_value(all_neigh_not_zero)./2;
end