function [normals] = calc_normals_from_points(p1,p2,p3)
% compute_normals based on points p1 p2 p3
%   normals is an unitary vector
vec1 = p2-p1;
vec2 = p3-p1;
normals = cross(vec1,vec2,2);
normals = normals ./ vecnorm(normals,2,2);
end


