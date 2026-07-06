function [out] = is_split(p1,p2, coords, normal)
tol = 10^-13;
vec = p2 -p1;
pc = 0.5*(p2+p1);
%d = vecnorm(vec,1,2);
plane_normal = cross(vec,normal);
plane_normal = repelem(plane_normal/norm(plane_normal),size(coords,1),1);
vec_test = coords - repelem(pc,size(coords,1),1);
dprod = dot(plane_normal,vec_test,2);
dprod(abs(dprod) < tol) = 0;
out = ones(size(dprod));
out(dprod>0) = 4;
out(dprod<0) = 8;
out(dprod==0) = 2;
end