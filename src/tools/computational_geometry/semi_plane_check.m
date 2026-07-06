function [ref] = semi_plane_check(list_coord, normal, fcenter, el_center)
% 
%   returns 1 if points are in the same same plane than el_center
tol = 10^-8;
n = size(list_coord,1);
d = dot(normal, fcenter);
flag = dot(normal,el_center) -  d;

dt = ones(n,1)*d;
rep_normal= repelem(normal,n,1);
val = dot(rep_normal,list_coord,2) - dt;

ref = compare_signs(val, flag*ones(size(val)),tol);
% if flag  < 0 
%     ref = val < 0;
% else
%     ref = val > 0;
% end
end
