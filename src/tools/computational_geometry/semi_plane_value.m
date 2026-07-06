function [val] = semi_plane_value(list_coord, normal, fcenter)
% 
%   returns the value of the points if points are in the same same plane than el_center
% n = size(list_coord,1);
dt = dot(normal, fcenter,2);
% flag = dot(normal,el_center) -  d;

%dt = ones(n,1)*d;
%rep_normal= repelem(normal,n,1);
val = dot(normal,list_coord,2) - dt;
% if flag  < 0 
%     ref = val < 0;
% else
%     ref = val > 0;
% end
end
