function [p_intersect ] = point_plane_intersection(pv, vecv, pp, vecp)
% pv -> point in line, vecv-> vector in the direction of the line
% pp -> point in plane, vecp-> normal vector
pv = repmat(pv,size(pp,1),1);
vecv = repmat(vecv,size(pp,1),1);
t = dot(vecp,pp-pv,2) ./ dot(vecp,vecv,2);
p_intersect = pv + t.*vecv;
end