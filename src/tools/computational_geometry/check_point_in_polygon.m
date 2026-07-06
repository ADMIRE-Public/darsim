function [ref] = check_point_in_polygon(poly_coord, check_coord,varargin)
% returns true if it is inide polygon and false otherwise
if iscell(check_coord)
    test_coord = cell2mat(check_coord);
else
    test_coord = check_coord;
end
tol = 10^-10;
p_center = mean(poly_coord);
p1 = poly_coord;
p2 = circshift(p1,-1,1);
if isempty(varargin)
    p3 = repmat(p_center, size(p1,1),1);
    Atotal = sum(calc_area_triang_using_coords(p1,p2,p3));
else
    Atotal = cell2mat(varargin);
end
Atotal = repelem(Atotal,size(test_coord,1),1);
p3 = repelem(test_coord,size(p1,1),1);
p1 = repmat(p1,size(test_coord,1),1);
p2 = repmat(p2,size(test_coord,1),1);
group_ref = repelem([1:size(test_coord,1)]',size(poly_coord,1),1);
A = calc_area_triang_using_coords(p1,p2,p3);
ref = abs(accumarray(group_ref,A) - Atotal) < tol;
end
