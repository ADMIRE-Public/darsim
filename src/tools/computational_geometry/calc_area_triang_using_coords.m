function [out] = calc_area_triang_using_coords(p1,p2,p3)
x1 = p1(:,1);
y1 = p1(:,2);
z1 = p1(:,3);
x2 = p2(:,1);
y2 = p2(:,2);
z2 = p2(:,3);
x3 = p3(:,1);
y3 = p3(:,2);
z3 = p3(:,3);
out = (abs((x1 - x2).*(y1 - y3) - (x1 - x3).*(y1 - y2)).^2 + ...
    abs((x1 - x2).*(z1 - z3) - (x1 - x3).*(z1 - z2)).^2 + ...
    abs((y1 - y2).*(z1 - z3) - (y1 - y3).*(z1 - z2)).^2).^(1/2);
out = 0.5 * out;
end