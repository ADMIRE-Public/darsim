function [area] = calc_area_convex_polygon2d(adj,coords)
% p1 = coords(adj(:,1),:);
% p2 = coords(adj(:,2),:);
% p3 = coords(adj(:,3),:);
% vec1 = p2 - p1;
% vec2 = p3 - p1;
% area = 0.5*(vecnorm(cross(vec1,vec2,2),2,2);
%%
area = zeros(size(adj,1),1);
num_col = size(adj,2);
p1 = coords(adj(:,1),:);
adj = adj(:,2:end);
for index = 0:(num_col-3) 
    cadj = circshift(adj,-index,2);
    ref = (cadj(:,1) ~= 0) & (cadj(:,2) ~= 0);
    p2 = coords(cadj(ref,1),:);
    p3 = coords(cadj(ref,2),:);
    area(ref) = area(ref) + calc_area_triang_using_coords(p1(ref,:),p2,p3);
end

end