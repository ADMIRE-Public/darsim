%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel;
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function boundbox = in_planar_bounding_box(coords,plane_normal, all_points)
% [~,fp2,fp3,fp4,fp5,fp6,fp7,~] = generate_bound_box_vertices(fp1,fp8);
% all_points = [fp1; fp2; fp3; fp4; fp5; fp6; fp7; fp8];
p1 = coords;
p2 = circshift(coords,-1,1);
p3 = p1 + repmat(plane_normal,size(p1,1),1);
p_edge = (1/3).*(p1 + p2 + p3);
pc = mean([p1;p3]);
[normals] = calc_normals_from_points(p1,p2,p3);
boundbox = false(size(all_points,1), size(normals,1));
for index = 1:size(normals,1)
    ref = semi_plane_check(all_points,normals(index,:),p_edge(index,:),pc);
    boundbox(:,index) = ref;
end
boundbox = all(boundbox,2);
%boundbox = any(reshape(all(boundbox,2),[],8),2);
end
