%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel;
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p1xyz,p2xyz,p3xyz,p4xyz] = create_3d_boundbox(coords, normal,N)
p1 = min(coords);
p2 = max(coords);
p_avg = 0.5*(p1+p2);
[ncoords] = proj_coplanar_coords_on_xy(p_avg,...
    normal,[p1;p2], 1);
pz = mean(ncoords(:,end));
p_xy = ncoords(:,1:2);
[p1xy,p2xy,p3xy,p4xy] = create_subbboxes(p_xy(1,:), p_xy(2,:), N);
p1xy = [p1xy pz*ones(size(p1xy,1),1)];
p2xy = [p2xy pz*ones(size(p1xy,1),1)];
p3xy = [p3xy pz*ones(size(p1xy,1),1)];
p4xy = [p4xy pz*ones(size(p1xy,1),1)];

[p1xyz] = proj_coplanar_coords_on_xy(p_avg,...
    normal,p1xy, 0);
[p2xyz] = proj_coplanar_coords_on_xy(p_avg,...
    normal,p2xy, 0);
[p3xyz] = proj_coplanar_coords_on_xy(p_avg,...
    normal,p3xy, 0);
[p4xyz] = proj_coplanar_coords_on_xy(p_avg,...
    normal,p4xy, 0);
end
