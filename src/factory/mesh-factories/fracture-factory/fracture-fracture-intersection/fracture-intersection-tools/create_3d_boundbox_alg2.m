%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel;
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p1xyz,p2xyz,p3xyz,p4xyz] = ...
    create_3d_boundbox_alg2(coords, normal,N)
p1xyz = [];
p2xyz = [];
p3xyz = [];
p4xyz = [];
pref = coords(1,:);
vec1 = (coords(2,:) - coords(1,:))/N;
vec2 = (coords(4,:) - coords(1,:))/N;
node_transfer = @(i,j)  pref + i*vec1 + j*vec2;
for ii = 0:N-1
    for jj =0:N-1
        coords_ii = [ii ii+1 ii+1 ii];
        coords_jj = [jj jj jj+1 jj+1];        
        p1xyz = [p1xyz; node_transfer(coords_ii(1),coords_jj(1))];
        p2xyz = [p2xyz; node_transfer(coords_ii(2),coords_jj(2))];
        p3xyz = [p3xyz; node_transfer(coords_ii(3),coords_jj(3))];
        p4xyz = [p4xyz; node_transfer(coords_ii(4),coords_jj(4))];
    end
end

% p1 = min(coords);
% p2 = max(coords);
% p_avg = 0.5*(p1+p2);
% [ncoords] = proj_coplanar_coords_on_xy(p_avg,...
%     normal,[p1;p2], 1);
% pz = mean(ncoords(:,end));
% p_xy = ncoords(:,1:2);
% [p1xy,p2xy,p3xy,p4xy] = create_subbboxes(p_xy(1,:), p_xy(2,:), N);
% p1xy = [p1xy pz*ones(size(p1xy,1),1)];
% p2xy = [p2xy pz*ones(size(p1xy,1),1)];
% p3xy = [p3xy pz*ones(size(p1xy,1),1)];
% p4xy = [p4xy pz*ones(size(p1xy,1),1)];
% 
% [p1xyz] = proj_coplanar_coords_on_xy(p_avg,...
%     normal,p1xy, 0);
% [p2xyz] = proj_coplanar_coords_on_xy(p_avg,...
%     normal,p2xy, 0);
% [p3xyz] = proj_coplanar_coords_on_xy(p_avg,...
%     normal,p3xy, 0);
% [p4xyz] = proj_coplanar_coords_on_xy(p_avg,...
%     normal,p4xy, 0);
end
