%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel;
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [num_intersection,p1_3d, p2_3d,v1,v2] = ...
     calculate_vol_edge_intersection(p1, p2, normal,fracmesh)
coords = fracmesh.vertex.coords;
adj_list = [fracmesh.faces.inner.vertices; fracmesh.faces.bound.vertices];
                % [num_int, int_p1,int_p2]  = ...
                %                 calculate_edge_intersection( ...
                %                  fracture_p1, fracture_p2,...
                %                  meshes(frac_mesh).vertex.coords, ...
                %                  geometry(frac_mesh).normal,...
                %                  meshes(frac_mesh).faces.inner.vertices,...
                %                  meshes(frac_mesh).faces.bound.vertices);
% 
% n_arg = size(varargin,2);
% size_input = ...
%    reshape(cell2mat(cellfun(@size,varargin,'UniformOutput',false))',[],2)';
%adj_list = cell2mat(varargin(:));
%num_intersection = sparse(size(adj_list,1),1);

%% projecting on 2d-space
center = mean(coords);
%% converts coords to 2d space
                % [num_int, int_p1,int_p2]  = ...
                %                 calculate_edge_intersection( ...
                %                  fracture_p1, fracture_p2,...
                %                  meshes(frac_mesh).vertex.coords, ...
                %                  geometry(frac_mesh).normal,...
                %                  meshes(frac_mesh).faces.inner.vertices,...
                %                  meshes(frac_mesh).faces.bound.vertices);
[xy_coords] = proj_coplanar_coords_on_xy(center,normal,coords,1);
[xy_p1] = proj_coplanar_coords_on_xy(center,normal,p1,1);
[xy_p2] = proj_coplanar_coords_on_xy(center,normal,p2,1);

[node_1,node_2,num_intersection] = ...
    edge_analysis_in_2d(adj_list, xy_coords, xy_p1, xy_p2);
%% finding if first and second nodes of the frac intersection line 
%  are contained by any volume
[v1,v2]= find_in_volume_frac_line(xy_p1, xy_p2, xy_coords, ...
     fracmesh.elements.faces, adj_list,...
     fracmesh.elements.vertices,fracmesh.elements.centroid);
%% projecting information back 
node_1(num_intersection > 0,3) = center(:,3);
node_2(num_intersection > 1,3) = center(:,3);
p1_3d = sparse(size(node_1,1),3);
p2_3d = sparse(size(node_2,1),3);
ref = num_intersection > 0;
[p1_3d(ref,:)] = proj_coplanar_coords_on_xy(center,normal,node_1(ref,:),0);
ref = num_intersection > 0;
[p2_3d(ref,:)] = proj_coplanar_coords_on_xy(center,normal,node_2(ref,:),0);
end

function [v1,v2]= find_in_volume_frac_line(p1, p2, ...
                    coords,vol_edges,edge_nodes,vol_nodes, cent_volume)
tol = 10^-12;
%cent_faces = (coords(edge_nodes(:,1),:) + coords(edge_nodes(:,2),:))*.5;
[min_vec, max_vec] = generate_boundbox(vol_nodes,coords, false);

%% bounding box to avoid uncessary calculations 
ref_p1_x = (min_vec(:,1) - tol <= p1(:,1)) & ...
           (max_vec(:,1) + tol >= p1(:,1));
ref_p1_y = (min_vec(:,2) - tol <= p1(:,2)) & ...
           (max_vec(:,2) + tol >= p1(:,2));
ref_p1 = ref_p1_x & ref_p1_y;
ref_p2_x = (min_vec(:,1) - tol <= p2(:,1)) & ...
           (max_vec(:,1) + tol >= p2(:,1));
ref_p2_y = (min_vec(:,2) - tol <= p2(:,2)) & ...
           (max_vec(:,2) + tol >= p2(:,2));
ref_p2 = ref_p2_x & ref_p2_y;
% %ref_test_vol = ref_p1 | ref_p2;
% [faces_edges1] = create_matedgezero(vol_nodes(ref_p1,:));
% [faces_edges2] = create_matedgezero(vol_nodes(ref_p2,:));
% 
% local_cent_volume_p1 = cent_faces(ref_p1,1:2);
% local_cent_volume_p2 = cent_faces(ref_p2,1:2);
% 1

%% testing volume 1
v1 = [];
for lvol = find(ref_p1)'
    lnodes = vol_nodes(lvol,:);
    plg = coords(lnodes,:);
    [flag_in, flag_on] = inpolygon(p1(:,1),p1(:,2), plg(:,1),plg(:,2));
    if flag_in || flag_on
        v1 = [v1; lvol];
    end
end
%% testing volume 2
v2 = [];
for lvol = find(ref_p2)'
    lnodes = vol_nodes(lvol,:);
    plg = coords(lnodes,:);
    [flag_in, flag_on] = inpolygon(p2(:,1),p2(:,2), plg(:,1),plg(:,2));
    if flag_in || flag_on
        v2 = [v2; lvol];
    end
end
%%
% rp1 = repmat(p1(:,1:2),size(local_cent_volume_p1,1),1);
% rp2 = repmat(p2(:,1:2),size(local_cent_volume_p2,1),1);
% out_ref1 = false(size(vol_nodes));
% out_ref2 = false(size(vol_nodes));

% for index = 1:(size(faces_edges1,2)-1)  
%     edge1 = faces_edges1(:,index:index+1);
%     edge2 = faces_edges2(:,index:index+1);
%     ref1 = ~any(edge1 == 0,2);
%     ref2 = ~any(edge2 == 0,2);
%     if any(ref1)
%         c1_2d = coords(edge1(ref1,1),1:2);
%         c2_2d = coords(edge1(ref1,2),1:2);
%         lref = is_above_2d(c1_2d,c2_2d,rp1,local_cent_volume_p1);
%         lbref = ref_p1;
%         lbref(lbref) = ref1;   
%         out_ref1(lbref,index) = lref;
%     end
%     if any(ref2)
%         c1_2d = coords(edge2(ref2,1),1:2);
%         c2_2d = coords(edge2(ref2,2),1:2);
%         lref = is_above_2d(c1_2d,c2_2d,rp2,local_cent_volume_p2);
%         lbref = ref_p2;
%         lbref(lbref) = ref2;
%         out_ref2(lbref,index) = lref;
%     end    
% end
% out_ref1(vol_edges == 0) = true;
% out_ref1 = all(out_ref1,2);
% out_ref2(vol_edges == 0) = true;
% out_ref2 = all(out_ref2,2);
%1

end
function [int_p1,int_p2,num_int] = edge_analysis_in_2d(adj,coords,fp1,fp2)
tol = 10^-13;
int_p1 = sparse(size(adj,1),2);
int_p2 = sparse(size(adj,1),2);
num_int = sparse(size(adj,1),1);
p1 = coords(adj(:,1),1:2);
p2 = coords(adj(:,2),1:2);
fp1 = fp1(:,1:2);
fp2 = fp2(:,1:2);
rfp1 = repmat(fp1,size(p1,1),1);
rfp2 = repmat(fp2,size(p2,1),1);
%% defining equation frac int equation
% eq1: X1 = p1 + t1(p2-p1)
%line_int(p1,p2 -p1 , rfp1,rfp2-rfp1,p1,p2)
%% defining system of equation
% edge parametric equation
a = p1;
vec1 = p2-p1;
% fractures intersection line parametric equation
b = rfp1;
vec2 = rfp2-rfp1;

%% finding systems that are not parallel
vecn1 = abs(vec1./vecnorm(vec1,2,2));
vecn2 = abs(vec2./vecnorm(vec2,2,2));
ref_parallel = all(abs(vecn1 - vecn2) < 10^-13,2);
ref_cross = ~ref_parallel;

%% adding the reference of edges crossing intfrac
[t1,t2] = solve_system_parametric(vec1(ref_cross,:),vec2(ref_cross,:),...
    b(ref_cross,:) -a(ref_cross,:));
ref1 = (t1 >= (0 - tol)) & (t1 <= (1 + tol));
ref2 = (t2 >= (0 - tol)) & (t2 <= (1 + tol));
ref = ref1 & ref2;
ref_cross(ref_cross) = ref;
num_int(ref_cross) = 1;
int_p1(ref_cross,:) = p1(ref_cross,:) + ...
    t1(ref) .*(p2(ref_cross,:) - p1(ref_cross,:));

%% adding influence of parallel edges
[t_p1] = interpolate_t(p1(ref_parallel,:),b(ref_parallel,:), ...
                                vec2(ref_parallel,:));
ref_tp1 = abs(t_p1(:,1) - t_p1(:,2)) < tol;
t_p1 = mean(t_p1,2);
[t_p2] = interpolate_t(p2(ref_parallel,:),b(ref_parallel,:),...
                                vec1(ref_parallel,:));
ref_tp2 = abs(t_p2(:,1) - t_p2(:,2)) < tol;
t_p2 = mean(t_p2,2);
ref_frac_ege = ref_tp1 & ref_tp2;
[t_p1(ref_frac_ege), t_p2(ref_frac_ege)] = ...
    sort_and_bound_interval(t_p1(ref_frac_ege), t_p2(ref_frac_ege));
% removing one node elements 
% ref_one_node = abs(t_p1(ref_frac_ege) - t_p2(ref_frac_ege)) > tol;edge_analysis_in_2d
% ref_frac_ege(ref_frac_ege) = ref_one_node;
ref_parallel(ref_parallel) = ref_frac_ege;
num_int(ref_parallel) = 2;
int_p1(ref_parallel,:) = rfp1(ref_parallel,:) + ...
    t_p1(ref_frac_ege) .*(rfp2(ref_parallel,:) - rfp1(ref_parallel,:));
int_p2(ref_parallel,:) = p1(ref_parallel,:) + ...
    t_p2(ref_frac_ege) .*(rfp2(ref_parallel,:) - rfp1(ref_parallel,:));

%% adding the influence of volumes that have either the initial or end point
%  of fracture intersection line
%    ______
%   |      |
% --x---E  |
%   |______|
%
%



end

function [nt1, nt2] = sort_and_bound_interval(t1, t2)
    auxmat = [t1 t2];
    auxmat(auxmat>1) = 1;
    auxmat(auxmat<0) = 0;
    auxmat = sort(auxmat,2);
    nt1 = auxmat(:,1);
    nt2 = auxmat(:,2);
end
function [t1,t2] =  solve_system_parametric(m,n,q)
m1 = m(:,1);
m2 = m(:,2);
n1 = n(:,1);
n2 = n(:,2);
q1 = q(:,1);
q2 = q(:,2);
t1 = (n2.*q1)./(m1.*n2 - m2.*n1) - (n1.*q2)./(m1.*n2 - m2.*n1);
t2 = (m2.*q1)./(m1.*n2 - m2.*n1) - (m1.*q2)./(m1.*n2 - m2.*n1);
end


function [t] = interpolate_t(x,a,b)
t = (x-a)./b;
end




