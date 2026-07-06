%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel;
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [int_rock_frac, pedfm_area, alpha_ind, alpha_area] = ...
   find_pdefm_projections_area(rock_frac_ind,int_area, is_in,...
                                    frac_el, test_points, rockmesh, eq)
%% finding pEDFM projection faces
all_normals = [rockmesh.face.inner.normal ; rockmesh.face.bound.normal];
all_normals = all_normals./vecnorm(all_normals,2,2);
all_face_neigh = [rockmesh.face.inner.upstream, ...
                  rockmesh.face.inner.downstream; ...
                  rockmesh.face.bound.upstream, ...
                  zeros(size(rockmesh.face.bound.upstream))];
all_fcenter = [rockmesh.face.inner.centroid ; ...
    rockmesh.face.bound.centroid];
all_fnodes = [rockmesh.face.inner.vertices ; rockmesh.face.bound.vertices];
fracture_plane_normal = [eq.A eq.B eq.C];

%% sorting is_in so unique will to ensure its order goes as rock_frac_ind
% sorting is wrong therefore removed
[~,ref] = sortrows([is_in(:,1), frac_el(is_in(:,2))]);
is_in = is_in(ref,:);
%frac_el = frac_el(ref);
%test_points = test_points(ref,:);
%% finding the indices of entry points in the same EDFM group
% group1 = Rock Cell I x Frac Cell J 
% group2 = Rock Cell I+1 x Frac Cell J
rock_cell_ind =  is_in(:,1);
test_node_ind = is_in(:,2);
frac_cell_ind = frac_el(test_node_ind);

[~,~,rock_frac_group] = unique([rock_cell_ind frac_cell_ind],"rows");

%% finding center area of  edfm GROUP
num_node_entries = accumarray(rock_frac_group,1);
gcenter = [accumarray(rock_frac_group,test_points(test_node_ind,1)), ...
           accumarray(rock_frac_group,test_points(test_node_ind,2)), ...
           accumarray(rock_frac_group,test_points(test_node_ind,3))];
rock_frac_center = gcenter./ num_node_entries;

% normal_frac_xy = dot(fracture_plane_normal,[0 0 1])./vecnorm([0 0 1]) * ...
%                                                              [0 0 1];
% normal_frac_xz = dot(fracture_plane_normal,[0 1 0])./vecnorm([0 1 0]) * ...
%                                                              [0 1 0];
% normal_frac_yz = dot(fracture_plane_normal,[1 1 0])./vecnorm([1 0 0]) * ...
%                                                              [1 0 0];
% % 
% w_xy = vecnorm(normal_frac_xy);
% w_xz = vecnorm(normal_frac_xz);
% w_yz = vecnorm(normal_frac_yz);
% w_sum = w_xy + w_xz + w_yz;
% w_xy = w_xy/w_sum;
% w_xz = w_xz/w_sum;
% w_yz = w_yz/w_sum;


%[~,frac_group] = sortrows(rock_frac_ind);
% [~,frac_group_big]  = ...The easiest way to install the Ubuntu packages is to first add the opm personal package archive (ppa). In order to do that we need to install the apt-add-repository command.


%             sortrows([is_in(:,1) frac_el(is_in(:,2))]);
% is_in = is_in(frac_group_big,:);
% [~,~,frac_group_big]  = unique([is_in(:,1) frac_el(is_in(:,2))],"rows");
% [~,~,fexd] = unique(rock_frac_ind(:,1));
% 
crossed_volume = unique(rock_frac_ind(:,1));
crossed_area = accumarray(rock_frac_ind(:,1), int_area);
crossed_area = crossed_area(crossed_volume);

possible_projected_faces = rockmesh.element.faces(rock_frac_ind(:,1),:);

dcenter_possible_faces = find_distance_node_to_row_faces(... 
                rock_frac_center, possible_projected_faces, all_fcenter);

dotproduc_fracnormal_faces = find_dotproduct_frac_normal_faces_normal( ...
             fracture_plane_normal, possible_projected_faces, all_normals);

not_possible_faces = abs(dotproduc_fracnormal_faces - 0) < 10^-14;

dotproduc_xy_faces = abs(find_dotproduct_frac_normal_faces_normal( ...
             [0 0 1], possible_projected_faces, all_normals));
dotproduc_xy_faces(not_possible_faces) = 0;

face_xy = find_possible_face_based_distance(possible_projected_faces,  ...
                            dotproduc_xy_faces, dcenter_possible_faces);


dotproduc_xz_faces = abs(find_dotproduct_frac_normal_faces_normal( ...
             [0 1 0], possible_projected_faces, all_normals));
dotproduc_xz_faces(not_possible_faces) = 0;

face_xz = find_possible_face_based_distance(possible_projected_faces,  ...
                            dotproduc_xz_faces, dcenter_possible_faces);


dotproduc_yz_faces = abs(find_dotproduct_frac_normal_faces_normal( ...
             [1 0 0], possible_projected_faces, all_normals));
dotproduc_yz_faces(not_possible_faces) = 0;


face_yz = find_possible_face_based_distance(possible_projected_faces,  ...
                            dotproduc_yz_faces, dcenter_possible_faces);


ref1 = face_xy ~=0;
ref2 = face_xz ~=0;
ref3 = face_yz ~=0;
 
projected_faces = [double(rock_frac_ind(ref1,:)) face_xy(ref1) int_area(ref1) ; ...
                   double(rock_frac_ind(ref2,:)) face_xz(ref2) int_area(ref2) ;
                   double(rock_frac_ind(ref3,:)) face_yz(ref3) int_area(ref3) ];

[~,ref] = sortrows(projected_faces(:,1:3));
projected_faces = projected_faces(ref,:);
[~,~,frac_proj_int_group] = unique(projected_faces(:,1:2),"rows");

faces_list = unique(frac_proj_int_group)';
array = cell(length(faces_list),1);
for index = faces_list
    ref_points = rock_frac_group == index;
    ref_frac_group = frac_proj_int_group ==  index;
    local_faces = unique(projected_faces(ref_frac_group,3));
    num_points = sum(ref_points);
    list_points = repmat(test_points(ref_points,:), ...
                                                size(local_faces,1),1) ;
    local_normals = repelem(all_normals(local_faces,:),num_points,1);
    local_frac_groups = repelem([1:size(local_faces,1)]',num_points,1);
    local_fcenter = repelem(all_fcenter(local_faces,:),num_points,1);
    local_fnodes= repelem(all_fnodes(local_faces,:),num_points,1);  
    [projected_points] = project_point_on_face(local_normals, ...
                            local_fcenter, list_points);
    [ref_pdefm] = find_projected_nodes_in_faces(projected_points, ...
        local_normals, local_fcenter, local_fnodes, rockmesh);
    point_count_percentage = accumarray(local_frac_groups,ref_pdefm)./...
                                num_points;    
    local_areas = repmat(int_area(index),size(point_count_percentage,1),1);
    local_projected_areas_full = abs(local_areas .* dot(...
            repmat(fracture_plane_normal, size(local_faces,1),1),...
                                        all_normals(local_faces,:),2));
    local_pedfm_area = point_count_percentage .* ...
                                        local_projected_areas_full;

    d = find_dist_point_plane_speed(projected_points, ...
                                fracture_plane_normal, list_points);
    local_num_points = accumarray(local_frac_groups, ref_pdefm);
    local_d_avg = accumarray(local_frac_groups, d.*ref_pdefm) ./ ...
                                    local_num_points;
    local_rock_vol = repmat(rock_frac_ind(index,1), ...
                    size(point_count_percentage,1),1);
    local_frac_vol = repmat(rock_frac_ind(index,2), ...
                        size(point_count_percentage,1),1);   

    local_neigh = all_face_neigh(local_faces,:);
    local_pedfm_vol =  sum(~(local_neigh == local_rock_vol) .* ...
                                            local_neigh,2);
    
    %% removing pedfm influence on boundary faces
    ref_boundary = local_pedfm_vol == 0;
    % array{index} = [double(local_pedfm_vol), double(local_frac_vol), ... 
    %                 double(local_rock_vol), local_pedfm_area,  ...
    %                 local_d_avg, local_faces];
    array{index} = [double(local_pedfm_vol(~ref_boundary,:)), ...
                    double(local_frac_vol(~ref_boundary,:)), ... 
                    double(local_rock_vol(~ref_boundary,:)), ...
                    local_pedfm_area(~ref_boundary,:),  ...
                    local_d_avg(~ref_boundary,:), ...
                    local_faces(~ref_boundary,:)];
end
array = vertcat(array{:});
%% merging pdfm
[rock_frac_ind,ref_bb,pedfm_group] = unique(array(:,1:2),"rows");
pedfm_area = accumarray(pedfm_group,array(:,4) );
int_rock_frac = rock_frac_ind;
pedfm_face = array(ref_bb,end);

[alpha_ind,~,alpha_group] = unique(array(:,end));
alpha_area = accumarray(alpha_group,array(:,4));


%% merging groups, summing areeas
% [int_rock_frac,~,ref] = unique(array(:,1:2),"rows");
% pedfm_area = accumarray(ref, array(:,4));
% % int_percentage = array(:,4)./ pedfm_area(ref);
% % int_avg_d = array(:,4) .*int_percentage;
% % int_avg_d = accumarray(ref,int_avg_d);
% pedfm_faces = array(ref,end);



% proj_faces_area = proj_faces_area .*proj_faces_area;
% dist_point_faces(dot_prod_result == 0 ) = Inf;
% 
% [projected_point] = project_point_on_face( ...
%       all_normals(faces_to_points,:), all_fcenter(faces_to_points,:),...
%       list_points);
% 
% [ref_pdefm] = find_projected_nodes_in_faces(projected_point,...
%     all_normals(faces_to_points,:), all_fcenter(faces_to_points,:), ...
%     all_fnodes(faces_to_points,:), rockmesh); 
% 
% face_projected_ind = full(trans_matr(rock_frac_ind(:,1))); 
% [proj_area] = find_all_projected_area(face_projected_ind, int_area, ...
%     all_normals, eq);
 %debug_cornerfrac
 % [projected_faces] = find_projected_face(rockmesh, eq, crossed_volume);

% pedfm_vol_face = [crossed_volume, projected_faces];
% [int_rock_frac_uniq, int_pedfm_area, int_rock_face] = ...
%                 calculate_int_area_pedfm(is_in, frac_el, ref_pdefm, ...
%              proj_area, faces_to_points, pedfm_vol_face, untitled27rockmesh,rock_frac_ind );

%% improved the choice of projected pedfm face
% now it tries to select the face with the least dot product and distance 
% [crossed_faces] = find_crossed_faces(crossed_volume,rockmeshrock_vol_frac,eq);
%% if issues in calculate_int_area_pedfm verify if cross_volumes always 
% is associated with always a different projected_faces
% [projected_faces] = find_projected_face(rockmesh, eq, crossed_volume);


%% finding groups
% 
% [~,~,frac_int_group] = unique([is_in(:,1) frac_el(is_in(:,end))],"rows");
% 
% %% crossed volumes are comprised of all CV crossed by the Fracture Plane
% crossed_volume = unique(rock_frac_ind(:,1));
% [possible_projected_faces, face_x,face_y, face_z] = ...
%     find_bad_projected_face_3d(rockmesh, eq, crossed_volume);
% 
% is_in_possible_faces = rockmesh.element.faces(crossed_volume(:,1),:);
% points_in_faces = zeros(size(crossed_volume,1), ...
%                 size(is_in_possible_faces,2));
% dist_point_faces = zeros(size(crossed_volume,1), ...
%             array    size(is_in_possible_faces,2));
% proj_faces_area = zeros(size(crossed_volume,1), ...
%                 size(is_in_possible_faces,2));
% 
% for index = 1:size(rockmesh.element.faces,2)
%     ref = is_in_possible_faces(:,index) ~= 0;
%     [projected_points] = project_point_on_face( ...
%       all_normals(is_in_possible_faces(ref,index),:), ...
%       all_fcenter(is_in_possible_faces(ref,index),:),...
%       test_points(is_in(ref,2),:));
%     [ref_pdefm] = find_projected_nodes_in_faces(projected_points,...
%          all_normals(is_in_possible_faces(ref,index),:), ...
%          all_fcenter(is_in_possible_faces(ref,index),:), ...
%          all_fnodes(is_in_possible_faces(ref,index),:), rockmesh); 
%     ref_tmp = ref;
%     ref_tmp(ref_tmp) = ref_pdefm;   
%     h1 = accumarray(frac_int_group,1);
%     h2 = accumarray(frac_int_group,ref_tmp);
%     points_in_faces(:,index) = h2./h1;
%     d = find_dist_point_plane(test_points(is_in(ref_tmp,2),:), ...
%        all_normals(is_in_possible_faces(ref_tmp,index),:), ...
%        all_fcenter(is_in_possible_faces(ref_tmp,index),:));
%     d_avg = accumarray(frac_int_group(ref_tmp),d);
%     group_ind = unique(frac_int_group(ref_tmp));
%     dist_point_faces(group_ind,index) = d_avg;
% end
% 
% for index = 1:size(rockmesh.element.faces,2)
% local_faces = rockmesh.element.faces(rock_frac_ind(:,1),index);
% ref = local_faces ~= 0;
% [total_new_area] = find_all_projected_area(local_faces(ref), ...
%             int_area, all_normals, eq);
% proj_faces_area(ref,index) = total_new_area;
% end
% 
% %% correcting area to remove elements whose projection are not in the faces
% proj_faces_area  = points_in_faces .*proj_faces_area;
% 
% dist_point_faces(possible_projected_faces == 0 ) = Inf;
% min_dif = max(proj_faces_area,[],2) -  min(proj_faces_area,[],2); 
% min_dif = percentage .*min_dif;
% proj_faces_area(proj_faces_area < min_dif) = 0; 
% dist_point_faces(proj_faces_area == 0) = Inf;
% ref = dist_point_faces == min(dist_point_faces,[],2);
% if sum(ref,2) > 1
%     error("fix here");
% end

% ref = dist_point_faces == min(dist_point_faces,[],2);
% pedfm_facesX = face_x
% [pedfm_faceX,pedfm_areaX, pedfm_volumeX, vol_indX] = ...
%         find_proj_face(rockmesh.element.faces(rock_frac_ind(:,1),:),...
%         face_x, dist_point_faces, proj_faces_area, ...
%         rock_frac_ind, rockmesh);
% 
% [pedfm_faceY,pedfm_areaY, pedfm_volumeY, vol_indY] = ...
%         find_proj_face(rockmesh.element.faces(rock_frac_ind(:,1),:),...
%         face_y, dist_point_faces, proj_faces_area, ...
%         rock_frac_ind, rockmesh);
% 
% [pedfm_faceZ,pedfmfa_areaZ, pedfm_volumeZ, vol_indZ] = ...
%         find_proj_face(rockmesh.element.faces(rock_frac_ind(:,1),:),...
%         face_z, dist_point_faces, proj_faces_area, ...face_center
%         rock_frac_ind, rockmesh);
% 
% pedfm_faces = [pedfm_faceX; pedfm_faceY; pedfm_faceZ];
% pedfm_area =  [pedfm_areaX; pedfm_areaY; pedfm_areaZ];
% pedfm_volume = [pedfm_volumeX; pedfm_volumeY; pedfm_volumeZ];
% int_rock_frac = [vol_indX; vol_indY; vol_indZ];
% int_rock_frac(:,1) = pedfm_volume;
% 
% pedfm_faces = sum(rockmesh.element.faces(rock_frac_ind(:,1),:) .*ref,2);
% pedfm_area = sum(proj_faces_area .*ref,2);
% [R_faces, R_facesX, R_facesY,R_facesZ]
% pedfm_volumes = find_opposite_volume(rock_frac_ind(:,1), ...
%                                              pedfm_faces, rockmesh);
% int_rock_frac_uniq = rock_frac_ind;
% int_rock_frac_uniq(:,1) = pedfm_volumes;

%% improved the choice of projected pedfm face
% now it tries to select the face with the least dot product and distance 
% [crossed_faces] = find_crossed_faces(crossed_volume,rockmesh,eq);
%% if issues in calculate_int_area_pedfm verify if cross_volumes always 
% % is associated with always a different projected_faces
%% projected faces finds the face of the pEDFM to be projected 
% [projected_faces] = find_projected_face(rockmesh, eq, crossed_volume);ar



% %% structure that associated volume and faces
% pedfm_vol_face = [crossed_volume, projected_faces];
% 
% trans_matr = sparse(double(crossed_volume),ones(size(crossed_volume)), ...
%                         projected_faces);
% faces_to_pointsfa = full(trans_matr(is_in(:,1)));
% %% faces_to points consists in all the pEDFM projected faces associated
% % with the volumes in is_in(:,1). 
% list_points = test_points(is_in(:,2),:);
%% list_points is a list points associated with the refrence points in
% is_in(:,2). if no points were removed the size of test_points would be 
% the same as is_in number of lines. also, to find the fracture element
% one has to use frac_el(is_in(:,))

%% project the points keeping the order and associated in is_in
% % points are projected on the associated faces
% [projected_point] = project_point_on_face( ...
%       all_normals(faces_to_points,:), all_fcenter(faces_to_points,:),...
%       list_points);
% %% verify if the projection points are inside the face, otherwise
% % they are excluded
% [ref_pdefm] = find_projected_nodes_in_faces(projected_point,...
%     all_normals(faces_to_points,:), all_fcenter(faces_to_points,:), ...
%     all_fnodes(faces_to_points,:), rockmesh); 
% %% project the EDFM area into the projection face plane.
% face_projected_ind = full(trans_matr(rock_frac_ind(:,1))); 
% [proj_area] = find_all_projected_area(face_projected_ind, int_area, ...
%     all_normals, eq);
%  %debug_cornerfrac
 %% uses the projected pedfm points to estimate area on the projected
 % face plane that lies inside the projected face. it removes the influence
 % of points not inside the plane.


 % verify which points are inside the faces and uses the proj_area to
 % estimate projection that lies inside 
% [int_rock_frac_uniq, int_pedfm_area, int_rock_face] = ...
%                 calculate_int_area_pedfm(is_in, frac_el, ref_pdefm, ...
%    proj_area, faces_to_points, pedfm_vol_face, rockmesh,rock_frac_ind );

end
function [chosen_face] =  find_possible_face_based_distance(faces, dot_prod, ...
    distances)
chosen_face = zeros(size(faces,1),1);
dot_ref = (dot_prod ~=0);
num_poss_faces = sum(dot_ref,2);
num_face_set = unique(num_poss_faces);
possible_faces = dot_ref .* faces;
dist_possible_faces = distances .* dot_ref;
if sum(num_face_set) ~= 0
    for index = setdiff(num_face_set,0)'
        ref = num_poss_faces == index;
        local_possible_faces = possible_faces(ref,:)';
        local_dist_possible_faces = dist_possible_faces(ref,:)';
        local_dot_ref = dot_ref(ref,:)';
        faces_candidates = ...
            reshape(local_possible_faces(local_dot_ref),index,[])';
        distance_candidates = ...
            reshape(local_dist_possible_faces(local_dot_ref),index,[])';
        equal_ref = distance_candidates == min(distance_candidates,[],2);
        equal_ref = process_boolean_matrix(equal_ref);
        chosen_face(ref,1) =  sum(equal_ref .* faces_candidates,2);
    end
end
end

function out = process_boolean_matrix(A)
% Find the indices of the first true value in each row
if ~isempty(A)
    [~, col_indices] = max(A, [], 2);
    % Create a column vector to index each row
    row_indices = (1:size(A, 1)).';
    % Create a linear index using row and column indices
    linear_indices = sub2ind(size(A), row_indices, col_indices);
    % Create a logical mask with true values at linear indices
    mask = false(size(A));
    mask(linear_indices) = true;
    % Output matrix with only the first true value in each row
    out = A & mask;
    % Handle the case where there are no true values in a row
    no_true_row = all(A == 0, 2);
    out(no_true_row, :) = 0;
else
    out = A;
end
end



function d = find_dist_point_plane_speed(test_points, normal, points)
%% we dont need to divide by sqrt(A^2 + B^2 +  C^2) the norm of the normals
% is 1
normal = repmat(normal,size(test_points,1),1);
d = abs(dot(test_points,normal,2) - dot(normal,points,2));
end


% function [filter] = apply_local_filter(input , table)
% input = double(input(:,1));
% ind = double(table(:,1));
% value = double(table(:,2));
% oind = ones(size(ind));
% max_value = max(ind);
% transf_matrix = sparse(ind,oind,value,max_value,1);
% filter = full(transf_matrix(input));
% end
% function [filter] = apply_pdefm_filter(input, ind,value)
%     ind = double(ind);    
%     oind = ones(size(ind));
%     max_value = max(ind);
%     transf_matrix = sparse(ind,oind,value,max_value,1);
%     filter = full(transf_matrix(input));
% end

function [projected_point] = ...
            project_point_on_face(normals, plane_points, points)    
    V = points - plane_points;
    d = dot(normals,V,2);
    projected_point = points - d.*normals;
end
