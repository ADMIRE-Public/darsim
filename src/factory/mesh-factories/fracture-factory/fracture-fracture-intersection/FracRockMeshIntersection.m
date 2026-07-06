%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel;
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [frac_element_ind,rock_volumes_ind,intersect_points ] = ...
    FracRockMeshIntersection(fracmesh,rockmesh, plane_geo, useHalf)
%% useHalf if a flag, if true the routines will use half precision routines
% if available
%% intersection between each fracmesh with the rockmatirx
nfrac = max(size(fracmesh));
rock_normals =  [rockmesh.face.inner.normal;...
    rockmesh.face.bound.normal];
rock_face_cent =  [rockmesh.face.inner.centroid;...
    rockmesh.face.bound.centroid];
for index = 1:nfrac
    tol = 10^-13;
    [test_points, frac_el] =...
        generate_testpoint_parameters(fracmesh(index), ...
        plane_geo(index).normal,60, useHalf);

    [test_volumes] = identify_plane_crossed_volumes(...
        rockmesh.element.vertices,...
        rockmesh.vertex.coord, plane_geo(index).eq);
    % check which cells of the adj matrix of the rock matrix have points
    % in the boundary box of the the tested points
    [min_vec, max_vec ] = ...
        generate_boundbox(rockmesh.element.vertices(test_volumes,:), ...
        rockmesh.vertex.coord, useHalf);
    in_frac = in_fracture_bounding_box(plane_geo(index), min_vec,max_vec);
    min_vec = min_vec(in_frac,:);
    max_vec = max_vec(in_frac,:);
    test_volumes(test_volumes) = in_frac;
    abs_min = min(min_vec);
    abs_max = max(max_vec);
    bound_points_ref = in_between(test_points,abs_min,abs_max,tol);
    % test_faces = generate_test_faces(rockmesh, test_volumes);
    sign_matrix = generate_test_volume_sign_matrix(rockmesh, test_volumes);
    % generate_test_volume_sign_matrix(rockmesh, test_volumes)

    tic, [ref_ij] = is_in_boundbox_hierarchically(test_points(...
        bound_points_ref,:), min_vec, max_vec,rockmesh.element.faces( ...
        test_volumes,:),rock_normals, rock_face_cent, sign_matrix); toc
    frac_element_ind = frac_el(bound_points_ref);
    frac_element_ind = frac_element_ind(ref_ij(:,1));
    intersect_points = test_points(bound_points_ref,:);
    intersect_points = intersect_points(ref_ij(:,1),:);
    rock_volumes_ind = find(test_volumes);
    rock_volumes_ind = rock_volumes_ind(ref_ij(:,2));
    1
end
end

function sign_matrix = generate_test_volume_sign_matrix(rockmesh, ...
    test_volumes)
sign_matrix = ...
    sparse_to_zero_ind_mat(rockmesh.face.summation(test_volumes,:));
end




function [is_in] = is_in_boundbox_hierarchically(points, min_vec, ...
    max_vec, rock_el_faces, rock_normals, rock_points, sign_matrix)
%% is_in_boundbox_hierarchically checks if points are in bound box
% is_in_boundbox_hierarchically uses splits the points in n x n x n
% bounding boxes.  if points is contained in a given bounding box as
% well as either one of min_vec or max_vec, these points are checked
% against the respective min_vec and max_vec. This allows for large
% to skip large amount of verifications decreasing memory and speed
% constrains.
% INPUT:
%       points   -> query points
%       frac_vol -> fracture indices associated with these points
%       min_vec and max_vec -> the boundary boxes of the volumes
% OUTPUT:min
%
%
% COMMENTS:
% the comparison between min_vec max_vec and bdiag needs to take in
% consideration not only the main diag but also the other nodes
%
% n = floor(0.0082*((size(points,1)* size(min_vec,1))^(1/3)));
% n = max(n,9);
% n = min(n,11);
n = 6;
% disp(n)
%tol = 0*eps(points(1,1));
tol = 10^-13;
[bdiag1, bdiag2]  = divide_box(min(min_vec),max(max_vec), n);
is_in = cell(size(bdiag1,1),1);
ind_points = uint32([1:size(points,1)]');
int_el_box = uint32([1:size(max_vec,1)]');
bcount = 1;
qq = [];
for index = 1:size(bdiag1,1)
    bound_box_local_min = bdiag1(index,:);
    bound_box_local_max = bdiag2(index,:);


    points_ref = ...
        in_between(points,bound_box_local_min, bound_box_local_max,tol);
    [~,p2,p3,p4,p5,p6,p7,~] = generate_bound_box_vertices(min_vec, ...
        max_vec);
    local_bound_ref = ...
        in_between(min_vec,bound_box_local_min,bound_box_local_max,tol)|...
        in_between(p2,bound_box_local_min,bound_box_local_max,tol) | ...
        in_between(p3,bound_box_local_min,bound_box_local_max,tol) | ...
        in_between(p4,bound_box_local_min,bound_box_local_max,tol) | ...
        in_between(p5,bound_box_local_min,bound_box_local_max,tol) | ...
        in_between(p6,bound_box_local_min,bound_box_local_max,tol) | ...
        in_between(p7,bound_box_local_min,bound_box_local_max,tol) | ...
        in_between(max_vec,bound_box_local_min,bound_box_local_max,tol);
    num_points_ref = sum(points_ref);
    num_local_bound_ref = sum(local_bound_ref);
    if (num_points_ref ~= 0) && (num_local_bound_ref ~= 0)
        %[num_points_ref num_local_bound_ref]
        local_ind_points = ind_points(points_ref);
        local_ind_el_box = int_el_box(local_bound_ref);
        num_points_ref = size((local_ind_points),1);
        num_local_bound_ref = size((local_ind_el_box),1);
        analyzed_points = ...
            repelem(points(points_ref,:),num_local_bound_ref,1);
        qq = [qq; analyzed_points];
        x_min = repmat(min_vec(local_ind_el_box,:),num_points_ref,1)...
            - tol;
        x_max = repmat(max_vec(local_ind_el_box,:),num_points_ref,1)...
            + tol;

        ii = repelem(local_ind_points, num_local_bound_ref,1);
        jj = repmat(local_ind_el_box, num_points_ref,1);
        %[analyzed_points x_min x_max]
        ref = all((analyzed_points - x_min) >= 0,2) & ...
            all((x_max - analyzed_points) >= 0,2);
        % if any(ref)
        %     local_ii = ii(ref);
        %     local_jj = jj(ref);
        % 
        %     %% end of bounding box test - check in polyhedra
        %     [~,~,lm_ref] = unique(local_jj);
        %     local_sign_matrix = sign_matrix(local_jj,:);
        %     testing_sign_matrix = local_sign_matrix(lm_ref,:);
        %     local_rock_el_faces = rock_el_faces(local_jj,:);
        %     %% for testing, need to be removed
        %     % c_faces = [c_faces ; c_faces];
        %     % c_faces(end,end) = 0;
        %     %%
        %     local_points_align = zeros(size(local_rock_el_faces));
        %     local_points = points(local_ii,:);
        % 
        %     for local_index = 1:size(local_rock_el_faces,2)
        %         %% extra attention here on  local_index and index
        %         cref = local_rock_el_faces(:,local_index) ~= 0;
        %         auxmat = local_rock_el_faces(cref,local_index);
        % 
        %         l_val = semi_plane_value(local_points, ...
        %             rock_normals(auxmat,:),  ...
        %             rock_points(auxmat,:));
        %         cref_negative = cref;
        %         cref(cref) = (l_val + tol >= 0);
        %         cref_negative(cref_negative) =  (l_val - tol <= 0);
        % 
        %         local_points_align(cref,local_index) = 1;
        %         local_points_align(cref_negative,local_index) = -1;
        %     end
        %     local_ref = all(testing_sign_matrix == local_points_align,2);
        %     ref(ref) = local_ref;
        %     is_in{bcount,1} = [ii(ref) jj(ref)];
        %     bcount = bcount + 1;
        % end
        if any(ref)
            is_in{bcount,1} = [ii(ref) jj(ref)];
            bcount = bcount + 1;
        end
    end
    %
1
end
is_in = cell2mat(is_in(1:bcount-1));
end





function [diag1,diag2] = divide_box(A, B, n)
    % Calculate the step size for each dimension
    step = (B - A) / n;
    % Initialize the boxes array
    boxes = zeros(n^3, 2, 3);
    % Loop through each dimension and calculate the coordinates of the boxes
    index = 1;
    for i = 0:n-1
        for j = 0:n-1
            for k = 0:n-1
                boxes(index, 1, :) = A + ...
                    [i*step(1) j*step(2) k*step(3)];
                boxes(index, 2, :) = A + ...
                    [(i+1)*step(1) (j+1)*step(2) (k+1)*step(3)];
                index = index + 1;
            end
        end
    end
    diag1 = squeeze(boxes(:, 1, :));
    diag2 = squeeze(boxes(:, 2, :));
end






