%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel;
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [d] = calculate_fracture_average_distance_MC(edfm, pedfm, ...
                                    rockmesh, frac_geo, settings)
    n = settings.d_avg_refinement;
    tol = 10^-13;
    [vol_list, vol_list_by_frac] = ...
                    union_pedfm_edfm_new(edfm, pedfm, settings.pedfm);
    crossed_rock_elem = vol_list;
    [cross_rock_elem_ind] = [1:size(crossed_rock_elem,1)]';
    rock_normals =  [rockmesh.face.inner.normal;...
        rockmesh.face.bound.normal];
    rock_face_cent =  [rockmesh.face.inner.centroid;...
        rockmesh.face.bound.centroid];
    [min_vec, max_vec ] = ...
        generate_boundbox(rockmesh.element.vertices(crossed_rock_elem,:), ...
        rockmesh.vertex.coord, settings.useHalf);
    %% generating n x n x n  nodes inside each boundary box of the rock cells
    rock_ind = uint32(repelem(cross_rock_elem_ind,n*n*n,1));
    test_points = zeros(size(rock_ind,1),3);
    test_points(:,1) = min_vec(rock_ind,1) + ...
                        (max_vec(rock_ind,1) - ...
                        min_vec(rock_ind,1)).* ...
                        rand(size(rock_ind,1),1);
    test_points(:,2) = min_vec(rock_ind,2) + ...
                        (max_vec(rock_ind,2) - ...
                        min_vec(rock_ind,2)).* ...
                        rand(size(rock_ind,1),1);
    test_points(:,3) = min_vec(rock_ind,3) + ...
                        (max_vec(rock_ind,3) - ...
                        min_vec(rock_ind,3)).* ...
                        rand(size(rock_ind,1),1);
    rock_cell_id = crossed_rock_elem(rock_ind);
    
    %% filtering nodes outside the elements
    faces_matrix = rockmesh.element.faces(rock_cell_id,:);
    is_in_bool = false(size(test_points,1),size(faces_matrix,2));
    el_center = rockmesh.element.centroid(crossed_rock_elem(rock_ind),:);
    for index = 1:size(faces_matrix,2)
        ref = faces_matrix(:,index) ~= 0;
        is_in_bool(~ref,index) = 1;
        fadj = faces_matrix(ref,index);
        value = semi_plane_value(test_points(ref,:), rock_normals(fadj,:), ...
            rock_face_cent(fadj,:) );
        value_center = semi_plane_value(el_center(ref,:), ...
            rock_normals(fadj,:), rock_face_cent(fadj,:) );    
        is_in_bool(ref,index) = compare_signs(value, value_center, tol);
    end
    
    fref = all(is_in_bool,2);
    points_in_elements_ind = rock_cell_id(fref);
    test_points = test_points(fref,:);
    d = cell(size(frac_geo,1),2);
    for index = 1:size(frac_geo,1)
        local_ind = vol_list_by_frac{index};
        local_ind = unique(local_ind(:,1));
        [ref_points,~] = ismember(points_in_elements_ind,local_ind(:,1));
        local_points = test_points(ref_points,:);
        local_volume_ind = points_in_elements_ind(ref_points);   
        dist = distance_points_plane(local_points,frac_geo(index).normal, ...
                                           frac_geo(index).coords(1,:));
        [~,~,group_ref] = unique(local_volume_ind);
        num_points = accumarray(group_ref,1);
        dist_group = accumarray(group_ref,dist);
        d{index,1} = local_ind;
        d{index,2} = dist_group./num_points;   
    end
end

function [vol_list, vol_list_by_frac] = union_pedfm_edfm_new(edfm, ...
                                                        pedfm, flag)
    nfrac = max(size(edfm.vol_ind));
    vol_list_by_frac = cell(nfrac,1);
    for index = 1:nfrac
        lvol = unique(edfm.vol_ind{index,1});
        if flag
            lvol = union(lvol,unique(pedfm.vol_ind{index,1}));
        end
        vol_list_by_frac{index,1} = lvol;
    end
    vol_list = unique(cell2mat(vol_list_by_frac));
end

function [vol_list, vol_list_by_frac] = union_pedfm_edfm(edfm, pedfm, flag)
    nfrac = max(size(edfm.vol_ind));
    vol_list_by_frac = cell(nfrac,1);
    for index = 1:nfrac
        lvol = unique(edfm.vol_ind{index,1});
        if flag
            lvol = union(lvol,pedfm.vol_ind{index,1});
        end
        vol_list_by_frac{index,1} = lvol;
    end
    vol_list = unique(cell2mat(vol_list_by_frac));
end

function dist = distance_points_plane(points,normal,plane_point)
    A = normal(:,1);
    B = normal(:,2);
    C = normal(:,3);
    D = -dot(normal,plane_point,2);
    num = (A.*points(:,1) + B.*points(:,2) + C.*points(:,3)) + D;
    num = abs(num);
    den = (A^2 + B^2 + C^2)^.5;
    dist = num./den;
end



