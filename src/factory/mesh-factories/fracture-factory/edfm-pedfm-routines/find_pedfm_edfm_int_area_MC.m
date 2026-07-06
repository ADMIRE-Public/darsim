%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel;
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [edfm, pedfm, frac_int_avg_d ] = ...
    find_pedfm_edfm_int_area_MC(fracmesh,rockmesh, ...
                                 frac_int, plane_geo, plane_int, settings)
%% this routine populates the fracmesh cells with a set of points and these
% points are used to approximate the intersection area between the
% rockmatrix mesh and fracture mesh. if pedfm flag is on these nodes are 
% projected on the faces of the pEDFM neighbors to also approximate this
% area.
%% rock_frac_int area is a cell is a cell array with the same number of
% entries as the the number of fractures in the fracture network. for each
% at each entrie you find the rock matrix indices and the fracm mesh
% indices and the corresponding intersection area
%% rock_frac_int area is a cell is a cell array with the same number of
% entries as the the number of fractures in the fracture network. for each
% at each entrie you find the rock matrix indices and the fracm mesh
% indices and the corresponding intersection area

%% is_in is a matrix that the first line consists on a volume of the rock 
% mesh while the second column is a reference to the test points generated
% inisde each fracture mesh element in frac_el. 
% is_in                                 frac_el
% [rock vol_ind test_point_ind]     [frac_el associated with test_point]
%% test_points is a matrix containing all the test_points used
% the montecarlos approximation points inside the corresponding rock
% volume. the EDFM algorithm eliminates the test points that are not inside
%% note:
% the routine will remove shorten the size of is_in but not the frac_el and
% test_points



edfm = struct();
pedfm = struct();
tol = 10^-13;
nfrac = max(size(fracmesh));
rock_ind_int_area = cell(nfrac,1);
frac_ind_int_area = cell(nfrac,1);
area = cell(nfrac,1);
rock_ind_int_area_pedfm = cell(nfrac,1);
rock_faces_ind_int_area_pedfm = cell(nfrac,1);
frac_ind_int_area_pedfm = cell(nfrac,1);
area_pedfm = cell(nfrac,1);
alpha_ind = cell(nfrac,1);
alpha_area = cell(nfrac,1);
[frac_avg_d_ind] = create_avg_d_frac_structure(frac_int);
frac_int_avg_d = cell(size(frac_int));

rock_normals =  [rockmesh.face.inner.normal;...
    rockmesh.face.bound.normal];
rock_face_cent =  [rockmesh.face.inner.centroid;...
    rockmesh.face.bound.centroid];
[min_vec, max_vec ] = ...
    generate_boundbox(rockmesh.element.vertices, ...
    rockmesh.vertex.coord, settings.useHalf);
for index = 1:nfrac
    [test_points, frac_el] =...
        generate_testpoint_parameters(fracmesh(index), ...
        plane_geo(index).normal,settings.num_area_points, ...
        settings.useHalf);
    % this test volume contain all the volumes cross by the fracture plane
    [test_volumes] = identify_plane_crossed_volumes(...
        rockmesh.element.vertices,...
        rockmesh.vertex.coord, plane_geo(index).eq);
    % preferably it would be idea to have the following function being
    % computed outside here. however the generation of points is the most
    % complex part. therefore adding it here can save a lot of time.
    if ~all(all(cellfun(@isempty,frac_int)))
      frac_int_avg_d = calculate_avg_frac_d(frac_int_avg_d, test_points,...
        frac_el, frac_avg_d_ind, frac_int, plane_int,index);
    end
    local_min_vec = min_vec(test_volumes,:);
    local_max_vec = max_vec(test_volumes,:);

    %%
    % from this point on test_volumes is comprised of only of the elements
    % that that are crossed by the fracture mesh
    %test_volumes(test_volumes) = in_frac;
    if sum(test_volumes) == 0
        warning('Fracture %d and rockmatrix do not intersect!\n', index);
        continue
    end
    abs_min = min(local_min_vec);
    abs_max = max(local_max_vec);
    bound_points_ref = in_between(test_points,abs_min,abs_max,tol);
    [boundp1,boundp2,boundp3,boundp4] = create_3d_boundbox_alg2(...
        plane_geo(index).coords, plane_geo(index).normal,...
        settings.num_frac_boundbox);
    rock_el_node_adj = rockmesh.element.vertices(test_volumes,:);
    rock_nodes = setdiff(unique(rock_el_node_adj),0);
    is_in = cell(size(boundp1,1),1);
    %% bounding box loop to manage memory issues
    for bound_index = 1:size(boundp1,1)
        lcoord = [boundp1(bound_index,:); boundp2(bound_index,:); ...
            boundp3(bound_index,:); boundp4(bound_index,:)];
        local_bound_points_ref = bound_points_ref;
        frac_nodes_in_boundbox =  ...
            in_planar_bounding_box(lcoord, plane_geo(index).normal, ...
            test_points(bound_points_ref,:));
        local_bound_points_ref(local_bound_points_ref) = ...
            frac_nodes_in_boundbox;
        local_test_volumes = test_volumes;
        
        %% check if the nodes of fracture boundbox are fully contained
        % inside the boundbox an element of rockmesh element
        rock_nodes_in_boundbox =  ...
            in_planar_bounding_box(lcoord, plane_geo(index).normal, ...
            rockmesh.vertex.coord(rock_nodes,:));
        local_test_volumes(local_test_volumes) = ...
            any(ismember(rock_el_node_adj,rock_nodes(...
            rock_nodes_in_boundbox)),2);

        %% rock volumes that are crossed by the fracture
        local_test_volumes = local_test_volumes;

        local_bound_point_ind = find_uint(local_bound_points_ref);
        local_test_volumes_ind = find_uint(local_test_volumes);
        ii = ...
            repmat(local_test_volumes_ind,size(local_bound_point_ind,1),1);
        jj = ...
            repelem(local_bound_point_ind,size(local_test_volumes_ind,1),1);
        %% apply bound box for better performance
        analyzed_points = test_points(jj,:);
        x_min = min_vec(ii,:) - tol;
        x_max = max_vec(ii,:) + tol;
        ref = all((analyzed_points - x_min) >= 0,2) & ...
            all((x_max - analyzed_points) >= 0,2);
        ii = ii(ref);
        jj = jj(ref);
        in_ref = is_in_vol3d(test_points(jj,:), ...
            rockmesh.element.faces(ii,:), ...
            rock_normals,rock_face_cent, ...
            rockmesh.element.centroid(ii,:));
        %% II -> RockCell Volume ID
        %% JJ -> Index to to the referred test_point
        %% to obtain frac_cell associated with point look for frac_el(JJ)
        is_in{bound_index} = [ii(in_ref) jj(in_ref)];
    end
    is_in = cell2mat(is_in);
    if isempty(is_in)
        warning('Fracture %d and rockmatrix do not intersect!\n', index);
        continue
    end
    %% core of the EDFM and pEDFM algorithms
    % the part above is used to generate the points
    [rock_frac_ind, int_area] = calculate_int_area_edfm(is_in, frac_el,...
        fracmesh(index).elements.area);
    rock_ind_int_area{index,1} = rock_frac_ind(:,1);
    frac_ind_int_area{index,1} = rock_frac_ind(:,2);
    area{index,1} = int_area;
    if settings.pedfm
        % RH: This is Artur's method for projection face selection,
        % replaced with Mousa's method for CPG. 
        % [rock_frac_ind_pdf,int_area_pedfm, rock_face_ind, rock_face_area]= ...
        %     find_pdefm_projections_area(rock_frac_ind, int_area, ...
        %     is_in,frac_el, test_points, rockmesh, plane_geo(index).eq);
        % RH: Mousa's method for projection face selection used in CPG. 
        [rock_frac_ind_pdf,int_area_pedfm, rock_face_ind, rock_face_area]= ...
            find_pdefm_projection_faces_centroid_intersect(rock_frac_ind, ...
            int_area, rockmesh, plane_geo(index).eq);
        rock_ind_int_area_pedfm{index,1} = rock_frac_ind_pdf(:,1); 
        frac_ind_int_area_pedfm{index,1} = rock_frac_ind_pdf(:,2);
        rock_faces_ind_int_area_pedfm{index} = rock_face_ind;
        area_pedfm{index,1} = int_area_pedfm;  
        alpha_ind{index,1} = rock_face_ind;
        alpha_area{index,1} = rock_face_area;
    end 
end
%% erro provavel aqui

alpha_area = cell2mat(alpha_area);
alpha_ind = cell2mat(rock_faces_ind_int_area_pedfm);
[alpha_ind,~,alpha_group] = unique(alpha_ind);
alpha_area = accumarray(alpha_group, alpha_area);


edfm.vol_ind = rock_ind_int_area;
edfm.frac_ind = frac_ind_int_area;
edfm.area = area;
if settings.pedfm
    pedfm.vol_ind = rock_ind_int_area_pedfm;
    pedfm.frac_ind = frac_ind_int_area_pedfm;
    pedfm.face_ind = rock_faces_ind_int_area_pedfm;
    pedfm.area = area_pedfm;
    pedfm.alpha_ind = alpha_ind;
    pedfm.alpha_area = alpha_area;
end
end


