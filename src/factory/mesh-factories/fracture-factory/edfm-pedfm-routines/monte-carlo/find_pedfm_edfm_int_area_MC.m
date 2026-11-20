%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel;
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [edfm, pedfm] = ...
    find_pedfm_edfm_int_area_MC(fracmesh, rockmesh, plane_geo, settings)
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

rock_normals =  [rockmesh.face.inner.normal;...
    rockmesh.face.bound.normal];
rock_face_cent =  [rockmesh.face.inner.centroid;...
    rockmesh.face.bound.centroid];
[min_vec, max_vec ] = ...
    generate_boundbox(rockmesh.element.vertices, rockmesh.vertex.coord);
if settings.useHalf
    min_vec = half(min_vec);
    max_vec = half(max_vec);
end
for index = 1:nfrac
    % this test volume contain all the volumes cross by the fracture plane
    [test_volumes] = identify_plane_crossed_volumes(...
        rockmesh.element.vertices,...
        rockmesh.vertex.coord, plane_geo(index).eq);
    if ~settings.pedfm
        %% EDFM areas only: test points are streamed in batches of
        % fracture cells and reduced to hit weights (bounded memory)
        [rock_frac_ind, int_area] = streamed_int_area_edfm( ...
            fracmesh(index), plane_geo(index), test_volumes, min_vec, ...
            max_vec, tol, rockmesh, rock_normals, rock_face_cent, settings);
        if isempty(rock_frac_ind)
            warning('Fracture %d and rockmatrix do not intersect!\n', index);
            continue
        end
        rock_ind_int_area{index,1} = rock_frac_ind(:,1);
        frac_ind_int_area{index,1} = rock_frac_ind(:,2);
        area{index,1} = int_area;
        continue
    end
    [test_points, frac_el] =...
        generate_testpoint_parameters(fracmesh(index), ...
        plane_geo(index).normal,settings.num_area_points, ...
        settings.useHalf);
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
        %% [RockCell Volume ID, test_point ID]; the frac cell of a point
        % is frac_el(test_point ID). Pairs are pruned per fracture cell
        % and tested in batches to bound memory.
        is_in{bound_index} = find_points_in_volumes(test_points, ...
            local_bound_point_ind, frac_el, local_test_volumes_ind, ...
            min_vec, max_vec, tol, rockmesh, rock_normals, rock_face_cent);
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
        [rock_frac_ind_pdf,int_area_pedfm, rock_face_ind, rock_face_area]= ...
            find_pdefm_projections_area(rock_frac_ind, int_area, ...
            is_in,frac_el, test_points, rockmesh, plane_geo(index).eq);
        rock_ind_int_area_pedfm{index,1} = rock_frac_ind_pdf(:,1); 
        frac_ind_int_area_pedfm{index,1} = rock_frac_ind_pdf(:,2);
        rock_faces_ind_int_area_pedfm{index} = rock_face_ind;
        area_pedfm{index,1} = int_area_pedfm;  
        alpha_ind{index,1} = rock_face_ind;
        alpha_area{index,1} = rock_face_area;
    end 
end
edfm.vol_ind = rock_ind_int_area;
edfm.frac_ind = frac_ind_int_area;
edfm.area = area;
if settings.pedfm
    % alpha faces only exist with pEDFM (EDFM-only runs failed here)
    alpha_area = cell2mat(alpha_area);
    alpha_ind = cell2mat(rock_faces_ind_int_area_pedfm);
    [alpha_ind,~,alpha_group] = unique(alpha_ind);
    alpha_area = accumarray(alpha_group, alpha_area, [numel(alpha_ind) 1]);
    pedfm.vol_ind = rock_ind_int_area_pedfm;
    pedfm.frac_ind = frac_ind_int_area_pedfm;
    pedfm.face_ind = rock_faces_ind_int_area_pedfm;
    pedfm.area = area_pedfm;
    pedfm.alpha_ind = alpha_ind;
    pedfm.alpha_area = alpha_area;
end
end


