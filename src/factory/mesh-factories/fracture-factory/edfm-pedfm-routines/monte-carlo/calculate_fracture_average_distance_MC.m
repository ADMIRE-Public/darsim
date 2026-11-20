%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel;
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [d] = calculate_fracture_average_distance_MC(edfm, pedfm, ...
                                    rockmesh, frac_geo, settings)
%% Average distance between every connected rock cell and the plane of
% its fracture.
% - Cells not crossed by the (infinite) fracture plane: the distance is a
%   linear function of position inside the cell, so its average is exactly
%   the distance of the cell centroid; no sampling needed.
% - Cells crossed by the plane: stratified Monte Carlo, one random point
%   in each of the n x n x n sub-boxes of the cell bounding box
%   (n = FRACTURE_AVERAGE_DISTANCE_REFINEMENT), keeping the points inside
%   the cell. Cells are processed in batches to bound memory.
    n = settings.d_avg_refinement;
    max_points = 2e6;
    [~, vol_list_by_frac] = ...
                    union_pedfm_edfm_new(edfm, pedfm, settings.pedfm);
    rock_normals =  [rockmesh.face.inner.normal;...
        rockmesh.face.bound.normal];
    rock_face_cent =  [rockmesh.face.inner.centroid;...
        rockmesh.face.bound.centroid];
    d = cell(numel(frac_geo),2);
    for index = 1:numel(frac_geo)
        local_ind = unique(vol_list_by_frac{index});
        local_ind = local_ind(:);
        normal = frac_geo(index).normal;
        plane_point = frac_geo(index).coords(1,:);
        avg_dist = distance_points_plane( ...
            rockmesh.element.centroid(local_ind,:), normal, plane_point);
        crossed = is_crossed_by_plane(rockmesh, local_ind, normal, ...
                                      plane_point);
        if any(crossed)
            [sampled, empty] = sampled_average_distance(rockmesh, ...
                local_ind(crossed), normal, plane_point, n, max_points, ...
                rock_normals, rock_face_cent);
            if any(empty)
                warning("DARSim:MonteCarlo:EmptyCell", ...
                    "%d rock cell(s) received no sample point for " + ...
                    "fracture %d; using the centroid distance. Increase " + ...
                    "FRACTURE_AVERAGE_DISTANCE_REFINEMENT.", ...
                    sum(empty), index);
                crossed_ind = find(crossed);
                sampled(empty) = avg_dist(crossed_ind(empty));
            end
            avg_dist(crossed) = sampled;
        end
        d{index,1} = local_ind;
        d{index,2} = avg_dist;
    end
end

function crossed = is_crossed_by_plane(rockmesh, cells, normal, plane_point)
% true when the cell has vertices strictly on both sides of the plane
    vert = rockmesh.element.vertices(cells,:);
    valid = vert ~= 0;
    vert(~valid) = 1;
    side = reshape((rockmesh.vertex.coord(vert,:) - plane_point) * ...
                   (normal(:) / norm(normal)), size(vert));
    side(~valid) = NaN;
    tol = 1e-10 * max(abs(side), [], 2, 'omitnan');
    crossed = max(side, [], 2, 'omitnan') > tol & ...
              min(side, [], 2, 'omitnan') < -tol;
end

function [avg_dist, empty] = sampled_average_distance(rockmesh, cells, ...
            normal, plane_point, n, max_points, rock_normals, rock_face_cent)
    tol = 10^-13;
    ncell = numel(cells);
    [min_vec, max_vec] = generate_boundbox( ...
        rockmesh.element.vertices(cells,:), rockmesh.vertex.coord);
    [gi, gj, gk] = ndgrid(0:n-1);
    strata = [gi(:) gj(:) gk(:)];
    dist_sum = zeros(ncell,1);
    num_points = zeros(ncell,1);
    chunk = max(1, floor(max_points / n^3));
    for first = 1:chunk:ncell
        local = (first:min(ncell, first + chunk - 1))';
        cell_ind = repelem(local, n^3, 1);
        test_points = min_vec(cell_ind,:) + ...
            (repmat(strata, numel(local), 1) + rand(numel(cell_ind),3)) ...
            / n .* (max_vec(cell_ind,:) - min_vec(cell_ind,:));
        inside = points_in_cells(test_points, cells(cell_ind), ...
                                 rockmesh, rock_normals, rock_face_cent, tol);
        dist = distance_points_plane(test_points(inside,:), normal, ...
                                     plane_point);
        dist_sum = dist_sum + ...
            accumarray(cell_ind(inside), dist, [ncell 1]);
        num_points = num_points + ...
            accumarray(cell_ind(inside), 1, [ncell 1]);
    end
    empty = num_points == 0;
    avg_dist = dist_sum ./ num_points;
end

function inside = points_in_cells(test_points, cell_id, rockmesh, ...
                                  rock_normals, rock_face_cent, tol)
    faces_matrix = rockmesh.element.faces(cell_id,:);
    el_center = rockmesh.element.centroid(cell_id,:);
    inside = true(size(test_points,1),1);
    for index = 1:size(faces_matrix,2)
        ref = faces_matrix(:,index) ~= 0;
        fadj = faces_matrix(ref,index);
        value = semi_plane_value(test_points(ref,:), ...
            rock_normals(fadj,:), rock_face_cent(fadj,:));
        value_center = semi_plane_value(el_center(ref,:), ...
            rock_normals(fadj,:), rock_face_cent(fadj,:));
        inside(ref) = inside(ref) & compare_signs(value, value_center, tol);
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



