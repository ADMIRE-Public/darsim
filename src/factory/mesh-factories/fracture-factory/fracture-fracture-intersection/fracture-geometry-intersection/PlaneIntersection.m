%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel;
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef PlaneIntersection
    % line_intersection is written as
    % x = p + v*t
    % p and v are stored in line_parameters [p;v]
    properties
        intersectWith        
    end
    properties (Access = private)
        not_parallel
        line_parameters
        edges_vec
        edges_intersection
        normals
        coords
        eqs
        num_frac
    end

    methods
        function obj = PlaneIntersection(normals, eqs, nodes_coord,...
                edges_vec)
            obj.normals = cell2mat(normals(:));
            obj.num_frac = size(obj.normals,1);
            obj.coords = (nodes_coord(:));
            obj.edges_vec = (edges_vec(:));
            obj.eqs = eqs;
            obj.not_parallel = obj.check_which_planes_are_not_parallel();
            obj.line_parameters = obj.find_plane_intersection();
            edges_intersection = obj.find_plane_edges_line_interesction();
            obj.edges_intersection = ...
               obj.verify_if_intersection_are_in_plane(edges_intersection);
            obj.intersectWith = ~cellfun(@isempty,obj.edges_intersection);
            obj.intersectWith  = obj.intersectWith | obj.intersectWith';
            
        end
        function [coords] = get_boundary_line_intersection(obj,ii,jj)
            coords = [];
            if ii > jj
                [jj, ii] = deal(ii,jj);
            end
            if obj.not_parallel(ii,jj)
                coords = obj.edges_intersection{ii,jj};
                if isempty(coords)
                    coords = obj.edges_intersection{jj,ii};
                end
            end
        end
        function [p,v] = get_line(obj,ii,jj)
            if ii > jj
                [jj, ii] = deal(ii,jj);
            end
            if obj.not_parallel(ii,jj)
                line_par = obj.line_parameters{ii,jj};
                p = line_par(1,:);
                v = line_par(2,:);
            else
                p = []; v = [];
            end
        end
    end
    methods (Access = private)
        function [edges_int_update] = ...
                        verify_if_intersection_are_in_plane(obj,edges_int)
            edges_int_update = edges_int;
            [ii,jj] = find(~cellfun(@isempty,edges_int));
            fract_int = [ii,jj];
            for index = 1:size(fract_int,1)           
                fracture_ii = fract_int(index,1);
                fracture_jj = fract_int(index,2);
                frac_ii_coords = obj.coords{fracture_ii};
                frac_jj_coords = obj.coords{fracture_jj};
                intersection_ii_jj_line = ...
                 edges_int{fracture_ii,fracture_jj};
                points_ref = ~cellfun(@isempty,intersection_ii_jj_line);
                ref_ii = check_point_in_polygon(frac_ii_coords, ...
                                             intersection_ii_jj_line);                
                ref_jj = check_point_in_polygon(frac_jj_coords, ...
                                             intersection_ii_jj_line);

                if sum(ref_jj) == 1
                    % The other end point of intersection might exist in
                    % other fracture
                    intersection_jj_ii_line = ...
                        edges_int{fracture_jj,fracture_ii};
                    points_ref_other = ~cellfun(@isempty,...
                        intersection_jj_ii_line);
                    ref_other = check_point_in_polygon(frac_ii_coords, ...
                                             intersection_jj_ii_line);
                    % To Do: Replace the one that is not in jj with the 
                    % one that is.
                    points_ref_other(points_ref_other) = ref_other;
                    replacement_point = ...
                        intersection_jj_ii_line(points_ref_other);
                    ref = ref_ii & ~ref_jj;
                    points_ref(points_ref) = ref;
                    intersection_ii_jj_line(points_ref) = ...
                        replacement_point;
                else
                    ref = ref_ii & ref_jj;
                    points_ref(points_ref) = ref;
                    intersection_ii_jj_line(~points_ref) = ...
                        cell(sum(~points_ref),1);                
                end
                if all(cellfun(@isempty,intersection_ii_jj_line))
                    edges_int_update{fracture_ii,fracture_jj} = [];
                else
                    edges_int_update{fracture_ii,fracture_jj} = ...
                                                   intersection_ii_jj_line;
                end

            end
        end
        function [plane_edge_line_mat] = ...
                find_plane_edges_line_interesction(obj)
            plane_edge_line_mat = cell(obj.num_frac,obj.num_frac);
            for index = 1:obj.num_frac
                P1 = obj.coords{index};
                P2 = circshift(P1,-1,1);
                V = obj.edges_vec{index};
                plane_int = find(obj.not_parallel(index,:));
                for neigh = plane_int
                    [pl,vl] = obj.get_line(index, neigh);
                    auxcell = find_line_intersection_one(P1,P2, V,pl,vl);
                    if ~all(cellfun(@isempty,auxcell))
                        plane_edge_line_mat(index, neigh) = {auxcell};
                    end
                end
            end
        end
        function int_mat = check_which_planes_are_not_parallel(obj)
            int_mat = false(obj.num_frac,obj.num_frac);
            frac_int = nchoosek(1:obj.num_frac,2);
            for index = 1:size(frac_int,1)
                ii = frac_int(index,1);
                jj = frac_int(index,2);
                isparallel = ...
                abs(dot(obj.normals(ii,:), obj.normals(jj,:)) - 1)> 10^-10;
                int_mat(ii,jj) = isparallel;
                int_mat(jj,ii) = isparallel;
            end
        end
        function line_intersect = find_plane_intersection(obj)
            [ii,jj] = find(obj.not_parallel);
            line_intersect = cell(obj.num_frac);
            fract_int = [ii,jj];
            for index = 1:size(fract_int,1)
                ii = fract_int(index,1);
                jj = fract_int(index,2);
                [p,v] = find_intersection_of_plane_eq(...
                    obj.eqs(ii),obj.eqs(jj));
                line_intersect{ii,jj} = [p;v];
            end
        end

    end
end

function intersection_nodes = find_line_intersection_one(P1,P2, V,pl,vl)
intersection_nodes = cell(size(P1,1),1);
    pl = pl';
    vl = vl';
    warning('off','MATLAB:rankDeficientMatrix')
    for index = 1:size(P1,1)
        p1 = P1(index,:)';
        p2 = P2(index,:)';
        v = V(index,:)';
        %Ax=b
        b = pl - p1; 
        A = [v, -vl];
        t = mldivide(A,b);
        pq = (p1 + t(1)*v)';
        pz = (pl + t(2)*vl)';
        if abs(pq - pz) < 10^-11
            d = vecnorm((p2' - p1'),1,2);
            d1 = vecnorm((pq - p1'),1,2);
            d2 = vecnorm((pq - p2'),1,2);
            same_points = all(abs(pq - pz) < 10^-12,2);
            inside_line_segment = abs((d1+d2) - d) < 10^-12;
            if inside_line_segment && same_points
                intersection_nodes{index} = pq;
            end
        end
    end
    warning('on','MATLAB:rankDeficientMatrix')
end

function [p,v] = find_intersection_of_plane_eq(eq1,eq2)
    LHS = [eq1.A eq1.B eq1.C; eq2.A eq2.B eq2.C];
    RHS = -[eq1.D;eq2.D];
    p = linsolve(LHS,RHS)';
    v = cross(LHS(1,:),LHS(2,:));
end
