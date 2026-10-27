%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel;
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef FractureMeshBuilder
    %FractureBuilder created for quadrilateral fractures discretized with structured grids only
    properties
        num_frac
        frac_mesh
        mesh_type
        dV
    end

    methods
        function obj = FractureMeshBuilder(frac_geometry, rockmesh)
            obj.dV = mean(rockmesh.element.volume);
            obj.num_frac = max(size(frac_geometry));
            n_xy= obj.define_num_mesh_elemements(frac_geometry);
            obj.frac_mesh = obj.build_fractures(n_xy,frac_geometry);
        end
    end

    methods (Access = private)
        function frac = build_fractures(obj, n_xy, frac_geo)
            frac(obj.num_frac) = struct("vertex","","faces","", ...
                                                            "elements","");
            for index = 1:obj.num_frac
             frac(index) = ...
                create_frac_mesh(n_xy(index,:), frac_geo(index).coords, ...
               frac_geo(index).normal,frac_geo(index).vec_quad_edges);            
            end
        end
        function n_xy = define_num_mesh_elemements(obj, fracs)
            n_xy= zeros(obj.num_frac,2);
            avgL = obj.dV ^(1/3);
            for index = 1:obj.num_frac
                if ~isnan(fracs(index).GridResRatio)
                    n_xy(index,:) = ...
                        round([fracs(index).edges(1:2)'] * avgL * ...
                        fracs(index).GridResRatio + ...
                        eps([fracs(index).edges(1:2)'] * avgL * ...
                        fracs(index).GridResRatio));
                else
                    n_xy(index,:) = [fracs(index).GridNumAlongL ...
                        fracs(index).GridNumAlongW];
                end
            end
        end
    end
end

function [mesh] = create_frac_mesh(n_xy,geo_coords, normal, geo_vec)
    mesh = struct();
    [coords, elem] = ...
        create_frac_coords_and_elem(n_xy,geo_coords, geo_vec, normal);
    faces = create_faces(elem,coords);
    [element_face_adj,element_face_adj_bool] = ...
                                  create_elem_face_adj(faces.inner.upstream,...
                    faces.inner.downstream,faces.bound.upstream, size(elem,2));
    mesh.vertex.coords = coords;
    mesh.faces = faces;
    mesh.elements.vertices = elem;
    mesh.elements.faces_bol = element_face_adj_bool;
    mesh.elements.faces = element_face_adj;
    mesh.elements.area = calc_area_convex_polygon2d(mesh.elements.vertices, ...
                                                           mesh.vertex.coords);
    
    mesh.elements.centroid = centroid_calc(coords,elem);

end

function [centroid] = centroid_calc(coords,adj)
    centroid = zeros(size(adj,1),3);
    num_verticies = sum(adj~=0,2);
    for index = 1:size(adj,2)
        ref = adj(:,index) ~= 0;
        centroid(ref,:) = centroid(ref,:) + coords(adj(ref,index),:);
    end
    centroid = centroid./num_verticies;
end

function [coords, elem] = ...
   create_frac_coords_and_elem(n_xy,geo_coords, geo_vec, normal)
        [coords,elem] = create_structured_2d_mesh(n_xy, geo_coords, ...
            geo_vec);
end



function [coords, elem] = create_structured_2d_mesh(n_xy, coords, vec)
    multiplierX = [0:n_xy(1)]';
    multiplierY = [0:n_xy(2)]';
    direction_vectors = vec ./ [n_xy n_xy]';
    length_coords = repelem(coords(1,:),n_xy(1)+1,1) + ...
        (repelem(direction_vectors(1,:),n_xy(1)+1,1) .* multiplierX);
    mesh_coords = repmat(length_coords, 1, n_xy(2)+1) + ...
        repmat(reshape((repelem(direction_vectors(2,:), ...
        n_xy(2)+1,1).*multiplierY)',1,[]),n_xy(1)+1,1);
    %lm = 1:3:(4*n_xy(2)-1);
    lm = 1:3:size(mesh_coords,2);
    X = mesh_coords(:,lm);
    Y = mesh_coords(:,lm+1);
    Z = mesh_coords(:,lm+2);
    elem_ind = 1:prod(n_xy);
    %% creating coords
    coords = [reshape(X,[],1) reshape(Y,[],1) reshape(Z,[],1)];
    
    %% convert and find the indices based on the prefered order JJ,II
    [ii,jj] = ind2sub([n_xy(1),n_xy(2)] ,elem_ind');
    
    %%  find the indices of the coords in terms of JJ,II
    sqr_node_1 = [ii, jj];
    sqr_node_2 = [ii+1, jj];
    sqr_node_3 = [ii+1, jj+1];
    sqr_node_4 = [ii, jj+1];
    
    %% convert it to a 1:ncoord
    node_ref_1 = ...
        sub2ind([n_xy(1)+1,n_xy(2)+1] ,sqr_node_1(:,1),sqr_node_1(:,2));
    node_ref_2 = ...
        sub2ind([n_xy(1)+1,n_xy(2)+1] ,sqr_node_2(:,1),sqr_node_2(:,2));
    node_ref_3 = ...
        sub2ind([n_xy(1)+1,n_xy(2)+1] ,sqr_node_3(:,1),sqr_node_3(:,2));
    node_ref_4 = ...
        sub2ind([n_xy(1)+1,n_xy(2)+1] ,sqr_node_4(:,1),sqr_node_4(:,2));
    %% create elem
    elem = [node_ref_1 node_ref_2 node_ref_3 node_ref_4];
end

function [coords, elem] = ...
    create_unstructured_2d_mesh_density(n_xy, geo_coords, geo_vec,normal)
    mesh_avg_size = vecnorm(geo_vec(1:2,:),2,2) ./ n_xy';
    mesh_avg_size = mesh_avg_size/3;
    center = mean(geo_coords);
    centerm = repmat(center,size(geo_coords,1),1);
    project = geo_coords - centerm;
    normal = normal';
    rotM = rotate_to_z_axis(normal);
    ncoord = (rotM*project')' + centerm;
    bshape = polyshape(ncoord(:,1), ncoord(:,2));
    tr = triangulation(bshape);
    tnodes = tr.Points';
    telements = tr.ConnectivityList';
    model = createpde;
    geometryFromMesh(model,tnodes,telements);
    mesh = generateMesh(model,"GeometricOrder","linear","Hgrad",1.75,...
       "Hmax",max(mesh_avg_size),"Hmin",min(mesh_avg_size));
    elem = mesh.Elements';
    ancoord = [mesh.Nodes' center(3)*ones(size(mesh.Nodes,2),1)];
    centerm = repmat(center,size(ancoord,1),1);
    projectb = ancoord - centerm;
    coords = ((rotM^-1)*projectb')' + centerm;
end

function faces = create_faces(elem, coords)
    faces = struct();
    elem_ind = [1:size(elem,1)]';
    faces_ind = [[1:size(elem,2)]', circshift([1:size(elem,2)]',-1)];
    elem_face_nodes = [];
    for index = 1:size(faces_ind,1)
        [ii] = faces_ind(index,1);
        [jj] = faces_ind(index,2);
        elem_face_nodes = [elem_face_nodes; ...
                        elem_ind sort([elem(:,ii) elem(:,jj)],2)];
    end
    [unique_faces_adj,~,faces_class] = unique(elem_face_nodes(:,2:3),'rows');
    face_elem = sparse(elem_face_nodes(:,1), faces_class,true)';
    ref = sum(face_elem,2) == 2;
    inner = struct();
    bound = struct();
    inner.vertices =  unique_faces_adj(ref,:);
    P1 = coords(inner.vertices(:,1),:);
    P2 = coords(inner.vertices(:,2),:);
    inner.centroid = 0.5*(P1+P2);
    inner.length = vecnorm(P2-P1,2,2);
    
    [neigh,~]= find(face_elem(ref,:)');
    auxmat = reshape(neigh,2,[])';
    inner.downstream =  auxmat(:,1);
    inner.upstream =  auxmat(:,2);
    bound.vertices =  unique_faces_adj(~ref,:);
    P1 = coords(bound.vertices(:,1),:);
    P2 = coords(bound.vertices(:,2),:);
    bound.centroid = 0.5*(P1+P2);
    bound.length = vecnorm(P2-P1,2,2);
    
    [neigh,~]= find(face_elem(~ref,:)');
    bound.upstream = neigh;
    faces.inner = inner;
    faces.bound = bound;
end
