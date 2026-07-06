classdef MeshCore < handle
    %MeshCore Class contains the decompressed MRST grid
    %   It is the core of the wrapper routines that dictate how information
    %   is accessed.
    properties      
        griddim
        operations
        faces
        nodes
        cells
        edges
        props        
        G
    end    
    properties (Access = private)
       is_cartgrid
       cart_dim
       store_G
    end
    methods
        function obj = MeshCore(streamG, varargin)
            %MeshCore Construct an instance of this class
            % Using MRST Routine to compute Geometry
            % compute geometry if not defined
            input_args = obj.processing_input(streamG,varargin);          
            G = input_args.G;
            obj.props = input_args.props;
            obj.store_G = input_args.storeG;            
            if ~any(strcmp(G.type,'computeGeometry'))
                G = computeGeometry(G);
            end
            if any(strcmp(G.type,"cartGrid"))
               obj.is_cartgrid = true;
               obj.cart_dim = G.cartDims;
            end
            obj.griddim = G.griddim;
            obj.operations = string(G.type);
            
            [obj.cells, obj.faces, obj.nodes, obj.edges] = ...
                                                grid_decompress(obj,G);  
            % adding face centroid
            % adding flag indicating boundary faces
            obj.faces.bound_faces = any(G.faces.neighbors == 0,2);
            % adding structured grid ijk references
            if obj.is_cartgrid
                obj.cells.indexMap = G.cells.indexMap;
            end
            if obj.store_G
                obj.G = G;
            end
        end

        function [ijk] = cartDim(obj,varargin)
            if obj.is_cartgrid
                if isempty(varargin)
                    [ijk{1:obj.griddim}] = ...
                        ind2sub(obj.cart_dim, obj.cells.indexMap);
                else
                    [ijk{1:obj.griddim}] = ...
                        ind2sub(obj.cart_dim, ...
                        obj.cells.indexMap(varargin{1}));
                end
                ijk = [ijk{:}];
            else
                disp("Mesh is not a cartesian grid.")
            end
        end
        function plot(obj,varargin)
            if obj.store_G
                plotGrid(obj.G, varargin{:});
                view(obj.griddim);
            else
                disp("MRST not stored.")
            end
        end
    end
    methods (Access = private)
        function [input_struct] =  processing_input(obj, G, input_stream)
            validationG = @(x) isstruct(x) | isstring(x);
            default_store_G = true;
            default_rock_prop = [];
            p = inputParser;
            addRequired(p,'G',validationG)
            addParameter(p, "props", default_rock_prop,@isstruct);
            addParameter(p, "extraProps", default_rock_prop,@isstruct);
            addParameter(p, "storeG", default_store_G, @islogical);
            parse(p, G, input_stream{:});
            input_struct = p.Results; 
            if isstring(input_struct.G)
                if isfile(input_struct.G)
                    [input_struct.G, input_struct.props] = ...
                        initialize_from_file(obj, input_struct.G);
                else
                    error("No MRST grid and Input Mehs Found.")
                end
            end
            if ~isempty(input_struct.extraProps)
                input_struct.props = combineStructs(input_struct.props, ...
                                            input_struct.extraProps);
            end
        
        end

        function [G, prop] = initialize_from_file(obj, file_path)
            [~,~,ext] = fileparts(file_path);
            switch ext
                case ".grdecl"
                    [G, prop] =  init_eclipse(obj,file_path);
                case ".msh"
                    error("To be added soon")
            end

        end
        function [G, prop] = init_eclipse(~,file_path)
            grdecl_file = readGRDECL(file_path);
            usys   = getUnitSystem('METRIC');
            grdecl_file = convertInputUnits(grdecl_file, usys);
            G = processGRDECL(grdecl_file, 'Verbose', true);
            prop = grdecl2Rock(grdecl_file, G.cells.indexMap);
        end

        function [cells, faces, nodes, edges] =  grid_decompress(obj, G)
            cells = struct();
            faces = struct();
            nodes = struct();
            edges = struct();
            [cells.element_faces_bool, cells.element_faces_flag] = ...
                                    decompress_cells_faces(obj,G);                        
            [faces.face_nodes_bool, faces.face_nodes] = ...
                decompress_faces_nodes(obj, G);
            cells.element_nodes_bool = logical(cells.element_faces_bool ...
                                        * faces.face_nodes_bool);            
            faces.centroid = calculate_face_centroid(G);
            [faces.face_neigh, faces.normals] = ...
                retrieve_face_neigh_normal(obj,G);
            nodes.coords = G.nodes.coords;
            if obj.griddim == 3
                [edges.edges_node, edges.centroid, faces.faces_edges, ...
                    faces.faces_edges_bool] =  generate_edges(obj , ...
                    faces.face_nodes, nodes.coords);
            end
            faces.areas = G.faces.areas;
            cells.volume = G.cells.volumes;
            cells.centroid = G.cells.centroids;
                   %% adding dimensions
            cells.num = size(cells.volume,1);
            faces.num = size(faces.normals,1);
            nodes.num = size(nodes.coords,1);
        end


        function [unique_edges, edges_centroid, faces_edges, ...
                faces_edges_bool] = generate_edges(~ , faces_nodes, coords)
                faces_ind = [1:size(faces_nodes,1)]';
                num_zeros = sum(faces_nodes == 0,2);
                all_face_types = unique(num_zeros)';
                edges_creation = [];
                faces_edges_ref = [];
                mlength = size(faces_nodes,2);
                for face_type = all_face_types
                    ref = num_zeros == face_type;
                    auxmat = faces_nodes(ref,1:(mlength-face_type));
                    auxmat = [auxmat auxmat(:,1)];
                    for col = 1:(size(auxmat,2)-1)
                        faces_edges_ref = ...
                            [faces_edges_ref; faces_ind(ref)];
                        edges_creation = [edges_creation; ...
                            auxmat(:,(col:col+1))];   
                    end                    
                end
                all_edges = sort(edges_creation,2);
                [unique_edges,~,edges_ind ] = unique(all_edges,"rows");
                faces_edges_bool = sparse(faces_edges_ref, edges_ind, ...
                    true, size(faces_ind,1), size(unique_edges,1));
                [ii, ~] = find(faces_edges_bool');
                faces_edges = sparse_to_zero_ind_mat( ...
                    reshape_jagged_array_in_matrix(faces_edges_bool, ii));
                edges_centroid = 0.5*(coords(unique_edges(:,1),:) + ...
                                            coords(unique_edges(:,2),:));
        end
        function [face_neigh, normals] = retrieve_face_neigh_normal(~,G)
                face_neigh = G.faces.neighbors;
                ref = face_neigh(:,1) == 0;
                auxvec = face_neigh(ref,2);
                face_neigh(ref,2) = face_neigh(ref,1);
                face_neigh(ref,1) = auxvec;
                normals = G.faces.normals;
                normals(ref,:) = -normals(ref,:);
        end

        function [face_nodes_bool, face_nodes] = ....
                decompress_faces_nodes(~, G)                 
            auxmat = ...
                [rldecode(1:G.faces.num, diff(G.faces.nodePos), 2)' ...
                G.faces.nodes ];    
           dim_cell_faces_nodes = max(auxmat(:,1:2));
            %up_lim = accumarray(auxmat(:,1),1);
            %low_lim = ones(size(up_lim));
            % saving the order of the nodes in the faces
            [~,ord_ref] = sortrows(auxmat);
            auxmat = auxmat(ord_ref,:);
            %face_nodes_order = create_multiple_range(low_lim,up_lim);
            face_nodes_bool = ...
                sparse(auxmat(:,1),auxmat(:,2), true,  ...
                dim_cell_faces_nodes(1),dim_cell_faces_nodes(2));
            face_nodes = ...
                sparse(auxmat(:,2),auxmat(:,1), G.faces.nodes,  ...
                dim_cell_faces_nodes(2),dim_cell_faces_nodes(1))'; 
            face_nodes = sparse_to_zero_ind_mat(face_nodes);
        end

        function [element_faces_bool, element_faces_flag] = ...
                decompress_cells_faces(~, G)       
            auxmat = [gridCellNo(G) G.cells.faces];
            dim_cell_faces = max(auxmat(:,1:2));
            element_faces_bool = ...
                sparse(auxmat(:,1),auxmat(:,2), true,  ...
                dim_cell_faces(1),dim_cell_faces(2));                
            element_faces_flag = [];
            if size(auxmat,2) > 1
                element_faces_flag = auxmat(:,end);
            end            
        end
    end
end

function [face_centroid] = calculate_face_centroid(G)
auxmat = ...
    [rldecode(1:G.faces.num, diff(G.faces.nodePos), 2)' ...
    G.faces.nodes ];
xx = accumarray(auxmat(:,1),G.nodes.coords(auxmat(:,2),1));
yy = accumarray(auxmat(:,1),G.nodes.coords(auxmat(:,2),2));
zz = [];
if G.griddim > 2
    zz = accumarray(auxmat(:,1),G.nodes.coords(auxmat(:,2),3));
end
n = accumarray(auxmat(:,1),1);
face_centroid = [xx yy zz] ./n ;
end

