%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Ryan Haagenson; Artur Castiel
%TU Delft
%
% FractureFactory class for constructing the fracture geometry, fracture 
% mesh, and fracture-rock and fracture-fracture intersections required for 
% (p)EDFM simulations.
%
% Input:
%   fracture_init:      struct
%                       Contains basic fracture information
%   rockmesh:           struct
%                       Contains relevant rock mesh information
%   settings:           struct
%                       Contains relevant simulation settings
%
% Output:
%   FractureFactory:    object
%                       Contains fracture geometry, mesh, intersection,
%                       cross connection, and physical property information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef FractureFactory
    properties
        geometry_blueprint
        geometry
        settings
        frac_mesh_intersection
        fracmesh
        file
        num_frac
        type 
        avg_d
        pedfm
        edfm
        frac_core
        K
        rockmesh
    end
    properties (Access = private)
        geo_intersection
        unstructured_fracmesh
    end

    methods 
        function obj = FractureFactory(fracture_init, rockmesh,...
                                                    settings, varargin) 

            % Initialize fracmesh type: ustructured (1) or structured (0)
            obj.unstructured_fracmesh =  fracture_init.UnsFracmesh;

            % Initialize General simulation settings
            obj.settings = settings;

            % Initialize rockmesh information
            obj.rockmesh = rockmesh;

            % Check and intialize fracture geometry blueprint if not provided
            if isstruct(fracture_init) || isobject(fracture_init)
                obj.geometry_blueprint = fracture_init.blueprint;
            elseif isfile(fracture_init)
                fractures_description = FracReader(fracture_init);
                obj.file = fracture_init;
                obj.geometry_blueprint = fractures_description.blueprint;                
            else
                error("ERROR: Fracture network failed to intialize. " + ...
                    " No fracture geometry blueprint found.")
            end

            % Initialize fracture simulation type: EDFM or pEDFM
            obj.type = obj.geometry_blueprint.Type;

            % Intialize fracture geometries based on the blueprint
            [obj.geometry,obj.num_frac] = obj.build_geometry();

            % Intialize fracture meshes based on geometry and blueprint
            [obj.fracmesh] = obj.build_fracmesh(rockmesh);

            % Find fracture intersections based on geometry 
            % If there are more than one fracture
            frac_mesh_dist_int = cell(1);
            if obj.num_frac > 1
                obj.geo_intersection = obj.find_geo_intersection();
                [frac_mesh_dist_int] = ...
                    obj.create_frac_mesh_intersection(obj.geo_intersection);
            end

            % Find intersections between fractures and rock matrix
            [edfm_geo, pedfm_geo, avg_d_rock_geo, avd_d_frac_geo] = ...
                        obj.create_frac_rock_intersection(rockmesh, ...
                                    frac_mesh_dist_int, settings);
            
            % Compute the physical properties of the (p)EDFM connections
            obj.frac_core = compute_physical_prop_pedfm(rockmesh,...
                        obj.fracmesh, obj.geometry, frac_mesh_dist_int, ...
                        edfm_geo, pedfm_geo, avg_d_rock_geo, ...
                        avd_d_frac_geo );               
        end
    end
    methods (Access = private)
        function [edfm_geo, pedfm_geo, d_rock,d_frac] = ...
                    create_frac_rock_intersection(obj,rockmesh,...
                                                frac_mesh_dist_int,settings)
            if isfield(settings,"edfm_method")
                if settings.edfm_method == "monte_carlo"
                    [edfm_geo, pedfm_geo, d_frac] = ...
                        find_pedfm_edfm_int_area_MC(...
                        obj.fracmesh, ...
                        rockmesh, ...
                        frac_mesh_dist_int, ...
                        obj.geometry, ...
                        obj.geo_intersection, ...
                        settings);
                    [d_rock] = ...
                        calculate_fracture_average_distance_MC(...
                        edfm_geo, ...
                        pedfm_geo, ...
                        rockmesh, ...
                        obj.geometry, ...
                        settings);

                elseif settings.edfm_method == "geometric"
                    [edfm_geo, pedfm_geo] = ...
                        find_pedfm_edfm_int_area_geometric(...
                        obj.fracmesh, ...
                        rockmesh, ...
                        obj.geometry, ...
                        settings);
                    [d_frac] = find_distance_frac_to_frac_int(...
                        obj.fracmesh, ...
                        frac_mesh_dist_int, ...
                        obj.geometry, ...
                        obj.geo_intersection,settings);
                    [d_rock] = ...
                        calculate_fracture_average_distance_grid(...
                        edfm_geo, ...
                        pedfm_geo, ...
                        rockmesh, ...
                        obj.geometry, ...
                        settings);
                else
                    error("ERROR: EDFM method not recognized. " + ...
                        "Available EDFM_METHOD options: " + ...
                        "monte_carlo, geometric ")
                end
            else
                error("ERROR: No EDFM method found. " + ...
                    "Include EDFM method desired in input file.")
            end

        end

        function [geo,num_frac] = build_geometry(obj)
            fractures_geometry = ...
                FractureGeometryBuilder(obj.geometry_blueprint);
            geo = fractures_geometry.fracs(:);
            num_frac = fractures_geometry.num_frac;            
        end
        function [mesh] = build_fracmesh(obj,rockmesh)
            mesh = FractureMeshBuilder(obj.geometry, ...
                obj.unstructured_fracmesh,rockmesh);
            mesh = mesh.frac_mesh(:);
        end    
        function geo_intersection = find_geo_intersection(obj)            
            geo_intersection = PlaneIntersection({obj.geometry.normal},...
                    [obj.geometry.eq], {obj.geometry.coords},...
                    {obj.geometry.vec_quad_edges});
        end
        function [frac_mesh_int] = ...
            create_frac_mesh_intersection(obj, geo_intersection)
            frac_int_obj = FracMeshIntersection(geo_intersection, ...
                obj.geometry, obj.fracmesh);
            frac_mesh_int = frac_int_obj.meshFracInterval;
        end

    end
end
