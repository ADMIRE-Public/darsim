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
        fracmesh
        file
        num_frac
        type 
        frac_core
        rockmesh
    end
    properties (Access = private)
        geo_intersection
    end

    methods 
        function obj = FractureFactory(fracture_init, rockmesh,...
                                                    settings, varargin) 

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

            % Find fracture-fracture intersections based on geometry 
            % If there are more than one fracture
            frac_frac_int_L_edfm = cell(1);
            if obj.num_frac > 1
                obj.geo_intersection = obj.find_geo_intersection();
                [frac_frac_int_L_edfm] = ...
                    obj.create_frac_frac_intersection(obj.geo_intersection);
                [frac_frac_int_d] = ...
                    find_distance_frac_frac_int(obj.fracmesh, ...
                                                   frac_frac_int_L_edfm, ...
                                                   obj.geometry, ...
                                                   obj.geo_intersection,...
                                                   settings);
            end

            % Find fracture-rock intersections based on geometry
            [frac_rock_int_A_edfm, frac_rock_int_A_pedfm] = ...
                obj.create_frac_rock_intersection(rockmesh, settings);
            [frac_rock_int_d] = ...
                find_distance_frac_rock_int(frac_rock_int_A_edfm, ...
                                               frac_rock_int_A_pedfm, ...
                                               rockmesh, ...
                                               obj.geometry, ...
                                               settings);
            
            % Compute the physical properties of the (p)EDFM connections
            obj.frac_core = compute_physical_prop_pedfm(rockmesh,...
                        obj.fracmesh, obj.geometry, frac_frac_int_L_edfm, ...
                        frac_rock_int_A_edfm, frac_rock_int_A_pedfm, ...
                        frac_frac_int_d, frac_rock_int_d );               
        end
    end

    methods (Access = private)
        function [frac_frac_int_L_edfm] = ...
            create_frac_frac_intersection(obj, geo_intersection)
            frac_int_obj = FracMeshIntersection(geo_intersection, ...
                obj.geometry, obj.fracmesh);
            frac_frac_int_L_edfm = frac_int_obj.meshFracInterval;
        end
        function [frac_rock_int_A_edfm, frac_rock_int_A_pedfm] = ...
                    create_frac_rock_intersection(obj,rockmesh,settings)
            [frac_rock_int_A_edfm, frac_rock_int_A_pedfm] = ...
                find_pedfm_edfm_int_area(...
                                         obj.fracmesh, ...
                                         rockmesh, ...
                                         obj.geometry, ...
                                         settings);
        end
        function [geo,num_frac] = build_geometry(obj)
            fractures_geometry = ...
                FractureGeometryBuilder(obj.geometry_blueprint);
            geo = fractures_geometry.fracs(:);
            num_frac = fractures_geometry.num_frac;            
        end
        function [mesh] = build_fracmesh(obj,rockmesh)
            mesh = FractureMeshBuilder(obj.geometry, rockmesh);
            mesh = mesh.frac_mesh(:);
        end    
        function geo_intersection = find_geo_intersection(obj)            
            geo_intersection = PlaneIntersection({obj.geometry.normal},...
                    [obj.geometry.eq], {obj.geometry.coords},...
                    {obj.geometry.vec_quad_edges});
        end

    end
end
