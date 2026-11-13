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
        Settings
        FracGeometryBlueprint
        FracGeometry
        NumFrac
        EDFMType
        FracMesh
        RockMesh
        FracCore
    end
    properties (Access = private)
        GeoIntersection
    end

    methods 
        function obj = FractureFactory(FractureInit, ...
                                       RockMesh, ...
                                       Settings) 

            % Initialize General simulation settings
            obj.Settings = Settings;

            % Check and intialize fracture geometry blueprint
            obj.FracGeometryBlueprint = FractureInit.blueprint;               

            % Intialize fracture geometries based on the blueprint
            [obj.FracGeometry,obj.NumFrac] = obj.build_frac_geometry();

            % Initialize fracture simulation type: EDFM or pEDFM
            obj.EDFMType = obj.FracGeometryBlueprint.Type;

            % Intialize fracture meshes based on geometry and blueprint
            obj.FracMesh = obj.build_fracmesh(RockMesh);

            % Initialize rockmesh information
            obj.RockMesh = RockMesh;

            % Find fracture-fracture intersections based on geometry 
            % If there are more than one fracture
            LFracFracEDFM = cell(1);
            LFracFracPEDFM = cell(1);
            DFracFrac = cell(1);
            if obj.NumFrac > 1
                obj.GeoIntersection = ...
                    PlaneIntersection({obj.FracGeometry.normal}, ...
                                      [obj.FracGeometry.eq], ...
                                      {obj.FracGeometry.coords},...
                                      {obj.FracGeometry.vec_quad_edges});
                [LFracFracEDFM, LFracFracPEDFM] = ...
                    FracFracIntersection(obj.FracMesh, ...
                                         obj.GeoIntersection, ...
                                         obj.Settings);
                [DFracFrac] = ...
                    find_distance_frac_frac(LFracFracEDFM, ...
                                            LFracFracPEDFM, ...
                                            obj.FracMesh, ...
                                            obj.GeoIntersection, ...
                                            obj.Settings);
            end

            % Find fracture-rock intersections based on geometry
            [AFracRockEDFM, AFracRockPEDFM] = ...
                FracRockIntersection(obj.FracMesh, ...
                                     obj.RockMesh, ...
                                     obj.FracGeometry, ...
                                     obj.Settings);
            [DFracRock] = ...
                find_distance_frac_rock(AFracRockEDFM, ...
                                        AFracRockPEDFM, ...
                                        obj.RockMesh, ...
                                        obj.FracGeometry, ...
                                        obj.Settings);
            
            % Compute the physical properties of the (p)EDFM connections
            obj.FracCore = ...
                FindIntersectionPhysicalProperties(obj.RockMesh, ...
                                                   obj.FracMesh, ...
                                                   obj.FracGeometry, ...
                                                   LFracFracEDFM, ...
                                                   LFracFracPEDFM, ...
                                                   AFracRockEDFM, ...
                                                   AFracRockPEDFM, ...
                                                   DFracFrac, ...
                                                   DFracRock, ...
                                                   obj.Settings);               
        end

    end

    methods (Access = private)
        function [FracGeo,NumFrac] = build_frac_geometry(obj)
            fractures_geometry = ...
                FractureGeometryBuilder(obj.FracGeometryBlueprint);
            FracGeo = fractures_geometry.fracs(:);
            NumFrac = fractures_geometry.num_frac;            
        end
        function [Mesh] = build_fracmesh(obj,RockMesh)
            Mesh = FractureMeshBuilder(obj.FracGeometry, RockMesh);
            Mesh = Mesh.frac_mesh(:);
        end
    end
end
