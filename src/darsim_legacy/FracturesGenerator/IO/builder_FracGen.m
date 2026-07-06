% Builder Builds all objects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 13 July 2016
%Last modified: 4 March 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef builder_FracGen < handle
    properties
        MainInputData
        plotting
    end
    methods
        function Geometry = BuildGeometry(obj)
            switch obj.MainInputData.ReservoirProperties.Discretization
                case('CartesianGrid')
                    Geometry = geometry_CartesianGrid_FracGen();
                case('CornerPointGrid')
                    Geometry = geometry_CornerPointGrid_FracGen();
                case('TetrahedralGrid')
                    Geometry = geometry_TetrahedralGrid_FracGen();
            end
            
            Geometry.FractureModelType  = obj.MainInputData.Type;
            Geometry.DiscretizationType = obj.MainInputData.ReservoirProperties.Discretization;
            Geometry.NumOfDimensions    = obj.MainInputData.NumOfDimensions;
            Geometry.ReservoirGrid      = obj.BuildReservoir();
            %% builds fractures mesh and geometry
            Geometry.FracturesGrid      = obj.BuildFractures(Geometry.ReservoirGrid);
        end
        function Reservoir = BuildReservoir(obj)
            switch obj.MainInputData.ReservoirProperties.Discretization
                case('CartesianGrid')
                    Reservoir = reservoir_CartesianGrid_FracGen(obj.MainInputData.ReservoirProperties);
                case('CornerPointGrid')
                    Reservoir = reservoir_CornerPointGrid_FracGen(obj.MainInputData.ReservoirProperties);
                    Reservoir.CornerPointGridData = obj.MainInputData.ReservoirProperties.CornerPointGridData;
                case('TetrahedralGrid')
                    Reservoir = reservoir_TetrahedralGrid_FracGen(obj.MainInputData.ReservoirProperties);
                    Reservoir.TetrahedralGridData = obj.MainInputData.ReservoirProperties.TetrahedralGridData;
            end
            Reservoir.Discretize( obj.MainInputData.compareAccuracy , obj.MainInputData.Type );
        end
        function FracturesGrid = BuildFractures(obj, ReservoirGrid)
            FracturesGrid = fracturesGrid_FracGen( ReservoirGrid.Epsilon );
            FracturesGrid.BuildFracturesGrid( ReservoirGrid, obj.MainInputData.FracturesProperties );
        end
    end
end