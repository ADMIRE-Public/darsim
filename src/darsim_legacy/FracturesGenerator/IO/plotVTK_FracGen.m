% VTK Plotter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 11 July 2016
%Last modified: 17 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef plotVTK_FracGen < handle
    properties
        FileName
        VTKindex 
        isBinary
        ReferenceDepth
    end
    methods
        function obj = plotVTK_FracGen(Directory, Problem)
            if ~exist(strcat(Directory,'/VTK/'), 'dir')
                mkdir(Directory,'/VTK');
            else
                delete(strcat(Directory,'/VTK/*.vtk'));
            end
            obj.FileName = strcat(Directory, '/VTK/', Problem);
            obj.VTKindex = 1;
        end
        function PlotSolution(obj, Geometry)
            obj.PlotReservoirSolution(Geometry.ReservoirGrid);
            for f = 1 : Geometry.FracturesGrid.N_Frac
                obj.PlotFracturesSolution(Geometry.FracturesGrid.Fracture(f), f);
            end
            obj.VTKindex = obj.VTKindex + 1;
        end
        function PrintScalar2VTK(obj, fileID, scalar, name)
            %Print a scalar in VTK format
            fprintf(fileID, '\n');
            fprintf(fileID, strcat('SCALARS  ', name,' double 1\n'));
            fprintf(fileID, 'LOOKUP_TABLE default\n');
            if obj.isBinary
                fwrite(fileID, scalar', 'double', 'b');
            else
            	fprintf(fileID,'%1.5e ', scalar);
            end
        end
    end
end