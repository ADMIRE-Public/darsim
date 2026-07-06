% Fracture Generator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Mousa HosseiniMehr
%TU Delft
%Created: 2017-03-10
%Last modified: 2017-03-10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef pEDFM_Generator < handle
    properties
        Reader
        Builder
        Geometry
        Writer 
    end
    methods
        function obj = pEDFM_Generator(Directory, File)
            obj.Reader  = reader_FracGen(Directory, File);
            obj.Builder = builder_FracGen();
            obj.Writer = writer_FracGen(Directory,'Fracture_Output.txt');
        end
        function ReadInputFile(obj)
            obj.Reader.ReadInputFile(obj.Builder);
        end
        function BuildObjects(obj)
            obj.Geometry = obj.Builder.BuildGeometry();
        end
        function PrintInfo(obj)
            obj.Geometry.ReservoirGrid.PrintInfo();
            obj.Geometry.FracturesGrid.PrintInfo();
        end
        function Run(obj)
            obj.Geometry.PrintInfo();
            obj.Geometry.Compute_EDFM_Connectivities();
            if strcmp(obj.Geometry.FractureModelType,'pEDFM')
                obj.Geometry.Add_pEDFM_Connectivities();
            end
            close all;
            if strcmp(obj.Geometry.DiscretizationType,'CartesianGrid')
                obj.Geometry.PlotFractures();
            end
        end
        function OutputResults(obj)
            warning("VTK FILE WRITER IS DEFINED HERE")
            switch obj.Geometry.DiscretizationType
                case('CartesianGrid')
                    obj.Writer.VTK_Plotter = plotVTK_CarstesianGrid_FracGen(strcat(obj.Writer.Directory,'/Output_FracGen'), '3D-pEDFM');
                case('CornerPointGrid')
                    obj.Writer.VTK_Plotter = plotVTK_CornerPointGrid_FracGen(strcat(obj.Writer.Directory,'/Output_FracGen'), '3D-pEDFM');
                case('TetrahedralGrid')
                    obj.Writer.VTK_Plotter = plotVTK_TetrahedralGrid_FracGen(strcat(obj.Writer.Directory,'/Output_FracGen'), '3D-pEDFM');
            end
            obj.Writer.VTK_Plotter.PlotSolution(obj.Geometry);
            %%
            warning("text file is generated here - summary of the data structure and all the information")
            obj.Writer.PrintInfo();
            obj.Writer.WriteOutputTextFile(obj.Geometry);
        end
    end
end
