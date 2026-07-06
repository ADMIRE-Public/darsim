% Output writer MMs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Mousa HosseiniMehr
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef output_writer_MMs < output_writer_FS
    properties
        PlotBasisFunctions
    end
    methods
        function obj = output_writer_MMs(dir, problem, n_inj, n_prod, n_timers,  n_stats, n_previous_reports, n_comp)
            obj@output_writer_FS(dir, problem, n_inj, n_prod, n_timers,  n_stats, n_previous_reports, n_comp);
            obj.PlotBasisFunctions = false;
        end
        function PlotSolution(obj, ProductionSystem, DiscretizationModel)
            obj.Plotter.PlotSolution(ProductionSystem, DiscretizationModel);
            obj.Plotter.PlotPermeability(ProductionSystem, DiscretizationModel);
            obj.Plotter.PlotPorosity(ProductionSystem, DiscretizationModel);
            obj.Plotter.PlotADMGrid(ProductionSystem, DiscretizationModel);
            if obj.PlotBasisFunctions
                obj.Plotter.PlotBasisFunctions(DiscretizationModel.FineGrid, ...
                    DiscretizationModel.CoarseGrid, DiscretizationModel.OperatorsHandler.ProlongationBuilders(1).P, DiscretizationModel.Nf, DiscretizationModel.Nc);
                obj.PlotBasisFunctions = false;
            end
            obj.Plotter.VTKindex = obj.Plotter.VTKindex + 1;
        end
    end
end
       