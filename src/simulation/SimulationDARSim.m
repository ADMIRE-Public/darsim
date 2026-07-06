%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel;
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef SimulationDARSim < handle
    properties
        Initializer
        ProductionSystem
        FluidModel
        DiscretizationModel
        Formulation
        TimeDriver
        Summary
        Writer
    end
    properties%(Access = protected)
        SimulatorSettings
        SimulatorFolder
        SimulationTitle        
    end
    methods
        function prepare_simulation(obj)
            obj.define_simulation_properties_to_production_system();
            obj.DiscretizationModel.Initialize(obj.ProductionSystem, ...
                obj.Formulation);
            obj.Initializer.InitializeProductionSystem(...
                obj.ProductionSystem, obj.FluidModel, ...
                obj.Formulation, obj.DiscretizationModel);
            obj.ProductionSystem.InitializeWells(obj.Formulation, ...
                obj.FluidModel, obj.Formulation.GravityModel, ...
                obj.DiscretizationModel);
            obj.DiscretizationModel.InitializeMapping(...
                obj.ProductionSystem, obj.FluidModel);
            % obj.DiscretizationModel.AddWellsToInitialPressure(
            % obj.ProductionSystem, obj.FluidModel);
            % For multiscale to avoid having well functions
            % It's not working as it should.
        end
        function set.SimulatorSettings(obj, settings)
            obj.SimulatorSettings = settings;
        end
        function set.SimulationTitle(obj, title)
            obj.SimulationTitle = title;
        end        
        function Run(obj)
            disp('BEGIN TIME-DEPENDENT SIMULATION');
            obj.WriteSimulationInitialSate();
            obj.PrintInfo();
            TotalStart = tic;
            obj.TimeDriver.SolveTimeDependentProblem(...
                obj.ProductionSystem, obj.FluidModel, ...
                obj.DiscretizationModel, obj.Formulation, obj.Summary, ...
                obj.Writer);
            obj.OutputResults();
            TotalTime = toc(TotalStart);
            display_total_time(TotalTime,TotalStart, obj.TimeDriver.Coupling.Converged);
        end
    end
    methods (Access = private)
        function OutputResults(obj)
            obj.Writer.WriteSummary(obj.Summary);
        end
        function PrintInfo(obj)
            isFractured = ~isempty(obj.DiscretizationModel.FracturesGrid);
            fprintf('SIMULATION OF PROBLEM %s HAS STARTED.\n', ...
                obj.SimulationTitle);
            disp('FORMATION CHARACTERISTICS:');
            disp('Reservoir geometry:');
            disp(['Lx: ', ...
                num2str(obj.ProductionSystem.Reservoir.Length), ' m']);
            disp(['Ly: ', ...
                num2str(obj.ProductionSystem.Reservoir.Width), ' m']);
            disp(['Thickness:  ', ...
                num2str(obj.ProductionSystem.Reservoir.Thickness), ' m']);
            if ~strcmp(obj.DiscretizationModel.type,'TetrahedralGrid') ...
               && ~strcmp(obj.DiscretizationModel.type,'UnstructuredGrid')
                Nx = obj.DiscretizationModel.ReservoirGrid.Nx;
                Ny = obj.DiscretizationModel.ReservoirGrid.Ny;
                Nz = obj.DiscretizationModel.ReservoirGrid.Nz;
                disp(['Total Grid Cells: ', num2str(Nx), ' x ', ...
                    num2str(Ny), ' x ', num2str(Nz), ' = ', ...
                    num2str(Nx*Ny*Nz), ', Number of Active Cells: ',...
                    num2str(obj.DiscretizationModel.ReservoirGrid.N)]);
            else
                disp(['Total Grid Cells: ', ...
                    num2str(obj.DiscretizationModel.ReservoirGrid.N),...
                    ', Number of Active Cells: ', ...
                    num2str(obj.DiscretizationModel.ReservoirGrid.N)]);
            end
            disp(['-------------------------' ...
                '--------------------------------']);
            if isFractured
                if ~strcmp(obj.DiscretizationModel.type,'UnstructuredGrid')
                    for f = 1:obj.DiscretizationModel.FracturesGrid.Nfrac
                        Nx = obj.DiscretizationModel. ...
                            FracturesGrid.Grids(f).Nx;
                        Ny = obj.DiscretizationModel. ...
                            FracturesGrid.Grids(f).Ny;
                        fprintf(['Fracture %2d: Grid= %3.0f x' ...
                            ' %3.0f = %3.0f\n'], f, Nx, Ny, Nx*Ny);
                    end
                else
                    for f = 1:obj.DiscretizationModel.FracturesGrid.Nfrac
                        N = obj.DiscretizationModel. ...
                            FracturesGrid.Grids(f).N;
                        fprintf(['Fracture %2d: Grid= %3.0f Cells: ' ...
                            '%3.0f\n'], f, N);
                    end
                end
                fprintf('Total fracture grids: %3.0f\n', ...
                    sum(obj.DiscretizationModel.FracturesGrid.N));
                disp(['-------------------------' ...
                    '--------------------------------']);
            end
            disp('Fluid Model:');
            disp(['Type: ', obj.FluidModel.name]);
            disp(['N of phases: ', num2str(obj.FluidModel.NofPhases)]);
            disp(['N of components:  ',...
                num2str(obj.FluidModel.NofComponents)]);
            disp(['-------------------------' ...
                '--------------------------------']);

            Time = obj.TimeDriver.Sec2DHMS(obj.TimeDriver.TotalTime);
            disp(['Total Simulation Time: '...
                num2str(obj.TimeDriver.TotalTime) ' seconds (', ...
                num2str(Time.Days), ' days : ', ...
                num2str(Time.Hours), ' hrs : ', ...
                num2str(Time.Minutes), ' mins : ', ...
                num2str(Time.Seconds), ' sec)']);
            disp(['Maximum Allowed Time-steps: ', ...
                num2str(obj.SimulatorSettings.MaxNumTimeSteps) ]);
            disp(['Minimum Size of Time-step : ', ...
                num2str(obj.SimulatorSettings.MinMaxdt(1)) ]);
            disp(['Maximum Size of Time-step : ', ...
                num2str(obj.SimulatorSettings.MinMaxdt(2)) ]);
            disp(['Number of Reports         : ', ...
                num2str(obj.SimulatorSettings.Reports)     ]);
            disp(['-------------------------' ...
                '--------------------------------']);
        end
        function WriteSimulationInitialSate(obj)
            obj.Writer.PlotSolution(obj.ProductionSystem, ...
                obj.DiscretizationModel);
            % Write initial state on a file
            obj.Writer.WriteSolutionOnFile(obj.ProductionSystem, obj.FluidModel)
            obj.Writer.Index = obj.Writer.Index + 1;
            % Plot Wells
            obj.Writer.Plotter.PlotWells(obj.ProductionSystem.Wells.Inj, ...
                obj.ProductionSystem.Wells.Prod, ...
                obj.DiscretizationModel.ReservoirGrid);
        end

        function obj = ...
                define_simulation_properties_to_production_system(obj)
            if ~obj.ProductionSystem.FracturesNetwork.Active
                obj.ProductionSystem.Reservoir.State.AddProperties(...
                  obj.FluidModel, obj.DiscretizationModel.ReservoirGrid.N);
                obj.ProductionSystem.Reservoir.State_old.AddProperties(...
                  obj.FluidModel, obj.DiscretizationModel.ReservoirGrid.N);
            else
                obj.ProductionSystem.Reservoir.State.AddProperties(...
                  obj.FluidModel, obj.DiscretizationModel.ReservoirGrid.N);
                obj.ProductionSystem.Reservoir.State_old.AddProperties(...
                  obj.FluidModel, obj.DiscretizationModel.ReservoirGrid.N);
                for f=1:obj.ProductionSystem.FracturesNetwork.NumOfFrac
                    obj.ProductionSystem.FracturesNetwork. ...
                        Fractures(f).State.AddProperties(obj.FluidModel,...
                        obj.DiscretizationModel.FracturesGrid.N(f));
                    obj.ProductionSystem.FracturesNetwork.Fractures(f). ...
                        State_old.AddProperties(obj.FluidModel, ...
                        obj.DiscretizationModel.FracturesGrid.N(f));
                end
            end
        end
    end
end

function display_total_time(TotalTime, TotalStart, Converged)
%% Display elapsed time
disp(newline);
disp(['The Total Simulation time is '...
            num2str(TotalTime) ' [s]']);
%disp('end of the simulation');

if Converged ~= 0
    disp('Simulation successfully finished!');
end
end

