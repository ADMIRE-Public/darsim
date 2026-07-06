% FIM coupling strategy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef FIM_Strategy < Coupling_Strategy
    properties
        NLSolver
        chops
        MaxChops
        CFL
    end
    methods
        function obj = FIM_Strategy(name, NONLinearSolver)            
            obj@Coupling_Strategy(name);
            obj.NLSolver = NONLinearSolver;
            obj.chops = 0;
            obj.MaxChops = 15;
        end
        function [dt, End] = SolveTimeStep(obj, ProductionSystem, FluidModel, DiscretizationModel, Formulation, Time, Ndt)
            
            % Initialise
            dt = obj.TimeStepSelector.ChooseTimeStep(Time);
            obj.Converged = 0;
            obj.chops = 0;
            End = 0;
            
            DT = obj.TimeStepSelector.Sec2DHMS(dt);
            disp(['dt= ' num2str(dt) ' sec (', num2str(DT.Days), ' days : ', num2str(DT.Hours), ' hrs : ', num2str(DT.Minutes), ' mins : ', num2str(DT.Seconds), ' sec)']);
            
            % Save state of current time-step (it's useful for ADM to update based on time change)
            if (Ndt == 1)
                ProductionSystem.Reservoir.State_old_old.CopyProperties(ProductionSystem.Reservoir.State);
                ProductionSystem.Reservoir.State_old.CopyProperties(ProductionSystem.Reservoir.State);
                
                for f = 1:ProductionSystem.FracturesNetwork.NumOfFrac
                    ProductionSystem.FracturesNetwork.Fractures(f).State_old_old.CopyProperties(ProductionSystem.FracturesNetwork.Fractures(f).State);
                    ProductionSystem.FracturesNetwork.Fractures(f).State_old.CopyProperties(ProductionSystem.FracturesNetwork.Fractures(f).State);
                end
                
            else
                ProductionSystem.SavePreviousState(); % obj.Medium.State_old & obj.Medium.State_old_old
            end
            if FluidModel.HysteresisActive
                FluidModel.RelPermModel.ComputeTransportLabel(FluidModel.Phases, ProductionSystem.Reservoir); % determine the history for hysteresis
            end
            % Set Up non-linear solver
            % useless line
           obj.NLSolver.SetUpLinearSolver(ProductionSystem, DiscretizationModel);
           obj.NLSolver.SetUp(Formulation, ProductionSystem, FluidModel, DiscretizationModel, dt);
            
            while (obj.Converged == 0 && obj.chops < obj.MaxChops)
                % Print some info to the screen
                if (obj.chops > 0)
                    disp('Max num. of iterations reached or stagnation detected: Time-step was chopped');
                    disp(newline);
                    disp(['Restart Newton loop dt = ', num2str(dt)]);
                end
                
                % NL solver call
                obj.NLSolver.Solve(ProductionSystem, FluidModel, DiscretizationModel, Formulation, Time, dt);
                
                % Chop time-step if it failed to converge
                if obj.NLSolver.Converged == 0
                    dt = dt/2;
                    obj.chops = obj.chops + 1;
                    Formulation.Reset();
                    obj.NLSolver.SystemBuilder.SetInitalGuess(ProductionSystem); % copy SystemBuilder.status to ProductionSystem.Reservoir.State
                    FluidModel.RelPermModel.ResetRelPermModel(ProductionSystem.Reservoir.State); % reset relative permeability model
                    % FluidModel.ComputePhaseMobilitiesHysteresis(ProductionSystem.Reservoir); % recompute phase mobility
                    ProductionSystem.Wells.UpdateState(ProductionSystem, Formulation, FluidModel, Time);
                    obj.NLSolver.SetUp(Formulation, ProductionSystem, FluidModel, DiscretizationModel, dt); % copy ProductionSystem.Reservoir.State to SystemBuilder.status 
                end
                obj.Converged = obj.NLSolver.Converged;
            end
            Formulation.ComputeTotalFluxes(ProductionSystem, DiscretizationModel);
            % Formulation.CheckMassConservation(DiscretizationModel);
            obj.CFL = Formulation.ComputeCFLNumber(ProductionSystem, DiscretizationModel, dt);
            % disp(['CFL = ', num2str(obj.CFL)]);
            obj.TimeStepSelector.Update(dt, obj.NLSolver.itCount - 1, obj.chops);
        end
        function Summary = UpdateSummary(obj, Summary, Wells, Ndt, dt)
            %% Stats, timers and Injection/Production data
            Summary.CouplingStats.SaveStats(Ndt, dt, obj.NLSolver.itCount-1, obj.chops, obj.CFL);
            Summary.CouplingStats.SaveTimers(Ndt, obj.NLSolver.TimerConstruct, obj.NLSolver.TimerSolve, obj.NLSolver.TimerInner);
            Summary.SaveWellsData(Ndt+1, Wells.Inj, Wells.Prod, dt);
        end
    end
end