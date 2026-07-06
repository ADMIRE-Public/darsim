%  FIM System Builder base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 19 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef fim_system_builder < system_builder
    properties
        NumberOfEq  
    end
    methods
        function ComputePropertiesAndDerivatives(obj, Formulation, ProductionSystem, FluidModel, DiscretizationModel)
            Formulation.ComputePropertiesAndDerivatives(ProductionSystem, FluidModel);
            Formulation.UpWindAndPhaseRockFluxes(DiscretizationModel, FluidModel.Phases, ProductionSystem);
        end
        function [Residual, RHS] = BuildResidual(obj, ProductionSystem, DiscretizationModel, Formulation, Time, dt)
           [Residual, RHS] = Formulation.BuildFullResidual(ProductionSystem, DiscretizationModel, Time, dt, obj.State);
           % [Residual, RHS] = Formulation.BuildFullResidual(ProductionSystem, DiscretizationModel, Time, dt, ProductionSystem.Reservoir.State_old);
        end
        function Jacobian = BuildJacobian(obj, ProductionSystem, Formulation, DiscretizationModel, Time, dt)
            Jacobian = Formulation.BuildFullJacobian(ProductionSystem, DiscretizationModel, Time, dt);
        end
        function SetUpSolutionChopper(obj, SolutionChopper, Formulation, ProductionSystem, DiscretizationModel)
            x = Formulation.GetPrimaryUnknowns(ProductionSystem, DiscretizationModel);
            SolutionChopper.DefineMaxDelta(x);
        end
        function delta = UpdateState(obj, delta, ProductionSystem, Formulation, FluidModel, DiscretizationModel, Time)
            % Update State
            delta = Formulation.UpdateState(delta, ProductionSystem, FluidModel, DiscretizationModel);
            
            % Update Wells
            % ProductionSystem.Wells.UpdateState(ProductionSystem.Reservoir, Formulation, FluidModel);
            ProductionSystem.Wells.UpdateState(ProductionSystem, Formulation, FluidModel, Time);
        end
    end
end