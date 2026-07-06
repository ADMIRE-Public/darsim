%  Pressure System Builder base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef pressure_system_builder < system_builder
    properties
    end
    methods 
        function ComputePropertiesAndDerivatives(obj, Formulation, ProductionSystem, FluidModel, DiscretizationModel)
            Formulation.ComputePropertiesAndDerivatives(ProductionSystem, FluidModel);
            Formulation.UpWindAndPhaseRockFluxes(DiscretizationModel, FluidModel.Phases, ProductionSystem);
        end
        function [Residual, RHS] = BuildResidual(obj, ProductionSystem, DiscretizationModel, Formulation, Time, dt)
           [Residual, RHS] = Formulation.BuildFullResidual(ProductionSystem, DiscretizationModel, Time, dt, obj.State);
        end
        function Jacobian = BuildJacobian(obj, ProductionSystem, Formulation, DiscretizationModel, Time,dt)
            Jacobian = Formulation.BuildPressureMatrix(ProductionSystem, DiscretizationModel, dt, obj.State);
        end
        function SetUpSolutionChopper(obj, SolutionChopper, Formulation, ProductionSystem, N)
            x = Formulation.GetPrimaryPressure(ProductionSystem, N);
            SolutionChopper.DefineMaxDelta(x);
        end
        function delta = UpdateState(obj, delta, ProductionSystem, Formulation, FluidModel, DiscretizationModel)
            % Update Reservoir State
            Formulation.UpdatePressure(delta, ProductionSystem, FluidModel, DiscretizationModel);
            % UpdateWells
            % ProductionSystem.Wells.UpdateState(ProductionSystem.Reservoir, FluidModel);
            ProductionSystem.Wells.UpdateState(ProductionSystem, Formulation, FluidModel);
        end
    end
end