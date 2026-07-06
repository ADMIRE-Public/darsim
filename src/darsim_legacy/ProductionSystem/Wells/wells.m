% Wells 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef wells < handle
    properties
        Inj
        Prod
        NofInj
        NofProd
    end
    methods
        function AddInjector(obj, Injector)
            obj.Inj{end+1} = Injector;
        end
        function AddProducer(obj, Producer)
            obj.Prod{end+1} = Producer;
        end
        function InitializeFluxes(obj, n_phases, n_comp)
            % Injectors
            for i=1:obj.NofInj
                obj.Inj{i}.QPhases = zeros(length(obj.Inj{i}.Cells), n_phases);
                obj.Inj{i}.Qh = zeros(length(obj.Inj{i}.Cells), n_phases);
                obj.Inj{i}.QComponents = zeros(length(obj.Inj{i}.Cells), n_comp);
            end
            % Producers
            for i=1:obj.NofProd
                obj.Prod{i}.QPhases = zeros(length(obj.Prod{i}.Cells), n_phases);
                obj.Prod{i}.Qh = zeros(length(obj.Prod{i}.Cells), n_phases);
                obj.Prod{i}.QComponents = zeros(length(obj.Prod{i}.Cells), n_comp);
            end
        end
        function UpdateState(obj, ProductionSystem, Formulation, FluidModel, Time)
            K = ProductionSystem.Reservoir.K(:,1);
            switch FluidModel.name 
                case {'Geothermal_SinglePhase','Geothermal_MultiPhase'} 
                    Formulation.Mob = FluidModel.ComputePhaseMobilities(ProductionSystem.Reservoir.State);
                otherwise
                    Formulation.Mob = FluidModel.ComputePhaseMobilities(ProductionSystem.Reservoir);
                    % if FluidModel.HysteresisActive
                    %     Formulation.Mob = FluidModel.ComputePhaseMobilitiesHysteresis(ProductionSystem.Reservoir);
                    % else
                    %     Formulation.Mob = FluidModel.ComputePhaseMobilities(ProductionSystem.Reservoir);
                    % end
            end
            
            % Injectors
            for i=1:obj.NofInj
                obj.Inj{i}.UpdateState(ProductionSystem.Reservoir.State, K, FluidModel, Time);         
            end
            % Producers
            for i=1:obj.NofProd
                obj.Prod{i}.UpdateState(ProductionSystem.Reservoir.State, K, Formulation.Mob, FluidModel, Time);        
            end
            
            % For fractures
            for f = 1:ProductionSystem.FracturesNetwork.NumOfFrac
                Formulation.Mob = [Formulation.Mob; FluidModel.ComputePhaseMobilities(ProductionSystem.FracturesNetwork.Fractures(f))];
            end            
        end
        function q = TotalFluxes(obj, Reservoir, Mobt)
            K = Reservoir.K(:,1);
            q = zeros(length(K), 1);
            
            % Injectors
            for i=1:obj.NofInj
                q = obj.Inj{i}.TotalFlux(q, Reservoir.State.p, K);         
            end
            % Producers
            for i=1:obj.NofProd
                q = obj.Prod{i}.TotalFlux(q, Reservoir.State.p, K, Mobt);        
            end
        end
    end
end