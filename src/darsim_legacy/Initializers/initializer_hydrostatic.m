% Initializer hydrostatic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef initializer_hydrostatic < initializer
    properties
        NLSolver
    end
    methods
        function obj = initializer_hydrostatic(names, values)
            obj@initializer(names, values);
        end
        function ComputeInitialState(obj, ProductionSystem, FluidModel, Formulation, DiscretizationModel)
            disp('Started Hydrostatic initialization');

            obj.Equilibrate(ProductionSystem, FluidModel, Formulation, DiscretizationModel);
            
             % Output initial status:      
            disp('Initial conditions:')
            disp(['reservoir pressure:' num2str(max(ProductionSystem.Reservoir.State.Properties('P_2').Value/1e5)), ' bar']);
            disp(['reservoir saturation:' num2str(max(ProductionSystem.Reservoir.State.Properties('S_1').Value))]);
            disp(['reservoir temperature: ', num2str(ProductionSystem.Reservoir.Temp)]);
            disp('---------------------------------------------------------');
            disp(newline);
        end
        function Equilibrate(obj, ProductionSystem, FluidModel, Formulation, DiscretizationModel)
            equilibrium = 0;
            % Save initial State
            %obj.NLSolver.SystemBuilder.SaveInitialState(ProductionSystem.Reservoir.State, Formulation);
            P_init = ProductionSystem.Reservoir.State.Properties('P_2').Value;
            P = ProductionSystem.Reservoir.State.Properties('P_2');
            
            % status_P_1 = ProductionSystem.Reservoir.State.Properties('P_1');
            % status_P_1.Value = P_init;
            
            while (~equilibrium)
                p0 = ProductionSystem.Reservoir.State.Properties('P_2').Value;
       
                % 1. Update Composition of the phases (Flash)
                Formulation.SinglePhase = FluidModel.Flash(ProductionSystem.Reservoir.State);
                status_SinglePhase = ProductionSystem.Reservoir.State.Properties('SinglePhase');
                status_SinglePhase.Value = Formulation.SinglePhase;
                
                % 2 Compute Phase Density
                FluidModel.ComputePhaseDensities(ProductionSystem.Reservoir.State, Formulation.SinglePhase);
                
                % 3. Update S based on components mass balance
                FluidModel.ComputePhaseSaturation(ProductionSystem.Reservoir.State, Formulation.SinglePhase);
                
                % 4. Total Density
                FluidModel.ComputeTotalDensity(ProductionSystem.Reservoir.State);
                
                % 5. Compute initial Pc
                FluidModel.ComputePc(ProductionSystem.Reservoir.State, DiscretizationModel, 1);
                
                % 6. Compute pressure 
                P.Value = P_init + ProductionSystem.Reservoir.State.Properties('rhoT').Value .* Formulation.GravityModel.g .* DiscretizationModel.ReservoirGrid.Depth;
                delta = ProductionSystem.Reservoir.State.Properties('P_2').Value - p0;

                if strcmpi(FluidModel.name,'Compositional_Energy')
                    
                    % 7. Compute Hr (Hr = Cpr*T)
                    FluidModel.ComputeRockEnthalpy(ProductionSystem.Reservoir);

                    % 8. Compute Deff
                    FluidModel.ComputeThermalConductivity(ProductionSystem.Reservoir);

                    % 9. Compute Enthalpy (H = H(p, T))
                    FluidModel.ComputePhaseEnthalpies(ProductionSystem.Reservoir.State);
                    
                end
                
                if (max(abs(delta)) < 1e-3) % If solution doesn t change I am at equilibrium
                    equilibrium = 1;
                end
            end
        end
    end
end

