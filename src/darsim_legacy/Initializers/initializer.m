% Initializer base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef initializer < handle
    properties
        VarNames
        VarValues
    end
    methods
        function obj = initializer(names, values)
            obj.VarNames = names;
            obj.VarValues = values;
        end
        function InitializeProductionSystem(obj, ProductionSystem, FluidModel, Formulation, DiscretizationModel)
            % Initialise state           
            % 1. Assign initial values 
            ProductionSystem.AssignInitialState(obj.VarNames, obj.VarValues);
            % Initialise capillary model
            FluidModel.CapillaryModel.Initialise(ProductionSystem);
            
            % 2. Compute initial phase and component distribution
            ProductionSystem.Reservoir.State.T = ProductionSystem.Reservoir.Temp; % For now it's fine like this
            if ProductionSystem.FracturesNetwork.Active
               for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                   ProductionSystem.FracturesNetwork.Fractures(f).State.T = ProductionSystem.Reservoir.Temp; % Using reservoir temperature for fractures as well
               end
            end
            
            % 3. Define initial state
            start = tic;
            obj.ComputeInitialState(ProductionSystem, FluidModel, Formulation, DiscretizationModel);
            initialization = toc(start);
            disp(['Initialization took ', num2str(initialization), ' s']);
            disp(newline);
        end
        function ComputeInitialState(obj, ProductionSystem, FluidModel, Formulation, DiscretizationModel)
            disp('Started simple initialization');
            %% Reservoir
            % 2. Update Composition of the phases (Flash)
            Index.Start = 1;
            Index.End = DiscretizationModel.ReservoirGrid.N;           
            
            % Initialize K-values based on Wilson's approximation for 3-comp 
            if FluidModel.NofComponents == 3
                FluidModel.FlashCalculator.K_Wilson(ProductionSystem.Reservoir.State, FluidModel.Thermo)
            end

            % To achive a stable initial conditions we do this in a loop because:
            % {P_1, z_i} > {V, x_i, y_i} > {rho_i} > {S_i} > {P_c} > {new P_1}
            P_0 = ProductionSystem.Reservoir.State.Properties('P_1').Value;
            equilibrium = false;
            while ~equilibrium
                switch FluidModel.name
                    case({'BlackOil','Compositional','Compositional_Energy'})
                        SinglePhase = FluidModel.Flash(ProductionSystem.Reservoir.State);
                        status_SinglePhase = ProductionSystem.Reservoir.State.Properties('SinglePhase');
                        status_SinglePhase.Value = SinglePhase;
                    otherwise
                        SinglePhase = FluidModel.Flash(ProductionSystem.Reservoir.State);
                end
                
                % 3 Compute Phase Density
                FluidModel.ComputePhaseDensities(ProductionSystem.Reservoir.State, SinglePhase);
                
                % 4. Update S based on components mass balance
                FluidModel.ComputePhaseSaturation(ProductionSystem.Reservoir.State, SinglePhase);
                
                % 6. Compute initial Pc
                FluidModel.ComputePc(ProductionSystem.Reservoir, DiscretizationModel);

                % 5. Total Density
                FluidModel.ComputeTotalDensity(ProductionSystem.Reservoir.State);

                if strcmpi(FluidModel.name,'Compositional_Energy')
                    
                    % 7. Compute Hr (Hr = Cpr*T)
                    FluidModel.ComputeRockEnthalpy(ProductionSystem.Reservoir);

                    % 8. Compute Deff
                    FluidModel.ComputeThermalConductivity(ProductionSystem.Reservoir);

                    % 9. Compute Enthalpy (H = H(p, T))
                    FluidModel.ComputePhaseEnthalpies(ProductionSystem.Reservoir.State);
                    
                end

                % Check equilibrium
                P = ProductionSystem.Reservoir.State.Properties('P_1').Value;
                if norm(P - P_0, Inf) < 1e-6
                    equilibrium = true;
                else
                    P_0 = P;
                end
            end
            
            % Output initial status:      
            disp('Initial conditions:')
            disp(['Pressure: ' num2str(max(ProductionSystem.Reservoir.State.Properties('P_2').Value/1e5)), ' bar']);
            disp(['Saturation: ' num2str(max(ProductionSystem.Reservoir.State.Properties('S_1').Value))]);
            disp(['Temperature: ', num2str(ProductionSystem.Reservoir.Temp), ' K']);
            disp('---------------------------------------------------------');
            %% Fractures
            if ProductionSystem.FracturesNetwork.Active
                for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                    Index.Start = DiscretizationModel.Index_Local_to_Global(DiscretizationModel.ReservoirGrid.N, f, 1);
                    Index.End = DiscretizationModel.Index_Local_to_Global(DiscretizationModel.ReservoirGrid.N, f, DiscretizationModel.FracturesGrid.Grids(f).N);

                    switch FluidModel.name
                        case({'BlackOil','Compositional','Compositional_Energy'})
                            SinglePhase = FluidModel.Flash(ProductionSystem.FracturesNetwork.Fractures(f).State);
                            status_SinglePhase = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('SinglePhase');
                            status_SinglePhase.Value = SinglePhase;
                        otherwise
                            SinglePhase = FluidModel.Flash(ProductionSystem.FracturesNetwork.Fractures(f).State);
                    end

                    % 3 Compute Phase Density
                    FluidModel.ComputePhaseDensities(ProductionSystem.FracturesNetwork.Fractures(f).State, SinglePhase);

                    % 4. Update S based on components mass balance
                    FluidModel.ComputePhaseSaturation(ProductionSystem.FracturesNetwork.Fractures(f).State, SinglePhase);

                    % 6. Compute initial Pc
                    FluidModel.ComputePc(ProductionSystem.FracturesNetwork.Fractures(f), DiscretizationModel);
                    % FluidModel.ComputePcHysteresis(ProductionSystem.FracturesNetwork.Fractures(f).State);

                    % 5. Total Density
                    FluidModel.ComputeTotalDensity(ProductionSystem.FracturesNetwork.Fractures(f).State);
                    
                    if strcmpi(FluidModel.name,'Compositional_Energy')

                        % 7. Compute Hr (Hr = Cpr*T)
                        FluidModel.ComputeRockEnthalpy(ProductionSystem.FracturesNetwork.Fractures(f));

                        % 8. Compute Deff
                        FluidModel.ComputeThermalConductivity(ProductionSystem.FracturesNetwork.Fractures(f));

                        % 9. Compute Enthalpy (H = H(p, T))
                        FluidModel.ComputePhaseEnthalpies(ProductionSystem.FracturesNetwork.Fractures(f).State);

                    end
                    
                end
            end
            disp(newline);
        end
    end
end

