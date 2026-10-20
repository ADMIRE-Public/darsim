%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel;
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef simulation_factory < handle
    properties
        SimulationInput
        SimulatorSettings
        SimulatorFolder
        incompressible = 0;
        capillarymodel = 0;
        NofEq
    end
    properties (Access = private)
        Directory
    end
    methods
        function obj = simulation_factory(darsim_settings)
                obj.SimulationInput = darsim_settings.SimulationInput;
                obj.SimulatorSettings = darsim_settings.SimulatorSettings;
                obj.SimulatorFolder = darsim_settings.SimulatorFolder;
        end
        function simulation = BuildSimulation(obj)
            simulation = SimulationDARSim();
            simulation.SimulatorFolder = obj.SimulatorFolder;
            %FractureMatrix = obj.SimulationInput.FracturesProperties;
            simulation.DiscretizationModel = obj.BuildDiscretization(); %
            % This line create the objects that are the grids with their properties.
            simulation.ProductionSystem = ...
                obj.BuildProductionSystemReservoir();
            simulation.FluidModel = obj.BuildFluidModel();
            simulation.Formulation = obj.BuildFormulation();
            simulation.TimeDriver = obj.BuildTimeDriver();
            simulation.Summary = obj.BuildSummary(simulation);
            simulation.Writer = build_visualization_tool(...
                                obj.SimulationInput,obj.SimulatorSettings); 
            simulation.SimulatorSettings = obj.SimulatorSettings;
            simulation.SimulationTitle = obj.SimulationInput.ProblemName;
            n_frac = obj.SimulationInput.FracturesProperties.NumOfFrac;      
            % Add gravity model
            simulation.Formulation.GravityModel = gravity_model(...
                simulation.DiscretizationModel, ...
                simulation.FluidModel.NofPhases, n_frac);
            if obj.SimulationInput.FluidProperties.Gravity
                simulation.Formulation.GravityModel.g = 9.806;
            else
                simulation.Formulation.GravityModel.g = 0;
            end            
            % Add hysteresis
            simulation.FluidModel.HysteresisActive = ...
                obj.SimulationInput.FluidProperties.Hysteresis;
            simulation.Initializer = ...
                build_simulator_InitalizerProcedure(obj, ...
                simulation.FluidModel.name, ...
                simulation.DiscretizationModel.N);            
        end 

        function Initializer = build_simulator_InitalizerProcedure(obj, ...
                                                            FluidName, N)
            %% Define Initialization procedure
            % Set the Variable Names and the Initialization Routine
            % Select FluidModel based on FluidModel keyword
            switch(FluidName)
                case('SinglePhase')
                    VarNames = {'P_1', 'S_1'};
                    if obj.SimulationInput.InitialConditions.Include
                        VarValues = ones(N,2) .* ...
                         obj.SimulationInput.InitialConditions.LoadedData;
                    else
                        VarValues = [ ones(N,1).* ...
                         obj.SimulationInput.InitialConditions.Pressure , ...
                         ones(N,1).*1];
                    end
                    Initializer = initializer_singlephase(VarNames, VarValues);
                case('Immiscible')
                    VarNames = {'P_2', 'S_1', 'S_2'};
                    if obj.SimulationInput.InitialConditions.Include
                        VarValues = ones(N,3) .* ...
                          obj.SimulationInput.InitialConditions.LoadedData;
                    else
                        VarValues = [ ones(N,1).* ...
                            obj.SimulationInput.InitialConditions.Pressure, ...
                            ones(N,1).*obj.SimulationInput.InitialConditions.Saturation_1 , ...
                            ones(N,1).*obj.SimulationInput.InitialConditions.Saturation_2   ];
                    end
                    Initializer = initializer(VarNames, VarValues);
                case('Geothermal_SinglePhase')
                    VarNames = {'P_1', 'hTfluid', 'T', 'S_1'};
                    if obj.SimulationInput.InitialConditions.Include
                        VarValues = ones(N,4) .* obj.SimulationInput.InitialConditions.LoadedData;
                    else
                        VarValues = [ ones(N,1).*obj.SimulationInput.InitialConditions.Pressure    , ...
                            ones(N,1).*obj.SimulationInput.InitialConditions.Enthalpy    , ...
                            ones(N,1).*obj.SimulationInput.InitialConditions.Temperature , ...
                            ones(N,1).*1                                                   ];
                    end
                    Initializer = Geothermal_SinglePhase_initializer(VarNames, VarValues);
                case('Geothermal_MultiPhase')
                    VarNames = {'P_2', 'hTfluid', 'T'};
                    if obj.SimulationInput.InitialConditions.Include
                        VarValues = ones(N,3) .* obj.SimulationInput.InitialConditions.LoadedData;
                    else
                        VarValues = [ ones(N,1).*obj.SimulationInput.InitialConditions.Pressure , ...
                            ones(N,1).*obj.SimulationInput.InitialConditions.Enthalpy , ...
                            ones(N,1).*obj.SimulationInput.InitialConditions.Temperature];
                    end
                    Initializer = Geothermal_MultiPhase_initializer(VarNames, VarValues);
                case{'Compositional','BlackOil'}
                    nc = obj.SimulationInput.FluidProperties.NofComponents;
                    VarNames = "P_2";
                    for i = 1:nc
                        VarNames = [VarNames, sprintf('z_%d', i)];
                    end

                    if obj.SimulationInput.InitialConditions.Include
                        VarValues = ones(N, nc+1) .* obj.SimulationInput.InitialConditions.LoadedData;
                    else
                    if nc == 2
                        VarValues = [obj.SimulationInput.InitialConditions.Pressure , ...
                                     obj.SimulationInput.InitialConditions.Composition_1 , ...
                                     obj.SimulationInput.InitialConditions.Composition_2];
                        VarValues = ones(N, 1) * VarValues;
                    elseif nc == 3
                        VarValues = [obj.SimulationInput.InitialConditions.Pressure, ...
                                     obj.SimulationInput.InitialConditions.Composition_1, ...
                                     obj.SimulationInput.InitialConditions.Composition_2, ...
                                     obj.SimulationInput.InitialConditions.Composition_3];
                        VarValues = ones(N, 1) * VarValues;
                    else
                        error('Max comp. number is 3.')
                    end
                    end
                    % To not have zero gas pressure
                    VarNames = [VarNames, "P_1"];
                    VarValues = [VarValues, VarValues(:, 1)]; 
                    Initializer = initializer(VarNames, VarValues); % simple versoin
                case{'Compositional_Energy'}
                    nc = obj.SimulationInput.FluidProperties.NofComponents;
                    VarNames = ["P_2", "T"];
                    for i = 1:nc
                        VarNames = [VarNames, sprintf('z_%d', i)];
                    end

                    if obj.SimulationInput.InitialConditions.Include
                        VarValues = ones(N, nc+1) .* obj.SimulationInput.InitialConditions.LoadedData;
                    else
                    if nc == 2
                        VarValues = [obj.SimulationInput.InitialConditions.Pressure , ...
                                     obj.SimulationInput.InitialConditions.Temperature, ...
                                     obj.SimulationInput.InitialConditions.Composition_1 , ...
                                     obj.SimulationInput.InitialConditions.Composition_2];
                        VarValues = ones(N, 1) * VarValues;
                    elseif nc == 3
                        VarValues = [obj.SimulationInput.InitialConditions.Pressure, ...
                                     obj.SimulationInput.InitialConditions.Temperature, ...
                                     obj.SimulationInput.InitialConditions.Composition_1, ...
                                     obj.SimulationInput.InitialConditions.Composition_2, ...
                                     obj.SimulationInput.InitialConditions.Composition_3];
                        VarValues = ones(N, 1) * VarValues;
                    else
                        error('Max comp. number is 3.')
                    end
                    end
                    % To not have zero gas pressure
                    VarNames = [VarNames, "P_1"];
                    VarValues = [VarValues, VarValues(:, 1)]; 
                    Initializer = initializer(VarNames, VarValues); % simple versoin
            
                otherwise
                    error(['DARSim Error: The fluid model keyword given in' ...
                        ' the input file is not valid. Please choose from ' ...
                        'the following options: SinglePhase,Immiscible,' ...
                        'Geothermal_SinglePhase,Geothermal_MultiPhase,' ...
                        'Compositional']);
            end

        end
        function Discretization = BuildDiscretization(obj)
            %% 1. Create fine-scale grids
            % 1. Reservoir Grid
            flag_unstructured_fracture = false;
            switch obj.SimulationInput.ReservoirProperties.Discretization
                case('CartesianGrid')
                    ReservoirGrid = cartesian_grid(...
                        obj.SimulationInput.ReservoirProperties.Grid.N);
                case('CornerPointGrid')
                    ReservoirGrid = corner_point_grid(...
                               obj.SimulationInput.ReservoirProperties);
                case('TetrahedralGrid')
                    ReservoirGrid = tetrahedral_grid(...
                                  obj.SimulationInput.ReservoirProperties);
                case('UnstructuredGrid')
                    ReservoirGrid = unstructured_grid(...
                     obj.SimulationInput.ReservoirProperties.UnstructuredGridData);
                    flag_unstructured_fracture = true;
                otherwise
                    error('DARSim Error: %s not yet implemented\n', ...
                        obj.SimulationInput.ReservoirProperties.Discretization);
            end
            % warning("Fracture mesh is finally ready here.")
            if obj.SimulationInput.FracturesProperties.isFractured && ...
                    ~flag_unstructured_fracture
                [FracturesGrid, CrossConnections] = ...
                    legacy_scanfracture_data(obj, ReservoirGrid);
                obj.SimulationInput.FracturesProperties.FracMesh_size = ...
                    FracturesGrid.N;
            elseif obj.SimulationInput.FracturesProperties.isFractured ...
                    && flag_unstructured_fracture
                [FracturesGrid, CrossConnections] = ...
                    scanfracture_data_uns(obj, ReservoirGrid);
                obj.SimulationInput.FracturesProperties.FracMesh_size = ...
                    FracturesGrid.N;
            end
            %warning("very important here")
            %% 2. Define your discretization Model (choose between FS and ADM)
            switch (obj.SimulatorSettings.DiscretizationModel)
                case('ADM')
                    % ADM discretization model (dynamic gridding)
                    ADMSettings = obj.SimulatorSettings.ADMSettings;
                    % Create the operatorshandler
                    operatorshandler = operators_handler_adm(ADMSettings.Coarsening(1,:,:));
                    
                    % a.1 Pressure prolongation builder
                    switch (ADMSettings.PInterpolator{1,1})
                        case ('Constant')
                            prolongationbuilder = prolongation_builder_constant(ADMSettings.maxLevel(1));
                        otherwise
                            prolongationbuilder = prolongation_builder_MSPressure(ADMSettings.maxLevel(1), ADMSettings.Coarsening(:,:,1));
                            if ~obj.SimulationInput.FracturesProperties.isFractured
                                switch obj.SimulatorSettings.Formulation
                                    case {'Geothermal_SinglePhase','Geothermal_MultiPhase','Molar_Energy'}
                                        prolongationbuilder.BFUpdater = bf_updater_ms_geothermal();
                                    otherwise
                                        prolongationbuilder.BFUpdater = bf_updater_ms();
                                end
                            else % Fractured
                                switch obj.SimulatorSettings.Formulation
                                    case {'Geothermal_SinglePhase','Geothermal_MultiPhase','Molar_Energy'}
                                        prolongationbuilder.BFUpdater = bf_updater_FAMS_geothermal();
                                    otherwise
                                        prolongationbuilder.BFUpdater = bf_updater_FAMS();
                                end
                                prolongationbuilder.BFUpdater.BFtype = ADMSettings.BFtype;
                            end
                            if strcmp(ADMSettings.PInterpolator, 'Homogeneous')
                                prolongationbuilder.BFUpdater.MaxContrast = 1;
                            else
                               prolongationbuilder.BFUpdater.MaxContrast = ADMSettings.BF_MaxContrast;
                            end
                            prolongationbuilder.BFUpdater.pEDFM_MaxContrast = ADMSettings.pEDFM_MaxContrast;
                            prolongationbuilder.alpha = ADMSettings.BF_alpha;
                    end
                    operatorshandler.AddProlongationBuilder(prolongationbuilder, 1);
                    
                    % a.2 Hyperbolic variables operators builder
                    n_comps = obj.SimulationInput.FluidProperties.NofComponents;
                    for i = 2:n_comps
                        switch(ADMSettings.HInterpolator)
                            case('Constant')
                                prolongationbuilder = prolongation_builder_constant(ADMSettings.maxLevel(1));
                            case('MS')
                                prolongationbuilder = prolongation_builder_MSHyperbolic(ADMSettings.maxLevel(1));
                                prolongationbuilder.key = strcat('S_',num2str(i-1));
                        end
                        operatorshandler.AddProlongationBuilder(prolongationbuilder, i);
                    end
                    if obj.SimulatorSettings.Formulation == "Geothermal_SinglePhase" || obj.SimulatorSettings.Formulation == "Geothermal_MultiPhase" || obj.SimulatorSettings.Formulation == "Molar_Energy"
                        switch(ADMSettings.HInterpolator)
                            case('Constant')
                                prolongationbuilder = prolongation_builder_constant(ADMSettings.maxLevel(1));
                            case('MS')
                                prolongationbuilder = prolongation_builder_MSHyperbolic(ADMSettings.maxLevel(1));
                                prolongationbuilder.key = 'T';
                        end
                        i = length(operatorshandler.ProlongationBuilders);
                        operatorshandler.AddProlongationBuilder(prolongationbuilder, i+1);
                    end

                    % b. Grid selection criterion (time\space based)
                    switch (ADMSettings.GridSelCriterion)
                        case('dfdx')
                            gridselector = adm_grid_selector_delta(ADMSettings.tol, ADMSettings.key);
                        case('dfdx2')
                            gridselector = adm_grid_selector_delta2(ADMSettings.tol, obj.SimulatorSettings.LTStol, ADMSettings.key);
                        case('dfdt')
                            gridselector = adm_grid_selector_time(ADMSettings.tol, ADMSettings.key, obj.SimulationInput.ReservoirProperties.Grid.N_ActiveCells, ADMSettings.maxLevel(1));
                        case('residual')
                            gridselector = adm_grid_selector_residual(ADMSettings.tol);
                    end

                    if strcmp(ADMSettings.Coupling,'COUPLED')
                        gridselector.isCoupled = 1;
                    else
                        gridselector.isCoupled = 0;
                    end
                    Discretization = ADM_Discretization_model(ADMSettings.maxLevel, ADMSettings.Coarsening);
                    Discretization.CoarseningSwitch = ADMSettings.CoarseningSwitch;
                    Discretization.AddADMGridSelector(gridselector);
                    Discretization.AddOperatorsHandler(operatorshandler);
                case ('MMs')
                    % Multilevel multiscale discretization model (static multiscale)
                    MMsSettings = obj.SimulatorSettings.MMsSettings;
                    % Create the operatorshandler
                    operatorshandler = operators_handler_MMs(MMsSettings.Coarsening(1,:,:));
                    switch (MMsSettings.PInterpolator)
                        case ('Constant')
                            prolongationbuilder = prolongation_builder_constant(MMsSettings.maxLevel(1));
                        otherwise
                            prolongationbuilder = prolongation_builder_MSPressure(MMsSettings.maxLevel(1), MMsSettings.Coarsening(:,:,1) );
                            if ~obj.SimulationInput.FracturesProperties.isFractured
                                switch obj.SimulatorSettings.Formulation
                                    case {'Geothermal_SinglePhase','Geothermal_MultiPhase'}
                                        prolongationbuilder.BFUpdater = bf_updater_ms_geothermal();
                                    otherwise
                                        prolongationbuilder.BFUpdater = bf_updater_ms();
                                end
                            else
                                switch obj.SimulatorSettings.Formulation
                                    case {'Geothermal_SinglePhase','Geothermal_MultiPhase'}
                                        prolongationbuilder.BFUpdater = bf_updater_FAMS_geothermal();
                                    otherwise
                                        prolongationbuilder.BFUpdater = bf_updater_FAMS();
                                end
                                prolongationbuilder.BFUpdater.BFtype = MMsSettings.BFtype;
                            end
                            % Reduce contrast for BF computation to remove peaks
                            if strcmp(MMsSettings.PInterpolatorBuildProductionSystem, 'Homogeneous')
                                prolongationbuilder.BFUpdater.MaxContrast = 1;
                            else
                                prolongationbuilder.BFUpdater.MaxContrast = MMsSettings.BF_MaxContrast;
                            end
                            prolongationbuilder.BFUpdater.pEDFM_MaxContrast = MMsSettings.pEDFM_MaxContrast;
                            prolongationbuilder.alpha = MMsSettings.BF_alpha;
                            
                            if MMsSettings.CorrectionFunctions
                                prolongationbuilder.BFUpdater.CorrectionFunctions = true;
                            end
                            
                    end
                    
                    % Add prolongation builder
                    operatorshandler.AddProlongationBuilder(prolongationbuilder, 1);
                    % Static Multiscale for flow solver
                    Discretization = Multiscale_Discretization_model(MMsSettings.maxLevel, MMsSettings.Coarsening);
                    Discretization.CoarseningSwitch = MMsSettings.CoarseningSwitch;
                    Discretization.AddOperatorsHandler(operatorshandler);
                case('FS')
                    % Fine-scale discretization model
                    switch obj.SimulationInput.ReservoirProperties.Discretization
                        case('CartesianGrid')
                            Discretization = Cartesian_Discretization_model();
                        case('CornerPointGrid')
                            Discretization = CornerPointGrid_Discretization_model();
                            Discretization.CornerPointGridData = obj.SimulationInput.ReservoirProperties.CornerPointGridData;
                        case('TetrahedralGrid')
                            Discretization = TetrahedralGrid_Discretization_model();
                            Discretization.TetrahedralGridData = obj.SimulationInput.ReservoirProperties.TetrahedralGridData;
                        case('UnstructuredGrid')
                            Discretization = UnstructuredGrid_Discretization_model();
                            Discretization.UnstructuredGridData = obj.SimulationInput.ReservoirProperties.UnstructuredGridData;                            
                        otherwise
                            error(['DARSim Error: At this moment,' ...
                                ' only "CartesianGrid", "CornerPointGrid",' ...
                                ' "UnstructuredGrid" discretization models' ...
                                ' are supported in DARSim!\n']);
                    end
            end
            
            %% 3. Add Grids to the Discretization Model
            Discretization.AddReservoirGrid(ReservoirGrid);
            Discretization.N = Discretization.ReservoirGrid.N;
            Discretization.type = ...
                obj.SimulationInput.ReservoirProperties.Discretization;
            if obj.SimulationInput.FracturesProperties.isFractured
                Discretization.AddFracturesGrid(FracturesGrid);
                Discretization.AddCrossConnections(CrossConnections);
                Discretization.N = Discretization.N + sum(Discretization.FracturesGrid.N);
            end
        end
        function ProductionSystem = BuildProductionSystemReservoir(obj)
            ProductionSystem = Production_System();
            %% RESERVOIR
            Lx = obj.SimulationInput.ReservoirProperties.size(1);       %Dimension in x-direction [m]
            Ly = obj.SimulationInput.ReservoirProperties.size(2);       %Dimension in y-direction [m]
            Lz  = obj.SimulationInput.ReservoirProperties.size(3);      %Reservoir thickness (z-direction) [m]
            Tres = obj.SimulationInput.ReservoirProperties.Temperature; %Res temperature [K]
            Reservoir = reservoir(Lx, Ly, Lz, Tres);
            RockCompressibility = obj.SimulationInput.RockProperties.Compressibility;
            RockSpecificHeat = obj.SimulationInput.RockProperties.SpecificHeat;
            RockDensity = obj.SimulationInput.RockProperties.Density;            
            % Permeability Data                      
            %%%% Adding permeability and Porosity info to the reservoir
            Reservoir.AddPermeabilityPorosity(obj.SimulationInput.ReservoirProperties.Perm, ...
                obj.SimulationInput.ReservoirProperties.phi, obj.SimulationInput.ReservoirProperties.KRotAngle);
            K_Coarse = obj.BuildProductionSystemCoarsePermeability();
                % Save them in ProductionSystem.
            Reservoir.AddCoarsePermeability(K_Coarse);
            % this function you have to create it            
            % Adding rock properties to the reservoir
            Reservoir.Cr = RockCompressibility;
            Reservoir.Cpr = RockSpecificHeat;
            Reservoir.Rho = RockDensity;
            Reservoir.P0 = obj.SimulationInput.InitialConditions.Pressure; 
            Reservoir.Diff = obj.SimulationInput.RockProperties.Diffusion;
            Reservoir.Disp = obj.SimulationInput.RockProperties.Dispersivity;
            if isfield(obj.SimulationInput.ReservoirProperties, "FaciesIndex")
                Reservoir.FaciesIndex = obj.SimulationInput.ReservoirProperties.FaciesIndex;
            else
                Reservoir.FaciesIndex = ones(length(Reservoir.Por), 1);
                obj.SimulationInput.ReservoirProperties.FaciesIndex = Reservoir.FaciesIndex;
            end            % Adding thermal conductivity to the reservoir (maybe find other place for this)
            if any(ismember(["Geothermal_SinglePhase", ...
                    "Geothermal_MultiPhase","Molar_Energy"], ...
                    obj.SimulatorSettings.Formulation))
                Reservoir.K_Cond_rock = ...
                    obj.SimulationInput.RockProperties.HeatConductivity;
            end
            %Add Reservoir to the ProductionSystem
            ProductionSystem.AddReservoir(Reservoir);
            Wells = obj.BuildProductionSystemWells();
            %ProductionSystem.Wells = Wells;
            ProductionSystem.AddWells(Wells);
            %ProductionSystem.AddReservoir(Wells);
            FracturesNetwork = ...
              obj.BuildProductionSystemFracture(Reservoir);
            ProductionSystem.AddFractures(FracturesNetwork);

            end
            function K_Coarse = BuildProductionSystemCoarsePermeability(obj)
                % This is for DLGR type ADM or MMs: it reads coarse permeabilities obtained by homogenization
                K_Coarse = [];
                switch obj.SimulatorSettings.DiscretizationModel
                    case{'ADM'}
                        if obj.SimulatorSettings.ADMSettings.DLGR
                            Add_Coarse_Permeability = 1;
                            MaxCoarseningLevel = obj.SimulatorSettings.ADMSettings.maxLevel(1);
                            CoarseningRatio = obj.SimulatorSettings.ADMSettings.Coarsening(1,:,1);
                        else
                            Add_Coarse_Permeability = 0;
                        end
                    case{'MMs'}
                        if obj.SimulatorSettings.MMsSettings.DLGR
                            Add_Coarse_Permeability = 1;
                            MaxCoarseningLevel = obj.SimulatorSettings.MMsSettings.maxLevel(1);
                            CoarseningRatio = obj.SimulatorSettings.MMsSettings.Coarsening(1,:,1);
                        else
                            Add_Coarse_Permeability = 0;
                        end
                    otherwise
                        Add_Coarse_Permeability = 0;
                end
                if Add_Coarse_Permeability
                    N_Finescale = obj.SimulationInput.ReservoirProperties.Grid.N;
                    Perm_CoarseScale_Homogenized = ...
                        DARSim2PermHomogenizer(N_Finescale, MaxCoarseningLevel, CoarseningRatio, K);
                    K_Coarse = cell(MaxCoarseningLevel + 1, 1);
                    K_Coarse{1} = K;
                    for L = 2 : MaxCoarseningLevel + 1
                        K_Coarse{L} = repmat(Perm_CoarseScale_Homogenized{L-1},1,3);
                    end
                end
            end
            function FracturesNetwork = BuildProductionSystemFracture(obj, Reservoir)
                %% FRACTURE NETWORK
                FracturesNetwork  = struct("Active",false,"NumOfFrac",0,...
                                                        "Fractures", []);            
                 
                if obj.SimulationInput.FracturesProperties.isFractured && ...
                        obj.SimulationInput.FracturesProperties.legacy
                 % STRUCTURED FRACTURE DARSIM 2.0
                    FracturesNetwork = ...
                                obj.BuildFractureNetworkLegacy(Reservoir);
                elseif obj.SimulationInput.FracturesProperties.isFractured ...
                    &&  ~obj.SimulationInput.FracturesProperties.legacy
                 % UNSTRUCTURED FRACTURE DARSIM 2.5
                    FracturesNetwork = ...
                            obj.BuildFractureNetworkLegacyUNS(Reservoir);
                end
            end


            function FracturesNetwork = BuildFractureNetworkLegacyUNS(obj, ...
                                        Reservoir)                      
                    num_frac = ...
                        obj.SimulationInput.FracturesProperties.NumOfFrac;
                    FracturesNetwork = fracture_system();
                    FracturesNetwork.Active = 1;
                    FracturesNetwork.NumOfFrac = num_frac;
                    FracturesNetwork.Fractures = fracture();
                    for f = 1 : num_frac
                        %DiscretizationModel.FracturesGrid.dfrac.geometry(f);
                        FracturesNetwork.Fractures(f).Length = ...
                            obj.SimulationInput.FracturesProperties.blueprint.FracturesProperties(f).Length;
                        FracturesNetwork.Fractures(f).Width = ...
                            obj.SimulationInput.FracturesProperties.blueprint.FracturesProperties(f).Width;
                        FracturesNetwork.Fractures(f).Thickness = ...
                            obj.SimulationInput.FracturesProperties.blueprint.FracturesProperties(f).Aperture;
                        FracturesNetwork.Fractures(f).Temp = obj.SimulationInput.InitialConditions.Temperature;     % Temperature of each fracture
                        FracturesNetwork.Fractures(f).Cpr = 0;
                        FracturesNetwork.Fractures(f).Cr  = 0;
                        FracturesNetwork.Fractures(f).Rho = 0;
                        
                        if obj.SimulationInput.InitialConditions.Include
                            FracturesNetwork.Fractures(f).P0 = obj.SimulationInput.InitialConditions.LoadedData(:,1);
                        else
                            FracturesNetwork.Fractures(f).P0 = obj.SimulationInput.InitialConditions.Pressure;
                        end
                        Porosity = obj.SimulationInput.FracturesProperties.blueprint.FracturesProperties(f).Porosity; % Porosity  of each fracture
                        Permeability = obj.SimulationInput.FracturesProperties.blueprint.FracturesProperties(f).Permeability;     % Permeability of each fracture
                        if obj.SimulationInput.FracturesProperties.blueprint.FracturesProperties(f).isActive
                            Kx = ones(obj.SimulationInput.FracturesProperties.FracMesh_size(f), 1)*Permeability;
                            Ky = Kx;
                            Kz = Kx;
                            K = [Kx, Ky, Kz];
                            FracturesNetwork.Fractures(f).AddPermeabilityPorosity(K, Porosity);                 % Adding porosity and permeability to the fracture
                            switch obj.SimulatorSettings.Formulation
                                case {'Geothermal_SinglePhase','Geothermal_MultiPhase','Molar_Energy'}
                                    FracturesNetwork.Fractures(f).K_Cond_rock = obj.SimulationInput.RockProperties.HeatConductivity;
                            end
                        end
                        FracturesNetwork.Fractures(f).Diff = zeros(length(K),2);
                        FracturesNetwork.Fractures(f).Disp = zeros(length(K),1);
                    end
            end


            function FracturesNetwork = BuildFractureNetworkLegacy(obj, ...
                                        Reservoir)
                    FractureMatrix = obj.SimulationInput. ...
                        FracturesProperties.FractureMatrix;    
                    FracturesNetwork = fracture_system();
                    FracturesNetwork.Active = 1;
                    FracMesh_size = ...
                     obj.SimulationInput.FracturesProperties.FracMesh_size;
                    temp = strfind(FractureMatrix, 'NUM_FRACS');
                    index = find(~cellfun('isempty', temp));
                    temp = strsplit(FractureMatrix{index},' ');
                    NumOfFrac = temp{end};
                    FracturesNetwork.NumOfFrac = str2double( NumOfFrac );
                    temp = strfind(FractureMatrix, 'PROPERTIES');
                    frac_index = find(~cellfun('isempty', temp));
                    FracturesNetwork.Fractures = fracture();
                    for f = 1 : FracturesNetwork.NumOfFrac
                        frac_info_split = strsplit(FractureMatrix{frac_index(f)},' ');                      % Splitted data for each fracture
                        FracturesNetwork.Fractures(f).Length = str2double( frac_info_split{3} );            % Length of each fracture
                        FracturesNetwork.Fractures(f).Width = str2double( frac_info_split{4} );             % Width of each fracture
                        FracturesNetwork.Fractures(f).Thickness = str2double( frac_info_split{5} );         % Thickness of each fracture
                        FracturesNetwork.Fractures(f).Temp = Reservoir.Temp; % Temperature of each fracture
                        FracturesNetwork.Fractures(f).Cpr = Reservoir.Cpr;
                        FracturesNetwork.Fractures(f).Cr  = Reservoir.Cr;
                        FracturesNetwork.Fractures(f).Rho = Reservoir.Rho;
                        if obj.SimulationInput.InitialConditions.Include
                            FracturesNetwork.Fractures(f).P0 = obj.SimulationInput.InitialConditions.LoadedData(:,1);
                        else
                            FracturesNetwork.Fractures(f).P0 = obj.SimulationInput.InitialConditions.Pressure;
                        end
                        Porosity = str2double( frac_info_split{6} );                                        % Porosity  of each fracture
                        Permeability = str2double( frac_info_split{7} );                                    % Permeability of each fracture
                        Kx = ones(FracMesh_size(f), 1)*Permeability;
                        Ky = Kx;
                        Kz = Kx;
                        K = [Kx, Ky, Kz];
                        FracturesNetwork.Fractures(f).AddPermeabilityPorosity(K, Porosity);                 % Adding porosity and permeability to the fracture
                        switch obj.SimulatorSettings.Formulation
                            case {'Geothermal_SinglePhase','Geothermal_MultiPhase','Molar_Energy'}
                                FracturesNetwork.Fractures(f).K_Cond_rock = obj.SimulationInput.RockProperties.HeatConductivity;
                        end

                        FracturesNetwork.Fractures(f).Diff = zeros(length(K),2);
                        FracturesNetwork.Fractures(f).Disp = zeros(length(K),1);
                    end
            end

        function Wells = BuildProductionSystemWells(obj)
            %% WELLS
            Lx = obj.SimulationInput.ReservoirProperties.size(1);     
            Ly = obj.SimulationInput.ReservoirProperties.size(2);      
            Lz  = obj.SimulationInput.ReservoirProperties.size(3); 
            nx = obj.SimulationInput.ReservoirProperties.Grid.N(1);
            ny = obj.SimulationInput.ReservoirProperties.Grid.N(2);
            nz = obj.SimulationInput.ReservoirProperties.Grid.N(3);
            N_ActiveCells = ...
                obj.SimulationInput.ReservoirProperties.Grid.N_ActiveCells;                        
            Wells = wells();
            Wells.NofInj = obj.SimulationInput.WellsInfo.NofInj;
            Wells.NofProd = obj.SimulationInput.WellsInfo.NofProd;
            n_phases = obj.SimulationInput.FluidProperties.NofPhases;
            switch obj.SimulationInput.ReservoirProperties.Discretization
                case('CartesianGrid')
                    % Useful information to compute PI
                    dx = Lx/nx; dy = Ly/ny; dz = Lz /nz;
                    GridVolume = dx*dy*dz * ones(N_ActiveCells,1);
                case('CornerPointGrid')
                    % Useful information to compute PI
                    dx = Lx/nx; dy = Ly/ny; dz = Lz /nz;
                    GridVolume = obj.SimulationInput.ReservoirProperties. ...
                        CornerPointGridData.Cells.Volume;
                case('TetrahedralGrid')
                    GridVolume = obj.SimulationInput.ReservoirProperties. ...
                        TetrahedralGridData.Cell.Volume;
                case('UnstructuredGrid')
                    GridVolume = obj.SimulationInput.ReservoirProperties. ...
                        UnstructuredGridData.element.volume;
                otherwise
                    Error("Well not defined for this " + ...
                        "                   discretization model.");
            end
            % Injectors
            for i=1:Wells.NofInj
                switch(obj.SimulationInput.WellsInfo.Inj{i}.Formula.Type)
                    case('DIRICHLET')
                        PI = dy*dz/(dx/2) .* ones(N_ActiveCells,1);
                    case('PI')
                        PI = obj.SimulationInput.WellsInfo.Inj{i}.Formula.Value .* ones(N_ActiveCells,1);
                    case('WI')
                        PI = obj.SimulationInput.WellsInfo.Inj{i}.Formula.Value .* GridVolume;
                    case('RADIUS')
                        %radius = obj.SimulationInput.WellsInfo.Inj{i}.Formula.Value;
                        error('DARSim Error: Radius calculation of PI is not implemented for now.')
                end
                Enthalpy = NaN;
                Tres = NaN;
                switch (obj.SimulationInput.FluidProperties.FluidModel)
                    case {'Geothermal_SinglePhase','Geothermal_MultiPhase','Compositional_Energy'}
                        switch obj.SimulationInput.WellsInfo.Inj{i}.BoundaryCondition.Name
                            case('TEMPERATURE')
                                Temperature = obj.SimulationInput.WellsInfo.Inj{i}.BoundaryCondition.Value;
                                BC_Formulation = 'Temperature';
                                Enthalpy = NaN;
                            case('ENTHALPY')
                                Temperature = NaN;
                                BC_Formulation = 'Enthalpy';
                                Enthalpy = obj.SimulationInput.WellsInfo.Inj{i}.BoundaryCondition.Value;
                            otherwise
                                error('DARSim Error: The boundary conditions for the injection wells are not defined properly in the input file.')
                        end
                    otherwise
                        Temperature = Tres;
                end
                Coordinate = obj.SimulationInput.WellsInfo.Inj{i}.Coordinate;
                switch (obj.SimulationInput.WellsInfo.Inj{i}.Constraint.Name)
                    case('PRESSURE')
                        pressure = obj.SimulationInput.WellsInfo.Inj{i}.Constraint.Value;
                        Injector = injector_pressure(PI, Coordinate, pressure, Temperature, n_phases);
                        switch (obj.SimulationInput.FluidProperties.FluidModel)
                            case {'Geothermal_SinglePhase','Geothermal_MultiPhase'}
                                Injector.BC_Formulation = BC_Formulation;
                        end
                    case('RATE')
                        rate = obj.SimulationInput.WellsInfo.Inj{i}.Constraint.Value;
                        p_init = obj.SimulationInput.InitialConditions.Pressure;
                        %rate = rate * TotalPV(1) / (3600 * 24); % convert pv/day to m^3/s
                        rate = rate * sum(GridVolume .* obj.SimulationInput.ReservoirProperties.phi) / (3600 * 24); % convert pv/day to m^3/s
                        Injector = injector_rate(PI, Coordinate, rate, p_init, Temperature, n_phases);
                    case('MASSRATE')
                        rate = obj.SimulationInput.WellsInfo.Inj{i}.Constraint.Value;
                        p_init = obj.SimulationInput.InitialConditions.Pressure;
                        %rate = rate * TotalPV(1) / (3600 * 24); % convert pv/day to m^3/s
                        rate = rate * 1000 / (3600 * 24); % convert ton/day to kg/s
                        Injector = injector_mass_rate(PI, Coordinate, rate, p_init, Temperature, n_phases);
                end
                Injector.h = Enthalpy;
                % Check to see if the well has any schedules
                try 
                    Injector.Schedule = obj.SimulationInput.WellsInfo.Inj{i}.Schedule * (3600 * 24);
                catch
                    % Always open well
                    Injector.Schedule = [0 Inf];
                end
                % Check to see if the well has any specific composition
                try 
                    Injector.Composition = obj.SimulationInput.WellsInfo.Inj{i}.Composition;
                catch
                    % Default injection is first composition
                    warning('Injector %d will be injecting component 1.', i);
                    Injector.Composition = 1;
                end
                Wells.AddInjector(Injector);
            end
            % Producers
            for i=1:Wells.NofProd
                switch(obj.SimulationInput.WellsInfo.Prod{i}.Formula.Type)
                    case ('DIRICHLET')
                        PI = dy*dz/(dx/2) .* ones(N_ActiveCells,1);
                    case ('PI')
                        PI = obj.SimulationInput.WellsInfo.Prod{i}.Formula.Value .* ones(N_ActiveCells,1);
                    case('WI')
                        PI = obj.SimulationInput.WellsInfo.Prod{i}.Formula.Value .* GridVolume;
                    case ('RADIUS')
                        %radius = obj.SimulationInput.WellsInfo.Prod{i}.Formula.Value;
                        error('DARSim2 Error: Radius calculation of PI is not implemented for now.')
                end
                Coordinate = obj.SimulationInput.WellsInfo.Prod{i}.Coordinate;
                switch (obj.SimulationInput.WellsInfo.Prod{i}.Constraint.Name)
                    case('PRESSURE')
                        pressure = obj.SimulationInput.WellsInfo.Prod{i}.Constraint.Value;
                        Producer = producer_pressure(PI, Coordinate, pressure);
                    case('RATE')
                        rate = obj.SimulationInput.WellsInfo.Prod{i}.Constraint.Value;
                        rate = rate * sum(GridVolume .* obj.SimulationInput.ReservoirProperties.phi) / (3600 * 24); % convert pv/day to m^3/s
                        Producer = producer_rate(PI, Coordinate, rate);
                    case('MASSRATE')
                        rate = obj.SimulationInput.WellsInfo.Inj{i}.Constraint.Value;
                        p_init = obj.SimulationInput.InitialConditions.Pressure;
                        %rate = rate * TotalPV(1) / (3600 * 24); % convert pv/day to m^3/s
                        rate = rate * 1000 / (3600 * 24); % convert ton/day to kg/s
                        Producer = producer_mass_rate(PI, Coordinate, rate, p_init, Temperature, n_phases);
                end
                % Check to see if the well has any schedules
                try 
                    Producer.Schedule = obj.SimulationInput.WellsInfo.Prod{i}.Schedule * (3600 * 24);
                catch
                    % Always open well
                    Producer.Schedule = [0 Inf];
                end
                Wells.AddProducer(Producer);
            end
        end
        function FluidModel = BuildFluidModel(obj)
            n_phases = obj.SimulationInput.FluidProperties.NofPhases;
            switch(obj.SimulationInput.FluidProperties.FluidModel)
                case('SinglePhase')
                    FluidModel = single_phase_fluid_model();
                    % Add phase
                    Phase = comp_phase();
                    % Gets all densities [kg/m^3]
                    Phase.rho0 = obj.SimulationInput.FluidProperties.Density;
                    % Gets the intial pressure as the reference
                    Phase.P0 = obj.SimulationInput.InitialConditions.Pressure;
                    % Gets all viscosities [Pa sec]
                    Phase.mu = obj.SimulationInput.FluidProperties.mu;
                    % Compressibility
                    Phase.cf = obj.SimulationInput.FluidProperties.Comp;
                    FluidModel.AddPhase(Phase, 1);
                    if Phase.cf == 0
                        obj.incompressible = 1;
                    end
                    obj.SimulatorSettings.CouplingType = 'SinglePhase';
                case('Immiscible')
                    FluidModel = Immiscible_fluid_model(n_phases);
                    % Add phases
                    for i = 1:FluidModel.NofPhases
                        Phase = comp_phase();
                        %Gets all densities [kg/m^3]
                        Phase.rho0 = obj.SimulationInput.FluidProperties.Density(i);
                        %Gets all viscosities [Pa sec]
                        Phase.mu = obj.SimulationInput.FluidProperties.mu(i);
                        % Compressibility
                        Phase.cf = obj.SimulationInput.FluidProperties.Comp(i);
                        FluidModel.AddPhase(Phase, i);
                    end
                    obj.SimulatorSettings.Formulation = 'Immiscible';
                case('BlackOil')
                    n_comp = 2; % Number of components for BlackOil model should be 2
                    FluidModel = BO_fluid_model(n_phases, n_comp);
                    switch(obj.SimulationInput.FluidProperties.ComponentList)
                        case('CO2-Brine')
                            % Black oil phases
                            % Gas
                            GasPhase = BO_gas_phase(obj.SimulationInput.FluidProperties.mu(1), obj.SimulationInput.InitialConditions.Temperature);
                            FluidModel.AddPhase(GasPhase, 1);
                            GasComponent = CO2_component(obj.SimulationInput.InitialConditions.Temperature);
                            FluidModel.AddComponent(GasComponent, 1);
                            % Oil
                            OilPhase = BO_oil_phase(obj.SimulationInput.FluidProperties.mu(2));
                            FluidModel.AddPhase(OilPhase, 2);
                            OilComponent = H2O_component();
                            FluidModel.AddComponent(OilComponent, 2);
                            Thermo = Thermo_CO2_Brine_BlackOil(obj.SimulationInput.InitialConditions.Temperature, obj.SimulationInput.FluidProperties.Salinity);
                            FluidModel.Thermo = Thermo;
                            FluidModel.Phases(2).Rs_p_table   = Thermo.Rs_p_table;
                            FluidModel.Phases(2).dRs_p_table  = Thermo.dRs_p_table;
                            FluidModel.Phases(1).zcf_p_table  = Thermo.zcf_p_table;
                            FluidModel.Phases(1).dzcf_p_table = Thermo.dzcf_p_table;
                        case('H2-Brine')
                            % Black oil phases
                            % Gas
                            GasPhase = BO_gas_phase(obj.SimulationInput.FluidProperties.mu(1), obj.SimulationInput.InitialConditions.Temperature);
                            FluidModel.AddPhase(GasPhase, 1);
                            GasComponent = H2_component();
                            FluidModel.AddComponent(GasComponent, 1);
                            % Oil
                            OilPhase = BO_oil_phase(obj.SimulationInput.FluidProperties.mu(2));
                            FluidModel.AddPhase(OilPhase, 2);
                            OilComponent = H2O_component();
                            FluidModel.AddComponent(OilComponent, 2);
                            Thermo = Thermo_H2_Brine_BlackOil(obj.SimulationInput.InitialConditions.Temperature, obj.SimulationInput.FluidProperties.Salinity);
                            FluidModel.Thermo = Thermo;
                            FluidModel.Phases(2).Rs_p_table   = Thermo.Rs_p_table;
                            FluidModel.Phases(2).dRs_p_table  = Thermo.dRs_p_table;
                            FluidModel.Phases(1).zcf_p_table  = Thermo.zcf_p_table;
                            FluidModel.Phases(1).dzcf_p_table = Thermo.dzcf_p_table;
                        otherwise
                            error('The components are not interpreted for BlackOil model.')
                            % if FluidModel.NofPhases == 3
                            %     Water = comp_phase();
                            %     Water.rho0 = 1000; % kg/m^3
                            %     Water.mu = 1e-3; % Pa s
                            %     Water.cf = 0;
                            %     FluidModel.AddPhase(Water, 3)
                            %     water = component();
                            %     FluidModel.AddComponent(water, 3);
                            % end
                    end
                            FluidModel.FlashCalculator = Rachford_Rice_flash_calculator();
                            % FlashCalculator = Standard_flash_calculator();
                            FluidModel.FlashCalculator.KvaluesCalculator = BO_Kvalues_calculator();
                case('Compositional')
                    n_comp = obj.SimulationInput.FluidProperties.NofComponents;
                    FluidModel = Compositional_fluid_model(n_phases, n_comp);
                    switch(obj.SimulationInput.FluidProperties.ComponentList)
                        case('CO2-Brine')
                            % Compositional phases
                            % Gas
                            GasPhase = Compositional_gas_phase(obj.SimulationInput.FluidProperties.mu(1), obj.SimulationInput.InitialConditions.Temperature);
                            FluidModel.AddPhase(GasPhase, 1);
                            GasComponent = CO2_component(obj.SimulationInput.InitialConditions.Temperature);
                            FluidModel.AddComponent(GasComponent, 1);
                            % Liquid
                            LiqPhase = Compositional_water_phase(obj.SimulationInput.FluidProperties.mu(2));
                            FluidModel.AddPhase(LiqPhase, 2);
                            LiqComponent = H2O_component();
                            FluidModel.AddComponent(LiqComponent, 2);
                            Thermo = Thermo_CO2_Brine_Compositional(obj.SimulationInput.InitialConditions.Temperature, obj.SimulationInput.FluidProperties.Salinity);
                            FluidModel.Thermo = Thermo;
                            FluidModel.Phases(1).zcf_p_table  = Thermo.zcf_p_table;
                            FluidModel.Phases(1).dzcf_p_table = Thermo.dzcf_p_table;
                            FlashCalculator = Rachford_Rice_flash_calculator_Compositional();
                            FlashCalculator.KvaluesCalculator = Compositional_Kvalues_calculator();
                            FluidModel.FlashCalculator = FlashCalculator;
                        case('H2-Brine')
                            % Compositional phases
                            % Gas
                            GasPhase = Compositional_gas_phase(obj.SimulationInput.FluidProperties.mu(1), obj.SimulationInput.InitialConditions.Temperature);
                            FluidModel.AddPhase(GasPhase, 1);
                            GasComponent = H2_component();
                            FluidModel.AddComponent(GasComponent, 1);
                            % Liquid
                            LiqPhase = Compositional_water_phase(obj.SimulationInput.FluidProperties.mu(2));
                            FluidModel.AddPhase(LiqPhase, 2);
                            LiqComponent = H2O_component();
                            FluidModel.AddComponent(LiqComponent, 2);
                            Thermo = Thermo_H2_Brine_Compositional(obj.SimulationInput.InitialConditions.Temperature, obj.SimulationInput.FluidProperties.Salinity);
                            FluidModel.Thermo = Thermo;
                            FluidModel.Phases(1).zcf_p_table  = Thermo.zcf_p_table;
                            FluidModel.Phases(1).dzcf_p_table = Thermo.dzcf_p_table;
                            FlashCalculator = Rachford_Rice_flash_calculator_Compositional();
                            FlashCalculator.KvaluesCalculator = Compositional_Kvalues_calculator();
                            FluidModel.FlashCalculator = FlashCalculator;
                        case('H2-CH4-Brine') %3Comp
                            % Compositional phases
                            % Gas
                            GasPhase = Compositional_gas_phase_H2_CH4(obj.SimulationInput.FluidProperties.mu(1), obj.SimulationInput.InitialConditions.Temperature);
                            FluidModel.AddPhase(GasPhase, 1);
                            GasComponent = H2_component();
                            FluidModel.AddComponent(GasComponent, 1);
                            GasComponent = CH4_component();
                            FluidModel.AddComponent(GasComponent, 2);
                            % Liquid
                            LiqPhase = Compositional_water_phase(obj.SimulationInput.FluidProperties.mu(2));
                            FluidModel.AddPhase(LiqPhase, 2);
                            LiqComponent = H2O_component();
                            FluidModel.AddComponent(LiqComponent, 3);
                            Thermo = Thermo_H2_CH4_Brine_Compositional(obj.SimulationInput.InitialConditions.Temperature, obj.SimulationInput.FluidProperties.Salinity);
                            FluidModel.Thermo = Thermo;
                            FluidModel.Phases(1).zcf  = Thermo.zcf;
                            FluidModel.Phases(1).dzcfdp = Thermo.dzcfdp;    
                            FluidModel.Phases(1).dzcfdy1 = Thermo.dzcfdy1;    
                            FlashCalculator = Rachford_Rice_3component();
                            FlashCalculator.KvaluesCalculator = Kvalues_calculator_3component();
                            FluidModel.FlashCalculator = FlashCalculator;
                    end
                    % FlashCalculator = Rachford_Rice_flash_calculator_Compositional();
                    % if n_comp > 2
                    %     FlashCalculator.KvaluesCalculator = Wilson_Kvalues_calculator(obj.SimulationInput.ReservoirProperties.Temperature);
                    % else
                    %     FlashCalculator.KvaluesCalculator = Compositional_Kvalues_calculator();
                    % end
                    % FluidModel.FlashCalculator = FlashCalculator;
                case('Compositional_Energy')
                    n_comp = obj.SimulationInput.FluidProperties.NofComponents;
                    FluidModel = Compositional_Energy_fluid_model(n_phases, n_comp);
                    switch(obj.SimulationInput.FluidProperties.ComponentList)
                        case('CO2-Brine')
                            % Compositional Energy phases
                            % Gas     
                            GasPhase = Compositional_Energy_gas_phase(obj.SimulationInput.FluidProperties.FluidConductivity(1));
                            FluidModel.AddPhase(GasPhase, 1);
                            GasComponent = CO2_component(obj.SimulationInput.InitialConditions.Temperature);
                            FluidModel.AddComponent(GasComponent, 1);
                            % Liquid
                            LiqPhase = Compositional_Energy_water_phase(obj.SimulationInput.FluidProperties.FluidConductivity(2));
                            FluidModel.AddPhase(LiqPhase, 2);
                            LiqComponent = H2O_component();
                            FluidModel.AddComponent(LiqComponent, 2);
                            Thermo = Thermo_CO2_Brine_Compositional_Energy(obj.SimulationInput.FluidProperties.Salinity);
                            FluidModel.Thermo = Thermo;
                            FlashCalculator = Rachford_Rice_flash_calculator_Compositional();
                            FlashCalculator.KvaluesCalculator = Compositional_Kvalues_calculator();
                            FluidModel.FlashCalculator = FlashCalculator;
                        case('H2-Brine')
                            % Compositional phases
                            % Gas
                            GasPhase = Compositional_gas_phase(obj.SimulationInput.FluidProperties.mu(1), obj.SimulationInput.InitialConditions.Temperature);
                            FluidModel.AddPhase(GasPhase, 1);
                            GasComponent = H2_component();
                            FluidModel.AddComponent(GasComponent, 1);
                            % Liquid
                            LiqPhase = Compositional_water_phase(obj.SimulationInput.FluidProperties.mu(2));
                            FluidModel.AddPhase(LiqPhase, 2);
                            LiqComponent = H2O_component();
                            FluidModel.AddComponent(LiqComponent, 2);
                            Thermo = Thermo_H2_Brine_Compositional(obj.SimulationInput.InitialConditions.Temperature, obj.SimulationInput.FluidProperties.Salinity);
                            FluidModel.Thermo = Thermo;
                            FluidModel.Phases(1).zcf_p_table  = Thermo.zcf_p_table;
                            FluidModel.Phases(1).dzcf_p_table = Thermo.dzcf_p_table;
                            FlashCalculator = Rachford_Rice_flash_calculator_Compositional();
                            FlashCalculator.KvaluesCalculator = Compositional_Kvalues_calculator();
                            FluidModel.FlashCalculator = FlashCalculator;
                        case('H2-CH4-Brine') %3Comp
                            % Compositional phases
                            % Gas
                            GasPhase = Compositional_gas_phase_H2_CH4(obj.SimulationInput.FluidProperties.mu(1), obj.SimulationInput.InitialConditions.Temperature);
                            FluidModel.AddPhase(GasPhase, 1);
                            GasComponent = H2_component();
                            FluidModel.AddComponent(GasComponent, 1);
                            GasComponent = CH4_component();
                            FluidModel.AddComponent(GasComponent, 2);
                            % Liquid
                            LiqPhase = Compositional_water_phase(obj.SimulationInput.FluidProperties.mu(2));
                            FluidModel.AddPhase(LiqPhase, 2);
                            LiqComponent = H2O_component();
                            FluidModel.AddComponent(LiqComponent, 3);
                            Thermo = Thermo_H2_CH4_Brine_Compositional(obj.SimulationInput.InitialConditions.Temperature, obj.SimulationInput.FluidProperties.Salinity);
                            FluidModel.Thermo = Thermo;
                            FluidModel.Phases(1).zcf  = Thermo.zcf;
                            FluidModel.Phases(1).dzcfdp = Thermo.dzcfdp;    
                            FluidModel.Phases(1).dzcfdy1 = Thermo.dzcfdy1;    
                            FlashCalculator = Rachford_Rice_3component();
                            FlashCalculator.KvaluesCalculator = Kvalues_calculator_3component();
                            FluidModel.FlashCalculator = FlashCalculator;
                    end

                case{'Geothermal_SinglePhase'}
                    % build the geothermal singlephase fluid model
                    FluidModel = Geothermal_SinglePhase_fluid_model();
                    Phase = Geothermal_SinglePhase_phase();
                    FluidModel.AddPhase(Phase, 1);
                    obj.SimulatorSettings.CouplingType = 'FIM';
                    %Gets all densities [kg/m^3]
                    Phase.rho0 = obj.SimulationInput.FluidProperties.Density;
                    %Gets all viscosities [Pa sec]
                    Phase.mu0 = obj.SimulationInput.FluidProperties.mu;
                    % Compressibility
                    Phase.cf0 = obj.SimulationInput.FluidProperties.Comp;
                    % Conductivity
                    Phase.Kf = obj.SimulationInput.FluidProperties.FluidConductivity;
                    % Specific Heat
                    Phase.Cp_std = obj.SimulationInput.FluidProperties.SpecificHeat;
                    FluidModel.AddPhase(Phase, 1);
                    if Phase.cf0 == 0
                        obj.incompressible = 1;
                    end
                case{'Geothermal_MultiPhase'}
                    % build the geothermal multiphase fluid model
                    FluidModel = Geothermal_Multiphase_fluid_model(n_phases);
                    % Add phases
                    for i = 1:FluidModel.NofPhases
                        Phase = Geothermal_Multiphase_phase();
                        % we are only using cp_std in the injection wells
                        Phase.Cp_std = obj.SimulationInput.FluidProperties.SpecificHeat(i); % Specific Heat
                        Phase.Kf = obj.SimulationInput.FluidProperties.FluidConductivity(i);   % Conductivity
                        FluidModel.AddPhase(Phase, i);
                    end
                    obj.SimulatorSettings.CouplingType = 'FIM';
                    
                otherwise
                    error('FluidModel is not defined!');
            end
            
            %% RelPerm model
            switch(obj.SimulationInput.FluidProperties.RelPerm.name)
                case('Linear')
                    FluidModel.RelPermModel = relperm_model_linear();
                case('Quadratic')
                    FluidModel.RelPermModel = relperm_model_quadratic();
                case('Foam')
                    FluidModel.RelPermModel = relperm_model_foam();
                case('Corey')
                    FluidModel.RelPermModel = relperm_model_brookscorey();
                case('Table')
                    FluidModel.RelPermModel = relperm_model_table(obj.SimulationInput);
                case('User-defined')
                    % Reading the original User-defined code
                    originalCode = fullfile(obj.SimulationInput.Directory, ...
                        obj.SimulationInput.FluidProperties.RelPerm.Files);
                    
                    txt = fileread(originalCode);
                    [~, originalClassName, ~] = fileparts(obj.SimulationInput.FluidProperties.RelPerm.Files);
                    
                    % Replace the classdef name to
                    % "relperm_model_user_defined"
                    txt = regexprep(txt, ...
                        "classdef\s+" + originalClassName, ...
                        "classdef " + "relperm_model_user_defined");
                    
                    % Copy the original file to \RelativePermeability
                    % folder with the name of "relperm_model_user_defined"
                    allPaths = split(path, pathsep);
                    matches = allPaths(endsWith(allPaths, filesep + "RelativePermeability") | allPaths == "RelativePermeability");
                    newFile = fullfile(matches{1}, "relperm_model_user_defined.m");
                    fid = fopen(newFile, "w");
                    fwrite(fid, txt);
                    fclose(fid);

                    FluidModel.RelPermModel = relperm_model_user_defined(obj.SimulationInput);

            end

            % Relperm of fractures
            if obj.SimulationInput.FracturesProperties.isFractured
                FluidModel.RelPermModelFractures = relperm_model_fractures(obj.SimulationInput.FluidProperties.RelPerm.srw_frac);
            end
            
            %% Capillary pressure model
            switch (obj.SimulationInput.FluidProperties.Capillarity.name)
                case('JLeverett')
                    obj.capillarymodel = true;
                    FluidModel.CapillaryModel = J_Function_model();
                    FluidModel.WettingPhaseIndex = ...
                        obj.SimulationInput.FluidProperties.Capillarity.wetting;
                case('Linear')
                    FluidModel.CapillaryModel = 'Not implemented';
                case('BrooksCorey')
                    obj.capillarymodel = true;
                    FluidModel.CapillaryModel = brookscorey_model();
                    FluidModel.WettingPhaseIndex = ...
                        obj.SimulationInput.FluidProperties.Capillarity.wetting;
                case('Table')
                    obj.capillarymodel = true;
                    FluidModel.CapillaryModel = ...
                                        capillary_pressure_model_table();
                case('User-defined')
                    obj.capillarymodel = true;

                    % Reading the original User-defined code
                    originalCode = fullfile(obj.SimulationInput.Directory, ...
                        obj.SimulationInput.FluidProperties.Capillarity.Files);
                    
                    txt = fileread(originalCode);
                    [~, originalClassName, ~] = fileparts(obj.SimulationInput.FluidProperties.Capillarity.Files);
                    
                    % Replace the classdef name to
                    % "capillary_pressure_model_user_defined"
                    txt = regexprep(txt, ...
                        "classdef\s+" + originalClassName, ...
                        "classdef " + "capillary_pressure_model_user_defined");
                    
                    % Copy the original file to \Capillarity
                    % folder with the name of "capillary_pressure_model_user_defined"
                    allPaths = split(path, pathsep);
                    matches = allPaths(endsWith(allPaths, filesep + "Capillarity") | allPaths == "Capillarity");
                    newFile = fullfile(matches{1}, "capillary_pressure_model_user_defined.m");
                    fid = fopen(newFile, "w");
                    fwrite(fid, txt);
                    fclose(fid);

                    FluidModel.CapillaryModel = capillary_pressure_model_user_defined();  

                otherwise
                    FluidModel.CapillaryModel = No_Pc_model();
            end
            
        end
        function Formulation = BuildFormulation(obj)
            switch(obj.SimulatorSettings.Formulation)
                case('Immiscible')
                    Formulation = Immiscible_formulation();
                    if obj.SimulatorSettings.LTS
                        Formulation.MatrixAssembler = LTS_matrix_assembler();
                    else
                        Formulation.MatrixAssembler = matrix_assembler();
                    end
                    obj.NofEq = obj.SimulationInput.FluidProperties.NofPhases;
                case('Natural')
                    Formulation = NaturalVar_formulation(sum(obj.SimulationInput.ReservoirProperties.Grid.N), obj.SimulationInput.FluidProperties.NofComponents);
                    obj.NofEq = obj.SimulationInput.FluidProperties.NofPhases + obj.SimulationInput.FluidProperties.NofComponents;
                case('Molar')
                    Formulation = Overall_Composition_formulation(obj.SimulationInput, obj.SimulationInput.FluidProperties.NofPhases, obj.SimulationInput.FluidProperties.NofComponents);
                    obj.NofEq = obj.SimulationInput.FluidProperties.NofComponents;
                case('Molar_Energy')
                    Formulation = Overall_Composition_Energy_formulation(obj.SimulationInput, obj.SimulationInput.FluidProperties.NofPhases, obj.SimulationInput.FluidProperties.NofComponents);
                    obj.NofEq = obj.SimulationInput.FluidProperties.NofComponents + 1; % one more energy balance
                case('OBL')
                    Formulation = OBL_formualtion();
                    Formulation.CreateTables();
                case('Geothermal_SinglePhase')
                    Formulation = Geothermal_SinglePhase_formulation();
                    Formulation.MatrixAssembler = matrix_assembler_geothermal();
                    obj.NofEq = 2;
                case('Geothermal_MultiPhase')
                    Formulation = Geothermal_MultiPhase_formulation(obj.SimulationInput.FluidProperties.NofPhases);
                    Formulation.MatrixAssembler = matrix_assembler_geothermal();
                    obj.NofEq = 2;
            end
            Formulation.NofPhases = obj.SimulationInput.FluidProperties.NofPhases;
        end
        function TimeDriver = BuildTimeDriver(obj)
            if obj.SimulatorSettings.LTSPlot == 1
                TimeDriver = LTSPrint_TimeLoop_Driver(obj.SimulatorSettings.Reports, obj.SimulationInput.TotalTime, obj.SimulatorSettings.MaxNumTimeSteps);
                % elseif obj.SimulatorSettings.ADT_SEQ == 1
                %   TimeDriver = LTS_TimeLoop_Driver(obj.SimulatorSettings.Reports, obj.SimulationInput.TotalTime, obj.SimulatorSettings.MaxNumTimeSteps);
            else
                TimeDriver = TimeLoop_Driver(obj.SimulatorSettings.Reports, obj.SimulationInput.TotalTime, obj.SimulatorSettings.MaxNumTimeSteps);
            end
            % Construct Coupling
            switch(obj.SimulatorSettings.CouplingType)
                case('FIM')
                    % FIM coupling
                    %%%%FIM settings
                    switch obj.SimulatorSettings.Formulation
                        case {'Geothermal_SinglePhase','Geothermal_MultiPhase'}
                            NLSolver = NL_Solver_geothermal();
                        otherwise
                            NLSolver = NL_Solver();
                    end
                    NLSolver.SystemBuilder = fim_system_builder();
                    switch obj.SimulatorSettings.DiscretizationModel
                        case ('ADM')
                            % Build a different convergence cheker and a proper LS for ADM
                            switch obj.SimulatorSettings.Formulation
                                case {'Geothermal_SinglePhase'}
                                    ConvergenceChecker = convergence_checker_ADM_geothermal_singlephase();
                                case {'Geothermal_MultiPhase'}
                                    ConvergenceChecker = convergence_checker_ADM_geothermal_multiphase();
                                case {'Molar_Energy'}
                                    ConvergenceChecker = convergence_checker_ADM_thermal();
                                otherwise
                                    ConvergenceChecker = convergence_checker_ADM();
                            end
                            
                            ConvergenceChecker.OperatorsAssembler = operators_assembler_fim(obj.NofEq);
                            NLSolver.LinearSolver = linear_solver_ADM(obj.SimulatorSettings.LinearSolver, 1e-6, 500);
                            NLSolver.LinearSolver.OperatorsAssembler = operators_assembler_fim(obj.NofEq);

                            if obj.SimulatorSettings.ADMSettings.DLGR
                                % it will change perm during ADM
                                % simulaiton to use upscaled ones
                                NLSolver.LinearSolver.DLGR = 1;
                            end
                        otherwise
                            switch (obj.SimulatorSettings.Formulation)
                                case('Molar')
                                    ConvergenceChecker = convergence_checker_finescale_molar();
                                case{'Geothermal_SinglePhase'}
                                    ConvergenceChecker = convergence_checker_finescale_geothermal_singlephase();
                                case{'Geothermal_MultiPhase'}
                                    ConvergenceChecker = convergence_checker_finescale_geothermal_multiphase();
                                case('Molar_Energy')
                                    ConvergenceChecker = convergence_checker_finescale_molar_energy();
                                otherwise
                                    ConvergenceChecker = convergence_checker_finescale();
                            end
                            NLSolver.LinearSolver = linear_solver(obj.SimulatorSettings.LinearSolver, 1e-6, 500);
                    end
                    NLSolver.SystemBuilder.NumberOfEq = obj.NofEq;
                    NLSolver.MaxIter = obj.SimulatorSettings.MaxIterations;
                    ConvergenceChecker.NumberOfEq  = obj.NofEq;
                    if length(obj.SimulatorSettings.ResidualTolerances) < obj.NofEq
                        obj.SimulatorSettings.ResidualTolerances = obj.SimulatorSettings.ResidualTolerances(1) * ones(1, obj.NofEq);
                    end
                    if length(obj.SimulatorSettings.SolutionTolerances) < obj.NofEq
                        obj.SimulatorSettings.SolutionTolerances = obj.SimulatorSettings.SolutionTolerances(1) * ones(1, obj.NofEq);
                    end
                    ConvergenceChecker.ResidualTol = obj.SimulatorSettings.ResidualTolerances;
                    ConvergenceChecker.SolutionTol = obj.SimulatorSettings.SolutionTolerances;
                    switch (obj.SimulatorSettings.Formulation)
                        case('Immiscible')
                            ConvergenceChecker.NormCalculator = norm_calculator_immiscible();
                        case{'Geothermal_SinglePhase'}
                            ConvergenceChecker.NormCalculator = norm_calculator_geothermal_singlephase();
                        case{'Geothermal_MultiPhase'}
                            ConvergenceChecker.NormCalculator = norm_calculator_geothermal_multiphase();
                        case{'Molar'}
                            ConvergenceChecker.NormCalculator = norm_calculator_compositional();
                        case{'Molar_Energy'}
                            ConvergenceChecker.NormCalculator = norm_calculator_compositional_energy();
                    end
                    NLSolver.AddConvergenceChecker(ConvergenceChecker);
                    % Build FIM Coupling strategy
                    if obj.SimulatorSettings.LTS
                        LTSNLSolver = NL_Solver();
                        LTSNLSolver.SystemBuilder = LTS_fim_system_builder();
                        LTSNLSolver.SystemBuilder.NumberOfEq = obj.NofEq;
                        LTSNLSolver.MaxIter = obj.SimulatorSettings.MaxIterations;
                        LTSNLSolver.LinearSolver = LTS_linear_solver(obj.SimulatorSettings.LinearSolver, 1e-6, 500);
                        LTSNLSolver.AddConvergenceChecker(ConvergenceChecker);
                        Coupling = LTS_FIM_Strategy('FIM', NLSolver, LTSNLSolver, obj.SimulatorSettings.LTStol);
                    else
                        Coupling = FIM_Strategy('FIM', NLSolver);
                    end
                    %% Sequential coupling
                case('Sequential')
                    % Sequential coupling
                    if obj.SimulatorSettings.LTS == 1 && strcmp(obj.SimulatorSettings.DiscretizationModel,'ADM') == 1
                        Coupling =  LTS_ADM_Adaptive_Sequential_Strategy('Sequential', obj.SimulatorSettings.LTStol);
                        
                        LTStransportsolver = LTS_NL_Solver();
                        LTStransportsolver.LinearSolver = LTS_linear_solver_ADM(obj.SimulatorSettings.LinearSolver, 1e-6, 500);
                        
                        LTStransportsolver.MaxIter = obj.SimulatorSettings.TransportSolver.MaxIter;
                        ConvergenceChecker = convergence_checker_transport();
                        %ConvergenceChecker = LTS_convergence_checker_transport_ADM();
                        ConvergenceChecker.ResidualTol = obj.SimulatorSettings.TransportSolver.Tol;
                        ConvergenceChecker.SolutionTol = 1e-3;
                        LTStransportsolver.AddConvergenceChecker(ConvergenceChecker);
                        LTStransportsolver.ConvergenceChecker.adm = 1;
                        LTStransportsolver.ConvergenceChecker.lts = 1;
                        if obj.capillarymodel == 1
                            LTStransportsolver.ConvergenceChecker.capillary = 1;
                        end
                        
                        LTStransportsolver.SystemBuilder = LTStransport_system_builder();
                        LTStransportsolver.LinearSolver.OperatorsAssembler = LTS_operators_assembler_seq(2, obj.NofEq);
                        Coupling.AddLTSTransportSolver(LTStransportsolver);
                        
                        % OUTER ITERATION (for LTS and not)
                        Coupling.ConvergenceChecker = convergence_checker_outer();
                        Coupling.ConvergenceChecker.ResidualTol = obj.SimulatorSettings.TransportSolver.Tol;
                        Coupling.MaxIter = obj.SimulatorSettings.MaxIterations;
                        % pressuresolver = incompressible_pressure_solver();
                        pressuresolver = NL_Solver();
                        pressuresolver.MaxIter = 15;
                        ConvergenceChecker = convergence_checker_pressure();
                        ConvergenceChecker.ResidualTol = 1e-6;
                        pressuresolver.AddConvergenceChecker(ConvergenceChecker);
                        pressuresolver.SystemBuilder = pressure_system_builder();
                        
                        
                        % at the moment we use a fine pressure solver
                        %                     switch obj.SimulatorSettings.DiscretizationModel
                        %                         case ('ADM')
                        %                             pressuresolver.LinearSolver = linear_solver_ADM(obj.SimulatorSettings.LinearSolver, 1e-6, 500);
                        %                             pressuresolver.LinearSolver.OperatorsAssembler = operators_assembler_seq(1, 1);
                        %                             pressuresolver.ConvergenceChecker.adm = 1;
                        %                             pressuresolver.ConvergenceChecker.OperatorsAssembler = operators_assembler_seq(1,1);
                        %                         case('MMs')
                        %                             pressuresolver.LinearSolver = linear_solver_MMs(obj.SimulatorSettings.LinearSolver, 1e-6, 500);
                        %                             MMsSettings = obj.SimulatorSettings.MMsSettings;
                        %                             pressuresolver.LinearSolver.MSFE = MMsSettings.MSFE;
                        %                             if MMsSettings.CorrectionFunctions
                        %                                 pressuresolver.LinearSolver.CorrectionFunctions = true;
                        %                             end
                        %                         otherwise
                        
                        pressuresolver.LinearSolver = linear_solver(obj.SimulatorSettings.LinearSolver, 1e-6, 500);
                        
                        if obj.incompressible
                            pressuresolver.ConvergenceChecker.Incompressible = 1;
                        end
                        Coupling.AddPressureSolver(pressuresolver);
                        
                        transportsolver = NL_Solver();
                        transportsolver.MaxIter = obj.SimulatorSettings.TransportSolver.MaxIter;
                        ConvergenceChecker = convergence_checker_transport();
                        ConvergenceChecker.ResidualTol = obj.SimulatorSettings.TransportSolver.Tol;
                        ConvergenceChecker.SolutionTol = 0.2;
                        
                        
                        transportsolver.AddConvergenceChecker(ConvergenceChecker);
                        transportsolver.SystemBuilder = transport_system_builder();
                        Coupling.ConvergenceChecker = convergence_checker_outer();
                        Coupling.ConvergenceChecker.SolutionTol = obj.SimulatorSettings.Tolerance;
                        
                        transportsolver.LinearSolver = linear_solver_ADM(obj.SimulatorSettings.LinearSolver, 1e-2, 500);
                        transportsolver.LinearSolver.OperatorsAssembler = LTS_operators_assembler_seq(2, obj.NofEq);
                        transportsolver.ConvergenceChecker.adm = 1;
                        if obj.capillarymodel == 1
                            transportsolver.ConvergenceChecker.capillary = 1;
                        end
                        transportsolver.ConvergenceChecker.OperatorsAssembler = LTS_operators_assembler_seq(2, obj.NofEq);
                        
                        Coupling.AddTransportSolver(transportsolver);
                    else
                        if obj.SimulatorSettings.LTS == 1
                            switch (obj.SimulatorSettings.LTSCriterion)
                                case('Fixed')
                                    Coupling = LTS_Sequential_Strategy('Sequential', obj.SimulatorSettings.LTStol);
                                case('Adaptive')
                                    Coupling = LTS_Adaptive_Sequential_Strategy('Sequential', obj.SimulatorSettings.LTStol);
                            end
                            LTStransportsolver = LTS_NL_Solver();
                            LTStransportsolver.MaxIter = obj.SimulatorSettings.TransportSolver.MaxIter;
                            ConvergenceChecker = convergence_checker_transport();
                            ConvergenceChecker.ResidualTol = obj.SimulatorSettings.TransportSolver.Tol;
                            ConvergenceChecker.SolutionTol = 2e-2;
                            LTStransportsolver.AddConvergenceChecker(ConvergenceChecker);
                            if obj.capillarymodel == 1
                                LTStransportsolver.ConvergenceChecker.capillary = 1;
                            end
                            LTStransportsolver.SystemBuilder = LTStransport_system_builder();
                            LTStransportsolver.LinearSolver = LTS_linear_solver(obj.SimulatorSettings.LinearSolver, 1e-6, 500);
                            Coupling.AddLTSTransportSolver(LTStransportsolver);
                        else
                            Coupling = Sequential_Strategy('Sequential');
                        end
                        % OUTER ITERATION (for LTS and not)
                        Coupling.ConvergenceChecker = convergence_checker_outer();
                        Coupling.ConvergenceChecker.ResidualTol = obj.SimulatorSettings.TransportSolver.Tol;
                        Coupling.MaxIter = obj.SimulatorSettings.MaxIterations;
                        
                        pressuresolver = NL_Solver();
                        pressuresolver.MaxIter = 15;
                        ConvergenceChecker = convergence_checker_pressure();
                        ConvergenceChecker.ResidualTol = 1e-6;
                        pressuresolver.AddConvergenceChecker(ConvergenceChecker);
                        pressuresolver.SystemBuilder = pressure_system_builder();
                        switch obj.SimulatorSettings.DiscretizationModel
                            case ('ADM')
                                pressuresolver.LinearSolver = linear_solver(obj.SimulatorSettings.LinearSolver, 1e-6, 500);
                                
                                %                                 pressuresolver.LinearSolver = linear_solver_ADM(obj.SimulatorSettings.LinearSolver, 1e-6, 500);
                                %                                 pressuresolver.LinearSolver.OperatorsAssembler = operators_assembler_seq(1, 1);
                                %                                 pressuresolver.ConvergenceChecker.adm = 1;
                                %                                 pressuresolver.ConvergenceChecker.OperatorsAssembler = operators_assembler_seq(1,1);
                            case('MMs')
                                pressuresolver.LinearSolver = linear_solver_MMs(obj.SimulatorSettings.LinearSolver, 1e-6, 500);
                                MMsSettings = obj.SimulatorSettings.MMsSettings;
                                pressuresolver.LinearSolver.MSFE = MMsSettings.MSFE;
                                if MMsSettings.CorrectionFunctions
                                    pressuresolver.LinearSolver.CorrectionFunctions = true;
                                end
                            otherwise
                                pressuresolver.LinearSolver = linear_solver(obj.SimulatorSettings.LinearSolver, 1e-6, 1);
                        end
                        if obj.incompressible
                            pressuresolver.ConvergenceChecker.Incompressible = 1;
                        end
                        Coupling.AddPressureSolver(pressuresolver);
                        switch (obj.SimulatorSettings.TransportSolver.Type)
                            case('IMPSAT')
                                transportsolver = NL_Solver();
                                transportsolver.MaxIter = obj.SimulatorSettings.TransportSolver.MaxIter;
                                ConvergenceChecker = convergence_checker_transport();
                                ConvergenceChecker.ResidualTol = obj.SimulatorSettings.TransportSolver.Tol;
                                ConvergenceChecker.SolutionTol = 1e-3;
                                transportsolver.AddConvergenceChecker(ConvergenceChecker);
                                if obj.capillarymodel == 1
                                    transportsolver.ConvergenceChecker.capillary = 1;
                                end
                                transportsolver.SystemBuilder = transport_system_builder();
                                Coupling.ConvergenceChecker = convergence_checker_outer();
                                Coupling.ConvergenceChecker.SolutionTol = obj.SimulatorSettings.Tolerance;
                                switch obj.SimulatorSettings.DiscretizationModel
                                    case('ADM')
                                        transportsolver.LinearSolver = linear_solver_ADM(obj.SimulatorSettings.LinearSolver, 1e-6, 500);
                                        transportsolver.LinearSolver.OperatorsAssembler = operators_assembler_seq(2, obj.NofEq);
                                        transportsolver.ConvergenceChecker.adm = 1;
                                        transportsolver.ConvergenceChecker.OperatorsAssembler = operators_assembler_seq(2, obj.NofEq);
                                    otherwise
                                        transportsolver.LinearSolver = linear_solver(obj.SimulatorSettings.LinearSolver, 1e-6, 500);
                                end
                            otherwise
                                transportsolver = explicit_transport_solver();
                                transportsolver.SystemBuilder = explicit_transport_system_builder();
                                Coupling.ConvergenceChecker = convergence_checker_impes();
                        end
                        Coupling.AddTransportSolver(transportsolver);
                    end
                case('SinglePhase')
                    % Single phase coupling
                    Coupling = SinglePhase_Strategy('SinglePhase');
                    pressuresolver = NL_Solver();
                    pressuresolver.MaxIter = 15;
                    ConvergenceChecker = convergence_checker_pressure();
                    ConvergenceChecker.ResidualTol = obj.SimulatorSettings.ResidualTolerances;
                    ConvergenceChecker.SolutionTol = obj.SimulatorSettings.SolutionTolerances;
                    ConvergenceChecker.NumberOfEq = obj.NofEq;
                    pressuresolver.AddConvergenceChecker(ConvergenceChecker);
                    pressuresolver.SystemBuilder = pressure_system_builder();
                    switch obj.SimulatorSettings.DiscretizationModel
                        case('ADM')
                            pressuresolver.LinearSolver = linear_solver_ADM(obj.SimulatorSettings.LinearSolver, 1e-6, 500);
                            pressuresolver.LinearSolver.OperatorsAssembler = operators_assembler_fim(obj.NofEq);
                            pressuresolver.ConvergenceChecker.adm = 1;
                            pressuresolver.ConvergenceChecker.OperatorsAssembler = operators_assembler_fim(obj.NofEq);
                        case ('MMs')
                            pressuresolver.LinearSolver = linear_solver_MMs(obj.SimulatorSettings.LinearSolver, 1e-6, 500);
                            pressuresolver.LinearSolver.OperatorsAssembler = operators_assembler_seq(1, 1);
                            pressuresolver.LinearSolver.MSFE = obj.SimulatorSettings.MMsSettings.MSFE;
                            if obj.SimulatorSettings.MMsSettings.CorrectionFunctions
                                pressuresolver.LinearSolver.CorrectionFunctions = true;
                            end
                        otherwise
                            pressuresolver.LinearSolver = linear_solver(obj.SimulatorSettings.LinearSolver, 1e-6, 500);
                    end
                    if obj.incompressible
                        Coupling.Incompressible = 1;
                        pressuresolver.ConvergenceChecker.Incompressible = 1;
                    end
                    Coupling.AddPressureSolver(pressuresolver);
            end
            if obj.SimulatorSettings.LTS == 1 || obj.SimulatorSettings.FixedStep == 1
                Coupling.TimeStepSelector = fixed_timestep_selector(obj.SimulatorSettings.cfl, obj.SimulatorSettings.MinMaxdt(1), obj.SimulatorSettings.MinMaxdt(2), obj.SimulatorSettings.MaxIterations);
            else
                Coupling.TimeStepSelector = timestep_selector(obj.SimulatorSettings.cfl, obj.SimulatorSettings.MinMaxdt(1), obj.SimulatorSettings.MinMaxdt(2), obj.SimulatorSettings.MaxIterations);
            end
            Coupling.TimeStepSelector.TotalTime = obj.SimulationInput.TotalTime;
            TimeDriver.AddCouplingStrategy(Coupling);
            if isempty(obj.SimulatorSettings.StopCriterion)
                obj.SimulatorSettings.StopCriterion = 'MAX TIME';
            end
            
            switch(obj.SimulatorSettings.StopCriterion)
                case('MAX TIME')
                    end_of_sim_eval = end_of_sim_evaluator(obj.SimulationInput.TotalTime, obj.SimulatorSettings.MaxNumTimeSteps);
                case('COMPONENT CUT')
                    end_of_sim_eval = end_of_sim_evaluator_gascut(obj.SimulationInput.TotalTime, obj.SimulatorSettings.MaxNumTimeSteps, 0.2);
                case('PV INJECTED')
                    end_of_sim_eval = end_of_sim_evaluator_PVInjected(obj.SimulationInput.TotalTime, obj.SimulatorSettings.MaxNumTimeSteps, 3);
            end
            TimeDriver.AddEndOfSimEvaluator(end_of_sim_eval);
        end
        function Summary = BuildSummary(obj, simulation)
            %%%%%%%%%%%%%%% BuildObjects for OUTPUT%%%%%%%%%
            switch(obj.SimulatorSettings.CouplingType)
                case('FIM')
                    if obj.SimulatorSettings.LTS == 1
                        CouplingStats = LTS_FIM_Stats(obj.SimulatorSettings.MaxNumTimeSteps, 'LTS_FIM');
                    else
                        CouplingStats = FIM_Stats(obj.SimulatorSettings.MaxNumTimeSteps, 'FIM');
                    end                    
                case('Sequential')
                    if obj.SimulatorSettings.LTS == 1
                        CouplingStats = LTS_Sequential_Stats(obj.SimulatorSettings.MaxNumTimeSteps, 'LTS_Sequential');
                    else
                        CouplingStats = Sequential_Stats(obj.SimulatorSettings.MaxNumTimeSteps, 'Sequential');
                    end
                case('SinglePhase')
                    CouplingStats = SinglePhase_Stats(obj.SimulatorSettings.MaxNumTimeSteps, 'SinglePhase');
            end
            wellsData = wells_data(obj.SimulatorSettings.MaxNumTimeSteps, simulation.FluidModel.NofPhases, simulation.FluidModel.NofComponents, simulation.ProductionSystem.Wells);
            switch (obj.SimulatorSettings.DiscretizationModel)
                case('ADM')
                    Summary = Run_Summary_ADM(obj.SimulatorSettings.MaxNumTimeSteps, CouplingStats, wellsData, simulation.DiscretizationModel.maxLevel(1)); % Only reservoir for now
                otherwise
                    Summary = Run_Summary(obj.SimulatorSettings.MaxNumTimeSteps, CouplingStats, wellsData);
            end
        end
    end
end
