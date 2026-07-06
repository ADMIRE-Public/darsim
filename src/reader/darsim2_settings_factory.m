%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DARSim
% Author: Artur Castiel Reis de Souza
% TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef darsim2_settings_factory < handle
    %darsim2_factory translate blueprint into darsim2 compatible input legacy code 
    % darsim2 requires the configuration settings to be organized in a specific way.
	% this class provides a compatibility layer to allow the use of the new reading routines.
	properties
        SimulationInput
        SimulatorSettings
        SimulatorFolder
    end
    properties (Access = private)
        blueprint
    end

    methods
        function obj = darsim2_settings_factory(blueprint)            
            obj.SimulatorFolder = blueprint.get_simulator_folder;
            % SimulationInput describes the parameters of the physics of
            % of the reservoir.
            obj.SimulationInput = ...
                obj.translate_simulation_input(blueprint);
            % SimulationSettings describes the paramters that simulator 
            % such as tolerances, residues, visualization tools and so on.
            obj.SimulatorSettings = ...
                obj.translate_simulation_settings(blueprint.SETTINGS,...
                     obj.SimulationInput.FracturesProperties.NumOfFrac);
        end
    end
    methods (Access = private)
        function SimulationSettings = ...
                translate_simulation_settings(obj,blueprint, num_frac)
            SimulationSettings = struct();
            SimulationSettings.MaxNumTimeSteps = blueprint.TIMESTEPS;
            SimulationSettings.MinMaxdt = get_variable(blueprint.MINMAXDT);
            SimulationSettings.FixedStep = blueprint.Fix;
            SimulationSettings.LinearSolver = blueprint.LINEARSOLVER;
            SimulationSettings.cfl = blueprint.CFL;
            SimulationSettings.Reports = blueprint.REPORTS;
            SimulationSettings.NumOfPreviousReports = ...
                                blueprint.NUM_PREVIOUS_REPORTS;            
            SimulationSettings.LTS = false;
            SimulationSettings.LTSPlot = false;            
            SimulationSettings.half_precision = blueprint.HALF_PRECISION;
            SimulationSettings.parallel = blueprint.PARALLEL;
            if iscell(blueprint.FIM)
                [SimulationSettings.CouplingType, ...
                    SimulationSettings.MaxIterations, ...
                    SimulationSettings.ResidualTolerances, ...
                    SimulationSettings.SolutionTolerances, ...
                    SimulationSettings.cfl] = build_FIM_config(...
                    blueprint.FIM, SimulationSettings.cfl);
            end
            [SimulationSettings.DiscretizationModel,...
                                    SimulationSettings.ADMSettings] = ...
                           build_DiscretizationModel(blueprint,num_frac,obj.SimulationInput.Directory);
            SimulationSettings.Formulation = ...
                build_FormulationModel_SettingsFile(blueprint, ...
                obj.SimulationInput.FluidProperties.FluidModel);                
            SimulationSettings.StopCriterion = blueprint.STOP_CRITERION;
            SimulationSettings.plotting = ...
                            build_VisualizationModel(blueprint);
            SimulationSettings.WriterSettings = build_WriterSettings( ...
                SimulationSettings.DiscretizationModel,...
                SimulationSettings.LTS);
        end


        function SimulationInput = ...
                            translate_simulation_input(~, blueprint)
            SimulationInput = build_SimulationInput(blueprint);
            SimulationInput.ReservoirProperties = ...
                                build_ReservoirProperties(blueprint);
            SimulationInput.FluidProperties = ...
                                build_FluidProperties(blueprint);
            SimulationInput.WellsInfo = build_well_info(...
                    blueprint.WELLS.WELL,...
                    SimulationInput.ReservoirProperties.Grid.N);
            SimulationInput.InitialConditions = ...
                                build_InitialConditions(blueprint);                        
            SimulationInput.RockProperties = ...
                                build_RockProperties(blueprint, SimulationInput.ReservoirProperties.Grid.N_ActiveCells, SimulationInput.FluidProperties.NofPhases);
            SimulationInput.FracturesProperties = ...
                    build_FracturesProperties(blueprint.FRACTURE, ...
                    blueprint.GENERAL.FRACTURED);            
        end
    end
end


function WriterSettings = build_WriterSettings(DiscretizationModel, LTS)
WriterSettings = struct();
switch DiscretizationModel
    case "FIM"
        WriterSettings.NTimer = 4;
        WriterSettings.NStats = 4;       
    case "Sequential"
        if LTS
            WriterSettings.NTimer = 2;
            WriterSettings.NStats = 5;
        else
            WriterSettings.NTimer = 2;
            WriterSettings.NStats = 2;
        end
    case "SinglePhase"
        WriterSettings.NTimer = 1;
        WriterSettings.NStats = 1; 
    otherwise
        WriterSettings.NTimer = 4;
        WriterSettings.NStats = 4; 
end

end


function plotting = build_VisualizationModel(blueprint)
    plotting = struct();
    plotting.Software = blueprint.OUTPUT{1,1};
    plotting.Format =  blueprint.OUTPUT{2,1};%'ASCII';
    plotting.PlotInterfaces = blueprint.PLOT_INTERFACES;
    if ~isfield(blueprint, 'PLOT_BASIS_FUNCTIONS')
        plotting.PlotBasisFunctions = true;
    else
        plotting.PlotBasisFunctions = false;
    end
end

function [DiscretizationModel,ADMSettings] = ...
                    build_DiscretizationModel(blueprint, num_frac, Directory)
    ADMSettings = struct();  
    if ~blueprint.ADM
        DiscretizationModel = "FS";
    else
        DiscretizationModel = "ADM";
        ADMSettings.CoarseningSwitch = zeros(1+num_frac, 1);
        ADMSettings.CoarseningSwitch(1) = 1;

        ADMSettings.DLGR = blueprint.DLGR;
        ADMSettings.HInterpolator = blueprint.HYPERBOLIC_INTERPOLATOR;
        ADMSettings.PInterpolator = blueprint.PRESSURE_INTERPOLATOR;
        ADMSettings.Coarsening = blueprint.COARSENING_RATIOS';
        ADMSettings.GridSelCriterion = blueprint.COARSENING_CRITERION;
        ADMSettings.tol = blueprint.TOLERANCE;
        ADMSettings.key = blueprint.VARIABLE;
        ADMSettings.BFtype = blueprint.PRESSURE_INTERPOLATOR{2};
        ADMSettings.BF_MaxContrast = blueprint.BASIS_FUNCTION_MAX_CONTRAST;
        ADMSettings.pEDFM_MaxContrast = blueprint.pEDFM_MAX_CONTRAST;
        ADMSettings.BF_alpha = ...
            blueprint.BASIS_FUNCTION_REGION_OF_INFLUENCE_CUT;
        cx = blueprint.COARSENING_RATIOS(1);
        cy = blueprint.COARSENING_RATIOS(2);
        cz = blueprint.COARSENING_RATIOS(3);

        ADMSettings.maxLevel = zeros(1 + num_frac,1);
        ADMSettings.maxLevel(1) = blueprint.LEVELS;
        for L = 1:ADMSettings.maxLevel(1)
                    ADMSettings.Coarsening(1,:,L) = [cx, cy, cz].^L; %Coarsening Factors: Cx1, Cy1; Cx2, Cy2; ...; Cxn, Cyn;
        end
        ADMSettings.BF_MaxContrast = blueprint.BASIS_FUNCTION_MAX_CONTRAST;
        ADMSettings.pEDFM_MaxContrast = blueprint.pEDFM_MAX_CONTRAST;
        ADMSettings.BF_alpha = 0;
        if ~isfield(blueprint, 'Coupling')
            ADMSettings.Coupling = 'DECOUPLED';
        end
        FractureFile = fullfile(Directory, 'Fracture_Output.txt');
        if exist(FractureFile, 'file') == 2
            fileID = fopen(FractureFile, 'r');
            matrix = textscan(fileID, '%s', 'Delimiter', '\n');
            FractureMatrix = matrix{1};
            fclose(fileID);
            temp = strfind(FractureMatrix, 'PROPERTIES');
            frac_index = find(~cellfun('isempty', temp));
            for f = 1 : num_frac
                frac_info_split = strsplit(FractureMatrix{frac_index(f)},' ');
                ADM_temp = regexprep(frac_info_split{9},' ' ,'');
                ADM_temp = strsplit(ADM_temp, { '[' , ',' , ']' });
                ADM_temp = [ str2double(ADM_temp(2)) , str2double(ADM_temp(3)) , str2double(ADM_temp(4)) , str2double(ADM_temp(5)) ];
                ADMSettings.CoarseningSwitch(1+f) = ADM_temp(1);
                ADMSettings.maxLevel(1+f) = ADM_temp(2);
                for L = 1:ADMSettings.maxLevel(1)
                    if L <= ADMSettings.maxLevel(1+f)
                        ADMSettings.Coarsening(1+f,:,L) = [ADM_temp(3), ADM_temp(4), 1].^L;
                    else
                        ADMSettings.Coarsening(1+f,:,L) = [ADM_temp(3), ADM_temp(4), 1].^ADMSettings.maxLevel(1+f);
                    end
                end
            end
        else
            return;
        end
    end
end

function [CouplingType, MaxIterations, ResidualTolerances, ... 
                SolutionTolerances, cfl] = build_FIM_config(FIM, cfl)
        CouplingType = 'FIM';
        MaxIterations = FIM{1};
        ResidualTolerances = get_variable(FIM{2});
        SolutionTolerances = get_variable(FIM{3});
        if size(FIM,1) > 3
           cfl = get_variable(FIM{4});
        end
end


function Formulation = ...
            build_FormulationModel_SettingsFile(blueprint, formulation)
    if isfield(blueprint,"FORMULATION")
       Formulation = blueprint.FORMULATION;
    else
        formulation = convertCharsToStrings(formulation);
        switch formulation
            case {"SinglePhase","Immiscible"}
                Formulation = "Immiscible";
            case "Geothermal_SinglePhase"
                Formulation = "Geothermal_SinglePhase";
            case "Geothermal_MultiPhase"
                Formulation = "Geothermal_MultiPhase";
            case {"Compositional","BlackOil"}
                Formulation = "Molar";
            case {"Compositional_Energy"}
                Formulation = "Molar_Energy";
            otherwise
                error("Error DARSim could not recognized formulation.")
        end
    end
end


function InitialConditions = build_InitialConditions(blueprint)
    InitialConditions = struct();
    if isfield(blueprint.BC, 'INITIAL_STATE')
        InitialConditions.Include = true;
        InitialConditions.LoadedData = load(fullfile(blueprint.get_folder_file(), blueprint.BC.INITIAL_STATE));
        InitialConditions.LoadedData(:, 1) =[];
    else
        InitialConditions.Include = false;
        InitialConditions.Saturation_1 = blueprint.BC.INITIAL_SATURATION_1;
        InitialConditions.Saturation_2 = blueprint.BC.INITIAL_SATURATION_2;
        InitialConditions.Composition_1 = blueprint.BC.INITIAL_COMPOSITION_1;
        InitialConditions.Composition_2 = blueprint.BC.INITIAL_COMPOSITION_2;
        if isfield(blueprint.BC, 'INITIAL_COMPOSITION_3')
            InitialConditions.Composition_3 = blueprint.BC.INITIAL_COMPOSITION_3;
        end
    end
        InitialConditions.Pressure = blueprint.BC.INITIAL_PRESSURE;
        InitialConditions.Enthalpy = blueprint.BC.INITIAL_ENTHALPY;
        InitialConditions.Temperature = blueprint.BC.INITIAL_TEMPERATURE;
end

function RockProperties = build_RockProperties(blueprint, N_active_cells, n_phases)
    RockProperties = struct();
    RockProperties.Density = blueprint.ROCK.ROCK_DENSITY;
    RockProperties.Compressibility = blueprint.ROCK.ROCK_COMPRESSIBILITY;
    RockProperties.HeatConductivity = blueprint.ROCK.ROCK_HEAT_CONDUCTIVITY;
    RockProperties.SpecificHeat = blueprint.ROCK.ROCK_SPECIFIC_HEAT;
    diffusion = blueprint.ROCK.DIFFUSION_CONSTANTS;
    dispersivity = blueprint.ROCK.DISPERSIVITY;
    if isnumeric(diffusion)
        RockProperties.Diffusion = diffusion' .* ones(N_active_cells, 1);
    elseif iscell(diffusion)
        for i = 1: n_phases
            temp = diffusion{i};
            if isnumeric(temp)
                RockProperties.Diffusion(:,i) = temp * ones(N_active_cells, 1);
            elseif ischar(temp) || (isstring(temp) && isscalar(temp))
                DiffusionFile = fullfile(blueprint.get_folder_file, char(temp));
                RockProperties.Diffusion(:,i) = read_unique_input_file_diffusion(DiffusionFile, N_active_cells);
            else
                error('Unsupported DIFFUSION_CONSTANTS entry type at position %d.', i);
            end
        end
    end

    if isnumeric(dispersivity)
        RockProperties.Dispersivity = dispersivity * ones(N_active_cells, 1);
    elseif ischar(dispersivity) || (isstring(dispersivity) && isscalar(dispersivity))
        DispersivityFile = fullfile(blueprint.get_folder_file, char(dispersivity));
        RockProperties.Dispersivity = read_unique_input_file_dispersion(DispersivityFile, N_active_cells);
    else
        error('Unsupported DISPERSIVITY type.')
    end

end

function diff = read_unique_input_file_diffusion(unique_input_file, n_total)
    diff = load(unique_input_file);

    if size(diff,1) ~= n_total
        error('The size of the diffusion file is not consistent with the mesh.');
    end    
end

function disp = read_unique_input_file_dispersion(unique_input_file, n_total)
    disp = load(unique_input_file);

    if size(disp,1) ~= n_total
        error('The size of the dispersion file is not consistent with the mesh.');
    end    
end

function SimulationInput = build_SimulationInput(blueprint)
SimulationInput = struct();
SimulationInput.Directory = blueprint.get_folder_file;
SimulationInput.ProblemName  = blueprint.GENERAL.TITLE;
SimulationInput.TotalTime = ...
    blueprint.GENERAL.TOTAL_TIME.total_time;
end







