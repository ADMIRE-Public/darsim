%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DARSim
% Author: Artur Castiel Reis de Souza
% TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef blueprint_factory < handle
	% blueprint_factory organize the input files processed with reader_darsim3 into
	% a collection of blueprints. 
	% Each blueprint is built based on the variables in the properties and the
	% default file it points.
	% 
	% Example:
	% 	properties
	% 	VAR = "var.bpdat"
	% 
	% After the intialization, the blueprint_factory object will have a properti called VAR containing the labels written in var.bpdat.	
	
    properties (Access=private)
        simulator_path = define_simulator_path
    end    
    properties
        %% to add new blueprint to be build, we only need to add the field
        % and the file containing the description added to the blueprint
        % configuratin folder
        GENERAL = "general.bpdat";
        GRID = "grid.bpdat";
        ROCK = "rock.bpdat";
        RELPERM_MODEL = "relperm_model.bpdat";
        FLUID = "fluid.bpdat";
        WELLS = "wells.bpdat";
        BC = "bc.bpdat"
        COMPONENTS = "components.bpdat";
        SETTINGS = "settings.bpdat";
        % Add other FIELDS only before OTHERS.
        % The order of the remaining fields is important.
        OTHERS
        FRACTURE
    end
    methods
        function obj = blueprint_factory(folder)
            fracture_folder = folder.FractureFile;
            folder.FractureFile = [];
            obj.initialize_fields(obj.simulator_path);
            all_fields = obj.generate_all_used_fields();
            obj.generate_blueprints(folder);
            if folder.isFractured
                obj.FRACTURE = fracture_folder;
            end
            obj.OTHERS = obj.generate_others(all_fields,folder);
            obj.get_folder_file(folder);
            obj.get_simulator_folder(folder);
        end
        function  [folder, file] = get_folder_file(~,varargin)
            persistent folder_memory 
            persistent file_memory
            if ~isempty(varargin)
              [folder_memory, file_memory] = varargin{1}.get_file_path();
            end
            folder = folder_memory;
            file = file_memory;
        end
        function  [folder] = get_simulator_folder(~,varargin)
            persistent folder_memory 
            if ~isempty(varargin)
              [folder_memory] = varargin{1}.SimulatorPath();
            end
            folder = folder_memory;
        end   
    
    end
    methods (Access = private)
        function initialize_fields(obj, simulator_path)
            simulator_path = fullfile(simulator_path, "reader", ...
                 "default-configuration-settings", "blueprint_settings");
            props = convertCharsToStrings(properties(obj)');
            props(end-1:end) = [];
            for label = props
                fields = obj.(label);
                obj.(label) = read_blueprint_fields(simulator_path,fields);
            end
        end
        function generate_blueprints(obj,folder)
            props = convertCharsToStrings(properties(obj)');
            props(end-1:end) = [];
            for label = props
                fields = obj.(label);
               obj.(label) = extract_blueprint(folder,fields);
            end
        end
        function others = generate_others(~, all_fields, folder)
            % generate OTHER blueprint with other tags not found in the
            % blueprint configuration
            all_fields_in_folder = [];
            all_upper_fields = convertCharsToStrings(properties(folder)');
            all_upper_fields = ...
                all_upper_fields(~ismember(...
                all_upper_fields,"FractureFile"));
            for index = 1:length(all_upper_fields)
                local_upper_field = all_upper_fields(index);
                local_upper_value = folder.(local_upper_field);
                if isstruct(local_upper_value)
                    fields = ...
                     convertCharsToStrings(fieldnames(local_upper_value)');
                    all_fields_in_folder = [all_fields_in_folder fields]; %#ok<*AGROW>
                end
            end
            other_fields = setdiff(all_fields_in_folder,all_fields);
            others = extract_blueprint(folder,other_fields);
        end

        function all_fields = generate_all_used_fields(obj)
            props = convertCharsToStrings(properties(obj)');
            props(end-1:end) = [];
            all_fields = [];
            for label = props
               fields = obj.(label);
               all_fields = [all_fields fields];
            end           
            all_fields = unique(all_fields);
        end
    end
end