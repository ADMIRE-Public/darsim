%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel;
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef reader_darsim3
    %reader_darsim3 base class that reads and transform input files into
    % matlab structures stores on its properties. The routines also add default 
	% values for missing labels.
    properties
        InputFile
        SettingsFile
        FractureFile
		% boolean that checks if the reservoir is fractured.
        isFractured
        SimulatorPath = define_simulator_path()
    end
    properties (Access=private)
        Directory
        File
        Path
    end
    methods
        function obj = reader_darsim3(dir, file)
            obj.Directory = dir;
            obj.File = file;
            obj.Path = fullfile(dir,file);
            [obj.InputFile, obj.SettingsFile, obj.FractureFile, ...
                obj.isFractured] = run(obj);
        end

        function [folder,file] = get_file_path(obj)
            folder = obj.Directory;
            file = obj.File;
        end
        function simulator_path = getter.SimulatorPath(obj)
            simulator_path = obj.SimulatorPath;
        end
    
    end
    methods (Access = private)
        function [InputBlueprint, SettingsBlueprint, ...
                FractureBlueprint, isFractured] = run(obj)
			% reads mains INPUT File
            InputBlueprint = obj.read_input_file;
            % reads  simulationg settings input file.
			SettingsBlueprint = obj.read_settings_file;
            isFractured = InputBlueprint.FRACTURED{1};
            % if fracture flag is raised, read fracture input file.
			if isFractured
               if contains(lower(InputBlueprint.FRACTURED{2}), ...
                                            lower("output"))
                    FractureBlueprint = legacy_fracture_reader(...
                        obj.Directory , InputBlueprint.FRACTURED{2});
                    FractureBlueprint.legacy = true;
               else
                   FractureBlueprint = ...
                       obj.read_fracture_file(InputBlueprint.FRACTURED{2});
                   FractureBlueprint.legacy = false;
               end
            else
                FractureBlueprint =[];
            end
       
	   end
        function FractureBlueprint = read_fracture_file(obj, frac_file)
            FracMatrix = TXTreaderplus(fullfile(obj.Directory, frac_file));
            FractureBlueprint = convert_to_struct_txt_object(FracMatrix);
            [defaultKeys, defaultValues] = generate_defaultKeyValues(...
                              obj.SimulatorPath, "defaultFracture.dat");              
            FractureBlueprint = verify_blueprint(FractureBlueprint, ...
                defaultKeys, defaultValues);            
        end


        function InputBlueprint = read_input_file(obj)            
            InputMatrix =  TXTreaderplus(obj.Path);
            InputStruct = convert_to_struct_txt_object(InputMatrix);
            %% defining default settings and values that are imposed if 
            % missing
            [defaultKeys, defaultValues] = generate_defaultKeyValues(...
                              obj.SimulatorPath, "defaultSimulation.dat");
            %% define if labels that contains time
            timeLabels = ["TOTAL_TIME"];
            InputBlueprint = verify_blueprint(InputStruct, defaultKeys, ...
                defaultValues, timeLabels);
        end
        function SettingsBlueprint = read_settings_file(obj)
            SettingsMatrix =  TXTreaderplus(fullfile(...
                obj.Directory, "SimulatorSettings.txt"));
            SettingsBlueprint = ...
                convert_to_struct_txt_object(SettingsMatrix);
%             file_folder = fullfile(obj.SimulatorPath,"reader", ...
%                             "default-configuration-settings");
            [defaultKeys, defaultValues] = generate_defaultKeyValues(...
                   obj.SimulatorPath,"defaultSettings.dat");  
            SettingsBlueprint = verify_blueprint(SettingsBlueprint, ...
                        defaultKeys, defaultValues);            
        end
    end
end
