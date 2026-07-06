function [fields] = read_blueprint_fields(simulator_folder, file)
%file_folder = fullfile(simulator_folder,"config","blueprint_settings");
file_path = fullfile(simulator_folder,file);
fields = read_file_to_cell_array(file_path);
%readcell(file_path, "FileType","text");
fields = convertCharsToStrings(fields)';
end

