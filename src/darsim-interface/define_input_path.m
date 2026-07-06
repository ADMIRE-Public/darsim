function [folder, file] = define_input_path(arg)
flag_empty_file = isempty(arg);
if flag_empty_file
    folder = uigetdir(pwd,"Select the simulation case folder.");
    file = uigetfile({'*.txt'},"Select the Input file.",folder);
else
    arg_size = length(arg);
    if arg_size == 1
        folder = arg{1};
        validate_folder_and_file(folder,[]);
        [file, folder] = uigetfile({'*.txt'},"Select the Input file.", ...
                    folder);
    elseif arg_size == 2
        folder = arg{1};
        file = arg{2};
    end
end
validate_folder_and_file(folder,file);
end


function validate_folder_and_file(folder,file)
if isnumeric(file) && isnumeric(folder)
    if ~file && ~folder
        error("Simulation cancelled by user.")
    end
end
arg = fullfile(folder,file);
[~, ~, fExt] = fileparts(arg);
folder_flag = isempty(fExt);
if folder_flag
    if ~isfolder(arg)
        error("The specified folder does not exist.");
    end
else
    if ~isfile(arg)
        error("The specified file does not exist.");
    end
end    
end
