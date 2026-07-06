function [defaultKeys, defaultValues] = generate_defaultKeyValues(...
                                                simulator_folder, file)
%% generate default keys and value from file config folder
file_folder = fullfile(simulator_folder,"reader", ...
                            "default-configuration-settings");
file_path = fullfile(file_folder,file);
defaultInput = convert_to_struct_txt_object(TXTreaderplus(file_path));
defaultKeys = convertCharsToStrings(fieldnames(defaultInput)');
defaultValues = struct2cell(defaultInput)';
end