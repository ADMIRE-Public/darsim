function current_folder = define_simulator_path
current_folder = pwd;
position = strfind(current_folder, "src");
if isempty(position)
    error("Please set src folder of DARSim as the main folder.");
end
if ~(length(current_folder)  == position +2)
    current_folder(position+3:end) = [];
end
end