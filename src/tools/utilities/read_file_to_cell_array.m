function lines = read_file_to_cell_array(file_path)
% Open the file for reading
fileID = fopen(file_path, 'r');
if fileID == -1
    error('Error: Unable to open the file.');
end
% Initialize the cell array to store lines
lines = cell(0, 1);
% Read the file line by line and store in the cell array
while ~feof(fileID)
    line = fgetl(fileID);
    lines{end+1, 1} = line;
end
% Close the file
fclose(fileID);
end