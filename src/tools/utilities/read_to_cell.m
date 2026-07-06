function [cellArray,lineCounter] = read_to_cell(file_path)
% file_path = fullfile(folder,'LookupTable.txt');
% Open the file for reading
fileID = fopen(file_path, 'r');
% Check if the file was successfully opened
if fileID == -1
    error('LookupTable.txt not found.');
end

% Initialize an empty cell array to store the file contents
cellArray = {};

% Initialize a line counter
lineCounter = 1;

% Read lines from the file and store them in the cell array
while ~feof(fileID)
    line = fgetl(fileID);
    cellArray{lineCounter} = line;
    lineCounter = lineCounter + 1;
end

% Close the file
fclose(fileID);
cellArray = cellArray';
end