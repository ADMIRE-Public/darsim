%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel;
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FractureMatrix] = legacy_fracmatrix_generator(file_path)
fileID = fopen(file_path, 'r');
%fprintf('Reading fracture file ...');
matrix = textscan(fileID, '%s', 'Delimiter', '\n');
FractureMatrix = matrix{1};
fclose(fileID);
%fprintf('--> Successful.\n\n');
end
