% Fracture_Generator Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Mousa HosseiniMehr
%TU Delft
%Created: 2017-03-10
%Last modified: 2017-03-10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pEDFM = DARSimFracGen(Directory, File)
close all; clc;

% Remove some warnings 
warning('off', 'MATLAB:singularMatrix');
warning('off', 'MATLAB:nearlySingularMatrix');

%% Set up Diary
if exist(   strcat(Directory, '/Output/Fracture_Output.txt'), 'file')
    delete( strcat(Directory, '/Output/Fracture_Output.txt')        );
end

%% Print title
disp('********************************************************************************');
disp('******************** Fracture Generator for DARSim2 Simulator ******************');
disp('********************************************************************************');
disp(char(5));
disp(['Reading input file ', File, ' from ', Directory]);
disp(char(5));

%% Build objects
% Construct pEDFM Generator class
pEDFM = pEDFM_Generator(Directory, File);
% Read Input File
pEDFM.ReadInputFile();
% Build objects

%% creating the fractures geometry e meshes
pEDFM.BuildObjects();
% Print info to screen
pEDFM.PrintInfo();

%% Run pEDFM Generator
TotalStart = tic;
%% here the intersection and the pEDFM and EDFM parameters are found
pEDFM.Run();
TotalTime = toc(TotalStart);

%% Display elapsed time
disp(['The total computation time is ' num2str(TotalTime) ' s']); fprintf('\n');

%% Output Results
pEDFM.OutputResults();

end