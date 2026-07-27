%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DARSim
% Author: Artur Castiel Reis de Souza
% TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ## Getting Started
% To begin using this release of DARSim, no complex installation process is required. Simply open the /src folder in MATLAB and type:
%
% DARSim
%
% DARsim will present a simple user interface, containing basic information about the project, tutorials (under development), and a link to <u>Find a case</u> that will open a Dialog Box where users can define the path to the Input File.
% This interface can bypassed by simply parsing as argument the folder and file name of the simulation input, as follows:
%
%    DARSim("simulation_path", "input_file.txt")
% or DARSim("simulation_path", "input_file.txt", N) where N is the number of CPU cores being used in linear solver
%
% Once the simulation is completed, DARSim will create folder called "Output" at the location of the input file. This directory contains the visualization files and production curves of the simulation.
%
% Requires Matlab 2022 or newer!
%
function DARSim(varargin)
startup_darsim

% Setup multicore for solver
% Written by Mathias Erdtmann
if length(varargin) > 2
    ncores = varargin{3};
    if ncores == 1
        delete(gcp('nocreate'));
    elseif ncores > 1 && isempty(gcp('nocreate'))
        myCluster = parcluster('Processes'); parpool(myCluster, ncores);
    elseif ncores > 1 && ~isempty(gcp('nocreate')) && ncores ~= gcp('nocreate').NumWorkers
        delete(gcp('nocreate'));
        myCluster = parcluster('Processes'); parpool(myCluster, ncores);
     end
     varargin = {varargin{1}, varargin{2}};
end

if ~isempty(varargin)
    [folder, file] = define_input_path(varargin);
    start_diary(folder);
    simulation_settings = create_darsim_settings(folder,file);
    simulator = create_simulator(simulation_settings);
    simulator.Run();
    diary off
    fclose('all');
else
	% If input is not provided, DARSim initialize user interface.
    darsim_menu();
end
end


%Mengjie was here

