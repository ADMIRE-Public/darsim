%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DARSim
% Author: Artur Castiel Reis de Souza
% TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [simulation_settings , blueprint] = ...
    create_darsim_settings(folder, file)
% this routine reads the simulation input file (folder + file)
% creates a folder containing all the information and adds default and
% missing parameters not included. this folder is split into blueprints
% that consist each in specific instructions that are used to create a 
% simulation. After that blueprints are translated into the simulation
% settings, the simulation factory can use it to create the simulation.
input_folder = reader_darsim3(folder, file);
blueprint = blueprint_factory(input_folder);
simulation_settings = darsim2_settings_factory(blueprint);
end
