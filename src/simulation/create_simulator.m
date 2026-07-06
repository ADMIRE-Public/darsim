%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel;
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function simulator = create_simulator(simulation_settings)
% create_simulator creates the simulation object based on the simulation_settings.
% intialize the simulation factory object.
factory = simulation_factory(simulation_settings);
% simulation factory object produces the simulator object.
simulator = factory.BuildSimulation();
% start up the simulator object
simulator.prepare_simulation();
end
