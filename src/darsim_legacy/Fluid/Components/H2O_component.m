% H2O component
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 28 April 2025
%Last modified: 22 July 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef H2O_component < component
    properties
        % rho = 1060; % It's the stock tank density in kg/m^3
    end
    methods
        function obj = H2O_component()
            obj.rho = 1060; % Set specific default value for rho in this subclass
            obj.Pc = 220.640; % Critical pressure [bar]
            obj.Tc = 647.096; % Critical temperature [K]
            obj.w = 0.3443;   % Acentric factor
            obj.Mw = 18.0152; % Molecular weight [g/mol]
            obj.Vc = 55.9478; % the molar volume in the critical point. cm3/mol
        end

       
    end    
end