% H2 component
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 28 April 2025
%Last modified: 06 May 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef H2_component < component
    properties
        rho_n_phi % Apparent density, Used in calculation of liquid phase density, Check SPE11 benchmark
        X % Constants fitting parameter for correct density in volume shift method
    end
    methods
        function obj = H2_component()
            obj.rho = 0.083725; % Set specific default value for rho in this subclass
            obj.rho_n_phi = 1e10; % We assume that the apparent density of H2 is Inf
                                  % This means that the dissolved H2
                                  % doesn't add to the volume of the
                                  % liquid. Assumed 1e10 instead of Inf
                                  % for stability reasons of drho of liquid
                                  % phase.
            obj.Pc = 12.964; % Critical pressure [bar]
            obj.Tc = 33.145; % Critical temperature [K]
            obj.w = -0.219;  % Acentric factor
            obj.Mw = 2.0157; % Molecular weight [g/mol]
            obj.X = 4.8070; % Constants fitting parameter for correct density in volume shift method
            obj.Vc = 64.4828; % the molar volume in the critical point. cm3/mol
        end

        
    end    
end