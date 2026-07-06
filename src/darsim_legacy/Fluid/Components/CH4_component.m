% CH4 component
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Amin MB
%TU Delft
%Created: 22 July 2025
%Last modified: 22 July 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef CH4_component < component
    properties
        rho_n_phi % Apparent density, Used in calculation of liquid phase density, Check SPE11 benchmark
        X % Constants fitting parameter for correct density in volume shift method
    end
    methods
        function obj = CH4_component()
            obj.rho = 0.66816; % Set specific default value for rho in this subclass
            obj.rho_n_phi = 1e10; % We assume that the apparent density of H2 is Inf
                                  % This means that the dissolved H2
                                  % doesn't add to the volume of the
                                  % liquid. Assumed 1e10 instead of Inf
                                  % for stability reasons of drho of liquid
                                  % phase.
            obj.Pc = 45.992;  % Critical pressure [bar]
            obj.Tc = 190.564; % Critical temperature [K]
            obj.w = 0.01142;  % Acentric factor
            obj.Mw = 16.0425; % Molecular weight [g/mol]
            obj.X = 6.404378; % Constants fitting parameter for correct density calculation for PR EOS.
            obj.Vc = 99.3; % the molar volume in the critical point. cm3/mol
        end

        
    end    
end