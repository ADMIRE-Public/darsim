% CO2 component
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 28 April 2025
%Last modified: 22 July 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef CO2_component < component
    properties
        rho_n_phi % Apparent density of CO2, Used in calculation of liquid phase density, Check SPE11 benchmark
    end
    methods
        function obj = CO2_component(T)
            obj.rho = 1.98; % Set specific default value for rho in this subclass
            T_C = T - 273.15;
            V_phi = 1e-6 * (37.51 - 9.585e-2 * T_C + 8.74e-4 * T_C^2 - 5.044e-7 * T_C^3);
            obj.rho_n_phi = 44.01e-3 / V_phi;  
            obj.Pc = 73.773; % Critical pressure [bar]
            obj.Tc = 304.1282; % Critical temperature [K]
            obj.w = 0.22394;  % Acentric factor
            obj.Mw = 44.009; % Molecular weight [g/mol]
            obj.Vc = 94.1168; % the molar volume in the critical point. cm3/mol
        end

        
    end    
end