% Compositional gas phase class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini, Yuhang Wang
%TU Delft
%Created: 28 September 2016
%Last modified: 06 May 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Compositional_gas_phase < phase
    properties
        % Pb = 2.5e7;
        % Bg = 0.0025;
        Pst = 1.01325e5;
        Tst = 293.15;
        T
        % alpha

        % gridded interpolant object for data
        % zcf is the compressibility factor of gas phase
        zcf_p_table; % attn. zcf: compressibility factor; z: mole(mass) fraction of a component;
        dzcf_p_table;
    end
    methods
        function obj = Compositional_gas_phase(mu, T) % N is number of matrix cells
            obj.mu = mu;
            obj.T = T;
        end
        function rho = ComputeDensity(obj, p, Components, ~, ~)            
            zcf = obj.zcf_p_table(p);
            Bg = obj.Pst/obj.Tst*zcf*obj.T./p;
            rho = Components(1).rho./Bg;
        end
        function drho = ComputeDrhoDp(obj, p, Components, ~, ~)
            zcf = obj.zcf_p_table(p);
            dzcf = obj.dzcf_p_table(p);
            Bg = obj.Pst/obj.Tst*zcf*obj.T./p;
            drho = - Components(1).rho./(Bg.^2).*obj.Pst/obj.Tst*obj.T.*(dzcf.*p - zcf)./(p.^2);
        end
        function drho = ComputeDrhoDz(obj, p, Components, ~, ~, ~)
            drho = zeros(length(p),1);
        end
        function mu = ComputeViscosity(obj, p, Components, Thermo, ~, T_in)
            if nargin < 6 || isempty(T_in)
                Tuse = obj.T;
            else
                Tuse = T_in;
            end

            if ~isempty(obj.mu_pT_table)
                mu = obj.mu_pT_table(p, Tuse);
                return;
            end
            
            mu = obj.mu; 
        end
    end    
end