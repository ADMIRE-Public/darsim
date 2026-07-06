% BO gas phase class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini, Yuhang Wang
%TU Delft
%Created: 28 September 2016
%Last modified: 06 May 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef BO_gas_phase < phase
    properties
        % Pb = 2.5e7;
        % Bg = 0.0025;
        Pst = 1.01325e5;
        Tst = 293.15;
        T
        % alpha

        % gridded interpolant object for data
        % Rs and zcf are pressure dependent
        % Rs is the dissolution of gas component in liquid phase and is
        % implemented in BO_oil_phase
        % zcf is the compressibility factor of gas phase and is implemented
        % in BO_gas_phase
        zcf_p_table; % attn. zcf: compressibility factor; z: mole(mass) fraction of a component;
        dzcf_p_table;
    end
    methods
        function obj = BO_gas_phase(mu, T) % N is number of matrix cells
            obj.mu = mu;
            obj.T = T;
            % obj.sr = 0; % it's assigned from input file
            % obj.alpha = log(obj.Bg)/(obj.Pb/obj.Pst - 1)*1/obj.Pst; 
            
            % allocate memory
            % obj.zcf = zeros(N,1);
            % obj.dzcf = zeros(N,1);
            % obj.Bg = zeros(N,1);
        end
        function rho = ComputeDensity(obj, p, Components, rs)
            % N = length(p);
            % Bg_r = 0.25 * ones(length(Bg_r),1);
            % Bg_r = exp(obj.alpha .* (p - obj.Pst));
            % Bg_r = 0.5*(p./obj.Pst).^-1;
            
            % Bg_r = 0.25;
            % rho = zeros(length(p),1)+Components(1).rho./Bg_r;
            
            zcf = obj.zcf_p_table(p);
            Bg = obj.Pst/obj.Tst*zcf*obj.T./p;
            rho = Components(1).rho./Bg;
        end
        function drho = ComputeDrhoDp(obj, p, Components)
            % Bg_r = exp(obj.alpha .* (p - obj.Pst));
            % Num = -obj.alpha * Bg_r;
            % drho = zeros(length(Bg_r), 1);
            % Bg_r = 0.5*(p./obj.Pst).^-1;
            % Num = 0.5*obj.Pst * p.^-2;
            % Den = (Bg_r).^2;
            % drho = Num ./ Den;
            % 
            % drho = 0;     
            
            zcf = obj.zcf_p_table(p);
            dzcf = obj.dzcf_p_table(p);
            Bg = obj.Pst/obj.Tst*zcf*obj.T./p;
            
            drho = - Components(1).rho./(Bg.^2).*obj.Pst/obj.Tst*obj.T.*(dzcf.*p - zcf)./(p.^2);
        end
        function [Rs, dRs] = ComputeRs(obj, p)
            Rs = 0;
            dRs = 0;
        end
        function [Rs, dRs] = RsOfUnderSaturatedPhase(obj, z, Components, Rs, dRs, SinglePhase)
            
        end
    end    
end