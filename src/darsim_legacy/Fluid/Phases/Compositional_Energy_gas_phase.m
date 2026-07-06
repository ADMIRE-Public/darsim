% Compositional Energy gas phase class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini, Yuhang Wang
%TU Delft
%Created: 28 September 2016
%Last modified: 06 May 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Compositional_Energy_gas_phase < phase
    properties
        % Pb = 2.5e7;
        % Bg = 0.0025;
        Pst = 1.01325e5;
        Tst = 273.15;
        T
        cond % Heat conductivity

        % alpha
    end
    methods
        function obj = Compositional_Energy_gas_phase(cond) % N is number of matrix cells
            obj.cond = cond;
        end

        function rho = ComputeDensity(obj, p, T, status, Thermo)
            SinglePhase = status.Properties('SinglePhase').Value;
            z = status.Properties('z_1').Value;
            rho = Thermo.rhoV_pT_table(p,T);
            rho(SinglePhase == 1) = Thermo.rhoV_pTz_table(p(SinglePhase==1),T(SinglePhase==1),z(SinglePhase==1)); 
        end

        function drho = ComputeDrhoDp(obj, p, T, status, Thermo)
            SinglePhase = status.Properties('SinglePhase').Value;
            z = status.Properties('z_1').Value;
            drho = Thermo.drhoV_dP_pT_table(p,T);
            drho(SinglePhase == 1) = Thermo.drhoV_dP_pTz_table(p(SinglePhase==1),T(SinglePhase==1),z(SinglePhase==1)); 
        end

        function drho = ComputeDrhoDT(obj, p, T, status, Thermo)
            SinglePhase = status.Properties('SinglePhase').Value;
            z = status.Properties('z_1').Value;
            drho = Thermo.drhoV_dT_pT_table(p,T);
            drho(SinglePhase == 1) = Thermo.drhoV_dT_pTz_table(p(SinglePhase==1),T(SinglePhase==1),z(SinglePhase==1)); 
        end
        
        function drho = ComputeDrhoDz(obj, p, T, status, Thermo)
            SinglePhase = status.Properties('SinglePhase').Value;
            z = status.Properties('z_1').Value;
            drho = zeros(size(p));
            drho(SinglePhase == 1) = Thermo.drhoV_dz_pTz_table(p(SinglePhase==1),T(SinglePhase==1),z(SinglePhase==1));
        end
        
        function h = ComputeEnthalpy(obj, p, T, status, Thermo)
            SinglePhase = status.Properties('SinglePhase').Value;
            z = status.Properties('z_1').Value;
            h = Thermo.hV_pT_table(p,T);
            h(SinglePhase == 1) = Thermo.hV_pTz_table(p(SinglePhase==1),T(SinglePhase==1),z(SinglePhase==1)); 
        end

        function dh = ComputeDhDp(obj, p, T, status, Thermo)
            SinglePhase = status.Properties('SinglePhase').Value;
            z = status.Properties('z_1').Value;
            dh = Thermo.dhV_dP_pT_table(p,T);
            dh(SinglePhase == 1) = Thermo.dhV_dP_pTz_table(p(SinglePhase==1),T(SinglePhase==1),z(SinglePhase==1)); 
        end

        function dh = ComputeDhDT(obj, p, T, status, Thermo)
            SinglePhase = status.Properties('SinglePhase').Value;
            z = status.Properties('z_1').Value;
            dh = Thermo.dhV_dT_pT_table(p,T);
            dh(SinglePhase == 1) = Thermo.dhV_dT_pTz_table(p(SinglePhase==1),T(SinglePhase==1),z(SinglePhase==1)); 
        end

        function dh = ComputeDhDz(obj, p, T, status, Thermo)
            SinglePhase = status.Properties('SinglePhase').Value;
            z = status.Properties('z_1').Value;
            dh = zeros(size(p));
            dh(SinglePhase == 1) = Thermo.dhV_dz_pTz_table(p(SinglePhase==1),T(SinglePhase==1),z(SinglePhase==1)); 
        end

        function mu = ComputeViscosity(obj, p, T, status, Thermo)
            mu = Thermo.muV_pT_table(p,T);
        end

        function dmu = ComputeDmuDp(obj, p, T, Thermo)
            dmu = Thermo.dmuV_dP_pT_table(p,T);
        end

        function dmu = ComputeDmuDz(obj, p, T, Thermo)
            dmu = zeros(size(p));
        end

        function dmu = ComputeDmuDT(obj, p, T, Thermo)
            dmu = Thermo.dmuV_dT_pT_table(p,T);
        end

    end    
end