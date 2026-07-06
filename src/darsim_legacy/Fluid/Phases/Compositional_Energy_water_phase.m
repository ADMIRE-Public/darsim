% Compositional Energy water phase class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Mengjie Zhao
%TU Delft
%Created: 09 Feb 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Compositional_Energy_water_phase < phase
    properties
        cond
    end
    methods
        function obj = Compositional_Energy_water_phase(cond)
            obj.cond = cond;
        end

        function rho = ComputeDensity(obj, p, T, status, Thermo)
            Tc = T - 273.15;
            x = status.Properties('x_1ph2').Value;
            rho_brine = Thermo.rhoAq_pT_table(p,T);

            V_app   = (37.51 - 9.585e-2 .* Tc + 8.740e-4 .* Tc .^ 2 - 5.044e-7.* Tc.^3) * 1e-6; % m^3/mol
            MW_CO2  = 44.01;                                                                    % g/mol
            rho_CO2 = MW_CO2 * 1e-3 ./ V_app;
            rho = 1.0 ./ ( (1.0 - x) ./ rho_brine + x ./ rho_CO2 );
        end

        function drho = ComputeDrhoDp(obj, p, T, status, Thermo, dxdp)
            Tc = T - 273.15;
            x  = status.Properties('x_1ph2').Value;
            rho_b    = Thermo.rhoAq_pT_table(p,T);
            drhob_dp = Thermo.drhoAq_dP_pT_table(p,T);

            V_app    = (37.51 - 9.585e-2 .* Tc + 8.740e-4 .* Tc .^ 2 - 5.044e-7.* Tc.^3) * 1e-6; % m^3/mol
            MW_CO2   = 44.01;                                                                    % g/mol
            rho_CO2  = MW_CO2 * 1e-3 ./ V_app;

            den      = (1 - x) ./ rho_b +  x ./ rho_CO2;
            rho      = 1 ./ den;

            drho     =  rho .^ 2 .* ((1-x)./rho_b.^2 .* drhob_dp + dxdp(:,2) .* (1./rho_b - 1./rho_CO2));
        end

        function drho = ComputeDrhoDT(obj, p, T, status, Thermo, dxdT)                        
            Tc = T - 273.15;
            x  = status.Properties('x_1ph2').Value;

            % brine
            rho_b   = Thermo.rhoAq_pT_table(p,T);
            drhob_dT  = Thermo.drhoAq_dT_pT_table(p,T);

            % CO2 (apparent molar volume correlation)
            a0 = 37.51;a1 = -9.585e-2;a2 =  8.740e-4;a3 = -5.044e-7;

            V_app    = (a0 + a1.*Tc + a2.*Tc.^2 + a3.*Tc.^3) * 1e-6;  % m^3/mol
            dVdT     = (a1 + 2*a2.*Tc + 3*a3.*Tc.^2) * 1e-6;          % dV_app/dT

            MW_CO2   = 44.01e-3;                                      % kg/mol
            rho_CO2  = MW_CO2 ./ V_app;
            drhoC_dT = -MW_CO2 ./ V_app.^2 .* dVdT;                   % drho_CO2/dT

            % mixture
            den   = (1 - x) ./ rho_b + x ./ rho_CO2;
            rho = 1 ./ den;

            drho = rho .^ 2 .* ( (1-x)./rho_b.^2 .* drhob_dT + dxdT(:,2) .* (1./rho_b - 1./rho_CO2) + x./rho_CO2.^2 .* drhoC_dT );
        end

        function drho = ComputeDrhoDz(obj, p, T, status, Thermo, dxdz)
            Tc = T - 273.15;
            x      = status.Properties('x_1ph2').Value;
            rho_b  = Thermo.rhoAq_pT_table(p,T);

            V_app    = (37.51 - 9.585e-2 .* Tc + 8.740e-4 .* Tc .^ 2 - 5.044e-7.* Tc.^3) * 1e-6; % m^3/mol
            MW_CO2   = 44.01;                                                                    % g/mol
            rho_CO2  = MW_CO2 * 1e-3 ./ V_app;

            den      = (1 - x) ./ rho_b +  x ./ rho_CO2;
            rho      = 1 ./ den;

            drho = rho .^ 2 .* dxdz(:,2) .* (1./rho_b - 1./rho_CO2);
        end

        function h = ComputeEnthalpy(obj, p, T, status, Thermo)
            x       = status.Properties('x_1ph2').Value;
            h_dis   = Thermo.hCO2_pT_table(p, T);
            h_brine = Thermo.hb_pT_table(p, T);

            h = (1-x) .* h_brine + x .* h_dis;
        end

        function dh = ComputeDhDp(obj, p, T, status, Thermo, dxdp)

            x       = status.Properties('x_1ph2').Value;
            h_dis   = Thermo.hCO2_pT_table(p, T);
            h_brine = Thermo.hb_pT_table(p, T);
            
            dh = - h_brine .* dxdp(:,2) + (1-x) .*  Thermo.dhb_dP_pT_table(p,T) + x .* Thermo.dhCO2_dP_pT_table(p,T) + dxdp(:,2) .* h_dis;
            
        end

        function dh = ComputeDhDT(obj, p, T, status, Thermo, dxdT)

            x       = status.Properties('x_1ph2').Value;
            h_dis   = Thermo.hCO2_pT_table(p, T);
            h_brine = Thermo.hb_pT_table(p, T);
            
            dh = - h_brine .* dxdT(:,2) + (1-x) .*  Thermo.dhb_dT_pT_table(p,T) + x .* Thermo.dhCO2_dT_pT_table(p,T) + dxdT(:,2) .* h_dis;
        end

        function dh = ComputeDhDz(obj, p, T, status, Thermo, dxdz)

            h_dis  = Thermo.hCO2_pT_table(p, T);
            h_brine = Thermo.hb_pT_table(p, T);
            
            dh = - h_brine .* dxdz(:,2) + dxdz(:,2) .*  h_dis;
        end

        function mu = ComputeViscosity(obj, p, T, status, Thermo)
            x    = status.Properties('x_1ph2').Value;
            Mw_CO2 = 44.01;        % g/mol
            Mw_H2O = 18.0152;      % g/mol

            x_molar = x / Mw_CO2 ./ (x/Mw_CO2 + (1-x)./Mw_H2O);
            mu_b = Thermo.muAq_pT_table(p,T);

            mu = mu_b .* (1.0 + 4.65 .* x_molar .^ 1.0134);
        end

        function dmu = ComputeDmuDp(obj, p, T, status, Thermo, dxdp)
            x    = status.Properties('x_1ph2').Value;
            dmu = zeros(length(x),1);
            Mw_CO2 = 44.01;        % g/mol
            Mw_H2O = 18.0152;      % g/mol

            x_molar = x / Mw_CO2 ./ (x/Mw_CO2 + (1-x)./Mw_H2O);

            dxm_dx = 1 / (Mw_CO2 * Mw_H2O) ./ (x/Mw_CO2 + (1-x)./Mw_H2O);
            dxm_dp = dxm_dx .* dxdp(:,2);

            mu_b = Thermo.muAq_pT_table(p,T);
            dmu_b = Thermo.dmuAq_dP_pT_table(p,T);

            dmu(x_molar > 0,:) =  dmu_b(x_molar>0) .* (1.0 + 4.65 .* x_molar(x_molar>0) .^ 1.0134) + mu_b(x_molar>0) .* 4.65 .*  1.0134 .* x_molar(x_molar>0) .^ 0.0134 .* dxm_dp(x_molar>0,:);
            dmu(x_molar ==0,:) =  dmu_b(x_molar==0);
        end

        function dmu = ComputeDmuDT(obj, p, T, status, Thermo, dxdT)
            x    = status.Properties('x_1ph2').Value;
            dmu = zeros(length(x),1);
            Mw_CO2 = 44.01;        % g/mol
            Mw_H2O = 18.0152;      % g/mol

            x_molar = x / Mw_CO2 ./ (x/Mw_CO2 + (1-x)./Mw_H2O);
            dxm_dx = 1 / (Mw_CO2 * Mw_H2O) ./ (x/Mw_CO2 + (1-x)./Mw_H2O);
            dxm_dT = dxm_dx .* dxdT(:,2);

            mu_b = Thermo.muAq_pT_table(p,T);
            dmu_b = Thermo.dmuAq_dT_pT_table(p,T);

            dmu(x_molar > 0) =  dmu_b(x_molar>0) .* (1.0 + 4.65 .* x_molar(x_molar>0) .^ 1.0134) + mu_b(x_molar>0) .* 4.65 .*  1.0134 .* x_molar(x_molar>0) .^ 0.0134 .* dxm_dT(x_molar>0,:);
            dmu(x_molar ==0) =  dmu_b(x_molar==0);
        end

        function dmu = ComputeDmuDz(obj, p, T, status, Thermo, dxdz)
            x    = status.Properties('x_1ph2').Value;
            Mw_CO2 = 44.01;        % g/mol
            Mw_H2O = 18.0152;      % g/mol

            x_molar = x / Mw_CO2 ./ (x/Mw_CO2 + (1-x)./Mw_H2O);
            dxm_dx = 1 / (Mw_CO2 * Mw_H2O) ./ (x/Mw_CO2 + (1-x)./Mw_H2O);
            dxm_dz = dxm_dx .* dxdz(:,2);

            mu_b = Thermo.muAq_pT_table(p,T);

            dmu = zeros(size(x));
            dmu(x_molar > 0) = mu_b(x_molar > 0) .* 4.65 .*  1.0134 .* x_molar(x_molar>0) .^ 0.0134 .* dxm_dz(x_molar>0,:);
        end

    end    
end