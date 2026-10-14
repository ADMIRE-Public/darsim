% Thermodynamics of CO2-Brine for Compositional Thermal Fluid Model class
%--------------------------------------------------------------------------
% DARSim 2 Reservoir Simulator
% Author: Mengjie Zhao
% TU Delft
% Created: 02 Feb 2026
%--------------------------------------------------------------------------
classdef Thermo_CO2_Brine_Compositional_Energy < matlab.mixin.Heterogeneous & handle
    % This MATLAB class calculates thermodynamic properties and phase
    % equilibria for a system containing CO2 and brine
    % (H2O with salinity). The code is based on the Peng-Robinson Equation of
    % State (PR EOS) and follows the methodology described by
    % Ziabakhsh-Ganji and Kooi (2012).
    %
    % References:
    %
    % 1. Ziabakhsh-Ganji, Z., & Kooi, H. (2012). An Equation of State for
    % thermodynamic equilibrium of gas mixtures and brines to allow simulation
    % of the effects of impurities in subsurface CO2 storage. International
    % Journal of Greenhouse Gas Control, 11, S21-S34.
    %
    % 2. Ziabakhsh-Ganji et al., (2014). Sensitivity of Joule–Thomson cooling
    % to impure CO2 injection in depleted gas reservoirs. Applied energy, 113,
    % 434-451.
    %
    % 3. Yongfan Guo et al., (2019). An enthalpy model of CO2-CH4-H2S-N2-brine
    % systems applied in simulation of non-isothermal multiphase and
    % multicomponent flow with 	high pressure, temperature and salinity".
    % Journal of CO2 Utilization, pp. 85–97.

    properties

        T_ref = 298.15;
        % gridded interpolant object for data

        % 2D tables
        muV_pT_table
        muAq_pT_table
        dmuV_dP_pT_table
        dmuV_dT_pT_table
        dmuAq_dP_pT_table
        dmuAq_dT_pT_table

        hV_pT_table
        dhV_dP_pT_table
        dhV_dT_pT_table

        hb_pT_table
        hCO2_pT_table
        dhb_dP_pT_table
        dhb_dT_pT_table
        dhCO2_dP_pT_table
        dhCO2_dT_pT_table

        rhoV_pT_table
        rhoAq_pT_table
        drhoV_dP_pT_table
        drhoV_dT_pT_table
        drhoAq_dP_pT_table
        drhoAq_dT_pT_table

        rhoV_pTz_table
        drhoV_dP_pTz_table
        drhoV_dT_pTz_table
        drhoV_dz_pTz_table

        hV_pTz_table          % Single phase vapor
        dhV_dP_pTz_table
        dhV_dT_pTz_table
        dhV_dz_pTz_table

        X_comp1_pT_table
        Y_comp1_pT_table
        X_comp2_pT_table
        Y_comp2_pT_table

        dX_comp1_dP_pT_table
        dX_comp1_dT_pT_table
        dY_comp1_dP_pT_table
        dY_comp1_dT_pT_table

        dX_comp2_dP_pT_table
        dX_comp2_dT_pT_table
        dY_comp2_dP_pT_table
        dY_comp2_dT_pT_table
    end
    methods

        function obj = Thermo_CO2_Brine_Compositional_Energy(Salinity)
            % Salinity (ppm:mg/kg)
            obj.ComputeLookupTable_2D(Salinity);
        end

        function ComputeLookupTable_2D(obj, Salinity)
            % Vapor:   PR EOS (with Mathias1989 volume shift) + enthalpy
            % Aqueous: Ziabakhsh2012 (phi_aq) + enthalpy
            % rhoAq:   Garcia2001 (Spivey2004 + dissolved CO2 correction)
            % muAq :   Islam2012  (MaoDuan2009 + CO2 correction)
            % muV  :   Fenghour1998 (CO2-rich gas viscosity)
            % Enthalpy: cp_R + PR departure + Guo2019

            P_grid = linspace(1e5, 4e7, 400)';          % Pa
            T_grid = linspace(293.15, 373.15, 81)';     % K
            z_grid = linspace(0, 1, 41)';
            nP = numel(P_grid);
            nT = numel(T_grid);
            nZ = numel(z_grid);

            % Salinity
            m_Na = Salinity/1e3/58.4427; % mol/kg
            m_Cl = Salinity/1e3/58.4427; % mol/kg
            m_NaCl = m_Na;

            % constants
            Mw_CO2 = 44.01;        % g/mol
            Mw_H2O = 18.0152;      % g/mol
            Mw_NaCl = 58.4427;     % g/mol
            Nw = 55.508;           % mol H2O per kg-water

            % allocate tables
            X_CO2 = nan(nP,nT);  Y_CO2 = nan(nP,nT);
            X_H2O = nan(nP,nT);  Y_H2O = nan(nP,nT);
            rhoV  = nan(nP,nT);  rhoAq = nan(nP,nT);
            muV   = nan(nP,nT);  muAq  = nan(nP,nT);
            hV    = nan(nP,nT);  h_b   = nan(nP,nT);
            h_dis = zeros(nP,nT);
            rhoV_3D = nan(nP,nT,nZ);
            hV_3D = nan(nP,nT,nZ);
            % iteration controls
            max_iter = 80;

            for j = 1 : nT
                T = T_grid(j);

                for i = 1 : nP
                    P = P_grid(i);

                    % aq fugacity coeffs at this (P,T)
                    [phi_aq_CO2, phi_aq_H2O] = obj.zia_phi_aq_CO2_H2O(P, T, m_Na, m_Cl);

                    y2 = 1e-4;

                    for it = 1 : max_iter
                        y1 = 1 - y2;
                        [~, phi_g_CO2, phi_g_H2O, ~] = obj.pr_phi_binary(P, T, y1, y2);

                        x_CO2 = (y1 * phi_g_CO2) / phi_aq_CO2;

                        x_H2O = 1 - x_CO2;

                        y2_new = x_H2O * (phi_aq_H2O / phi_g_H2O);

                        if abs(y2_new - y2) < 1e-12
                            y2 = y2_new;
                            break;
                        end
                        y2 = y2_new;
                    end

                    y_H2O = y2;
                    y_CO2 = 1 - y2;
                    x_CO2 = (y1 * phi_g_CO2) / phi_aq_CO2;
                    x_H2O = 1 - x_CO2;

                    % ---- store GAS mass fractions ----
                    Y_CO2(i,j) = (y_CO2*Mw_CO2) / (y_CO2*Mw_CO2 + y_H2O*Mw_H2O);
                    Y_H2O(i,j) = 1 - Y_CO2(i,j);

                    % ---- compute aqueous dissolved CO2 molality from x_CO2/x_H2O ----
                    mCO2 = x_CO2 / x_H2O * Nw;     % mol/kg-water
                    mass_brine = 1.0 + m_NaCl * Mw_NaCl * 1e-3;      % kg
                    mass_aq    = mass_brine + mCO2 * Mw_CO2 * 1e-3;  % kg

                    X_CO2(i,j) = mCO2 * Mw_CO2 * 1e-3 / mass_aq;                   % mass fraction in aqueous
                    X_H2O(i,j) = mass_brine / mass_aq;                      % mass fraction in aqueous

                    % densities
                    [rho_g, ~, ~, h_g] = obj.pr_phi_binary(P, T, y_CO2, y_H2O);
                    rhoV(i,j)  = rho_g;
                    rhoAq(i,j) = obj.garcia2001(P, T, m_NaCl);

                    % viscosities
                    muV(i,j)  = obj.mu_Fenghour1998(T, rho_g);
                    muAq(i,j) = obj.mu_Islam2012(P, T, m_NaCl);

                    % enthalpies
                    hV(i,j) = h_g;
                    [~, ~, h_brine, h_d] = obj.zia_phi_aq_CO2_H2O(P, T, m_Na, m_Cl, X_CO2(i,j), X_H2O(i,j));
                    h_b(i,j) = h_brine;
                    h_dis(i,j) = h_d;
                end
            end


            for k = 1 : nZ
                zCO2 = z_grid(k);
                zH2O = 1 - zCO2;

                denom = zCO2/Mw_CO2 + zH2O/Mw_H2O;
                yCO2  = (zCO2/Mw_CO2) / denom;      % mole fraction
                yH2O  = 1 - yCO2;

                for j = 1:nT
                    T = T_grid(j);

                    for i = 1:nP
                        P = P_grid(i);
                        [rho_g, ~, ~, h_C] = obj.pr_phi_binary(P, T, yCO2, yH2O);
                        rhoV_3D(i,j,k) = rho_g;
                        hV_3D(i,j,k) = h_C;
                    end
                end
            end

            d_dP = @(A) ...
                [ (A(2,:)-A(1,:)) ./ (P_grid(2)-P_grid(1)); ...
                (A(3:end,:)-A(1:end-2,:)) ./ (P_grid(3:end)-P_grid(1:end-2)); ...
                (A(end,:)-A(end-1,:)) ./ (P_grid(end)-P_grid(end-1)) ];

            d_dT = @(A) ...
                [ (A(:,2)-A(:,1)) ./ (T_grid(2)-T_grid(1)), ...
                (A(:,3:end)-A(:,1:end-2)) ./ (T_grid(3:end)-T_grid(1:end-2)).', ...
                (A(:,end)-A(:,end-1)) ./ (T_grid(end)-T_grid(end-1)) ];

            dmuV_dP = d_dP(muV);   dmuV_dT = d_dT(muV);
            dmuA_dP = d_dP(muAq);  dmuA_dT = d_dT(muAq);

            dhV_dP  = d_dP(hV);    dhV_dT  = d_dT(hV);
            dhb_dP  = d_dP(h_b);   dhb_dT  = d_dT(h_b);
            dhCO2_dP= d_dP(h_dis); dhCO2_dT= d_dT(h_dis);

            drV_dP  = d_dP(rhoV);  drV_dT  = d_dT(rhoV);
            drA_dP  = d_dP(rhoAq); drA_dT  = d_dT(rhoAq);

            dX1_dP  = d_dP(X_CO2); dX1_dT  = d_dT(X_CO2);
            dY1_dP  = d_dP(Y_CO2); dY1_dT  = d_dT(Y_CO2);

            dX2_dP  = d_dP(X_H2O); dX2_dT  = d_dT(X_H2O);
            dY2_dP  = d_dP(Y_H2O); dY2_dT  = d_dT(Y_H2O);

            drhoV_dP = obj.fd3D_axis(rhoV_3D, P_grid, 1);
            drhoV_dT = obj.fd3D_axis(rhoV_3D, T_grid, 2);
            drhoV_dz = obj.fd3D_axis(rhoV_3D, z_grid, 3);

            dhV_3D_dP = obj.fd3D_axis(hV_3D, P_grid, 1);
            dhV_3D_dT = obj.fd3D_axis(hV_3D, T_grid, 2);
            dhV_3D_dz = obj.fd3D_axis(hV_3D, z_grid, 3);

            obj.muV_pT_table      = griddedInterpolant({P_grid, T_grid}, muV,   'linear', 'nearest');
            obj.muAq_pT_table     = griddedInterpolant({P_grid, T_grid}, muAq,  'linear', 'nearest');
            obj.dmuV_dP_pT_table  = griddedInterpolant({P_grid, T_grid}, dmuV_dP,'linear', 'nearest');
            obj.dmuV_dT_pT_table  = griddedInterpolant({P_grid, T_grid}, dmuV_dT,'linear', 'nearest');
            obj.dmuAq_dP_pT_table = griddedInterpolant({P_grid, T_grid}, dmuA_dP,'linear', 'nearest');
            obj.dmuAq_dT_pT_table = griddedInterpolant({P_grid, T_grid}, dmuA_dT,'linear', 'nearest');

            obj.hV_pT_table      = griddedInterpolant({P_grid, T_grid}, hV,    'linear', 'nearest');
            obj.dhV_dP_pT_table  = griddedInterpolant({P_grid, T_grid}, dhV_dP,'linear', 'nearest');
            obj.dhV_dT_pT_table  = griddedInterpolant({P_grid, T_grid}, dhV_dT,'linear', 'nearest');

            obj.hb_pT_table     = griddedInterpolant({P_grid, T_grid}, h_b,   'linear', 'nearest');
            obj.dhb_dP_pT_table = griddedInterpolant({P_grid, T_grid}, dhb_dP,'linear', 'nearest');
            obj.dhb_dT_pT_table = griddedInterpolant({P_grid, T_grid}, dhb_dT,'linear', 'nearest');
            obj.hCO2_pT_table     = griddedInterpolant({P_grid, T_grid}, h_dis,   'linear', 'nearest');
            obj.dhCO2_dP_pT_table = griddedInterpolant({P_grid, T_grid}, dhCO2_dP,'linear', 'nearest');
            obj.dhCO2_dT_pT_table = griddedInterpolant({P_grid, T_grid}, dhCO2_dT,'linear', 'nearest');

            obj.rhoV_pT_table      = griddedInterpolant({P_grid, T_grid}, rhoV,   'linear', 'nearest');
            obj.rhoAq_pT_table     = griddedInterpolant({P_grid, T_grid}, rhoAq,  'linear', 'nearest');
            obj.drhoV_dP_pT_table  = griddedInterpolant({P_grid, T_grid}, drV_dP, 'linear', 'nearest');
            obj.drhoV_dT_pT_table  = griddedInterpolant({P_grid, T_grid}, drV_dT, 'linear', 'nearest');
            obj.drhoAq_dP_pT_table = griddedInterpolant({P_grid, T_grid}, drA_dP, 'linear', 'nearest');
            obj.drhoAq_dT_pT_table = griddedInterpolant({P_grid, T_grid}, drA_dT, 'linear', 'nearest');

            obj.rhoV_pTz_table     = griddedInterpolant({P_grid, T_grid, z_grid}, rhoV_3D, 'linear', 'nearest');
            obj.drhoV_dP_pTz_table = griddedInterpolant({P_grid, T_grid, z_grid}, drhoV_dP, 'linear','nearest');
            obj.drhoV_dT_pTz_table = griddedInterpolant({P_grid, T_grid, z_grid}, drhoV_dT, 'linear','nearest');
            obj.drhoV_dz_pTz_table = griddedInterpolant({P_grid, T_grid, z_grid}, drhoV_dz, 'linear','nearest');

            obj.hV_pTz_table     = griddedInterpolant({P_grid, T_grid, z_grid}, hV_3D, 'linear', 'nearest');
            obj.dhV_dP_pTz_table = griddedInterpolant({P_grid, T_grid, z_grid}, dhV_3D_dP, 'linear','nearest');
            obj.dhV_dT_pTz_table = griddedInterpolant({P_grid, T_grid, z_grid}, dhV_3D_dT, 'linear','nearest');
            obj.dhV_dz_pTz_table = griddedInterpolant({P_grid, T_grid, z_grid}, dhV_3D_dz, 'linear','nearest');

            obj.X_comp1_pT_table = griddedInterpolant({P_grid, T_grid}, X_CO2, 'linear', 'nearest');
            obj.Y_comp1_pT_table = griddedInterpolant({P_grid, T_grid}, Y_CO2, 'linear', 'nearest');
            obj.X_comp2_pT_table = griddedInterpolant({P_grid, T_grid}, X_H2O, 'linear', 'nearest');
            obj.Y_comp2_pT_table = griddedInterpolant({P_grid, T_grid}, Y_H2O, 'linear', 'nearest');

            obj.dX_comp1_dP_pT_table = griddedInterpolant({P_grid, T_grid}, dX1_dP, 'linear', 'nearest');
            obj.dX_comp1_dT_pT_table = griddedInterpolant({P_grid, T_grid}, dX1_dT, 'linear', 'nearest');
            obj.dY_comp1_dP_pT_table = griddedInterpolant({P_grid, T_grid}, dY1_dP, 'linear', 'nearest');
            obj.dY_comp1_dT_pT_table = griddedInterpolant({P_grid, T_grid}, dY1_dT, 'linear', 'nearest');

            obj.dX_comp2_dP_pT_table = griddedInterpolant({P_grid, T_grid}, dX2_dP, 'linear', 'nearest');
            obj.dX_comp2_dT_pT_table = griddedInterpolant({P_grid, T_grid}, dX2_dT, 'linear', 'nearest');
            obj.dY_comp2_dP_pT_table = griddedInterpolant({P_grid, T_grid}, dY2_dP, 'linear', 'nearest');
            obj.dY_comp2_dT_pT_table = griddedInterpolant({P_grid, T_grid}, dY2_dT, 'linear', 'nearest');
        end

        %% PR-EOS VAPOR DENSITY + FUGACITY + ENTHALPY (mole fractions)
        function [rho, phi1, phi2, h] = pr_phi_binary(obj, P, T, y1, y2)

            R = 8.314462618;  % Pa*m^3/(mol*K)

            % Critical points [CO2,H20]
            Tc    = [304.18 , 647.14];        % K
            Pc    = [73.75  , 220.64] * 1e5;  % Pa
            w     = [0.2273  , 0.3434];       %w     = [0.239  , 0.328];         % acentric
            kij12 = 0.1896;                   %kij12 = 0.19014;
            Tr    = T ./ Tc;

            kappa = 0.37464 + 1.54226 .* w - 0.26992 * w .^ 2;
            alpha = (1 + kappa .* (1 - sqrt(Tr))) .^ 2;

            a  =  0.45724 * R ^ 2 .* Tc .^ 2 ./ Pc .* alpha;
            b  =  0.07780 * R .* Tc ./ Pc;

            da_dT = 0.45724 * R ^ 2 .* Tc .^ 2 ./ Pc .*  2 .* (1 + kappa .* (1 - sqrt(Tr))) .* ( - kappa ./ (2 .* Tc .* sqrt(Tr)));

            % mixing
            a12   = sqrt(a(1) * a(2)) * (1 - kij12);
            a_mix = y1*y1*a(1) + 2*y1*y2*a12 + y2*y2*a(2);
            b_mix = y1*b(1) + y2*b(2);

            da12_dT = a12 * 0.5 * (da_dT(1)/a(1) + da_dT(2)/a(2));
            da_mix_dT = y1 * y1 * da_dT(1) + 2 * y1 * y2 * da12_dT + y2 * y2 * da_dT(2);

            A = a_mix * P / (R ^ 2 * T^2);
            B = b_mix * P / (R * T);

            % cubic in Z for PR
            c3 = 1;
            c2 = -(1 - B);
            c1 = A - 3*B^2 - 2*B;
            c0 = -(A*B - B^2 - B^3);

            rZ = roots([c3 c2 c1 c0]);
            rZ = rZ(abs(imag(rZ)) < 1e-12);
            Z  = max(real(rZ));   % vapor root

            M1 = 44.01e-3;   % CO2 kg/mol
            M2 = 18.0152e-3;  % H2O kg/mol
            Mmix = y1*M1 + y2*M2;

            % Mathias et al. (1989) density correction
            v = Z * R * T / P;                    % uncorrected molar volume (m^3/mol)
            D = v^2 + 2*b_mix*v - b_mix^2;
            dpdv = -R*T/(v - b_mix)^2 + 2*a_mix*(v + b_mix)/D^2;

            delta = -(v^2/(R*T)) * dpdv;          % dimensionless bulk modulus

            % component volume-shift parameters s_i (Table 1, cm^3/mol -> m^3/mol)
            s_i   = [ 1.56e-6 , -3.40e-6 ];       % [CO2, H2O]  (m^3/mol)
            s_mix = y1*s_i(1) + y2*s_i(2);

            % critical molar volumes from critical densities: v_c = M / rho_c
            rho_c  = [467.6, 322.0];              % [CO2, H2O] kg/m^3
            Vc_i   = [M1/rho_c(1), M2/rho_c(2)];  % m^3/mol
            Vc_mix = y1*Vc_i(1) + y2*Vc_i(2);

            % f_c chosen to match the true critical volume (v_c^PR = 3.946*b)
            f_c = Vc_mix - (3.946*b_mix + s_mix);

            % corrected molar volume and density
            v_corr = v + s_mix + f_c * (0.41/(0.41 + delta));
            rho = Mmix / v_corr;                  % kg/m^3

            % fugacity coefficients
            sqrt2 = sqrt(2);
            ln_term = log((Z + (1+sqrt2)*B) / (Z + (1-sqrt2)*B));

            % a_i_mix terms
            a1mix = y1 * a(1)  + y2 * a12;
            a2mix = y1 * a12 + y2 * a(2);

            % common
            ln_phi1 = b(1) / b_mix * (Z - 1) - log(Z - B) ...
                - A / (2 * sqrt2 * B) * ( 2 * a1mix / a_mix - b(1) / b_mix ) * ln_term;

            ln_phi2 = b(2) / b_mix * (Z - 1) - log(Z - B) ...
                - A / (2 * sqrt2 * B) * ( 2 * a2mix / a_mix - b(2) / b_mix ) * ln_term;

            phi1 = exp(ln_phi1);
            phi2 = exp(ln_phi2);

            % ---- ideal-gas enthalpy correlation ----
            T_star = 1000; % K
            T_refe = 298.15;
            tau = T / T_star;
            tau0 = T_refe / T_star;

            % CO2 coefficients from Guo2019 implementation
            aCO2 = [ -1.8188731, 12.903022, -9.6634864, 4.2251879, -1.042164, 0.12683515, -0.49939675E-02, 2.4950242, -0.8272375, 0.15372481, ...
                -0.015861243, 8.601715E-04, -0.19222165E-04, 2.10820700252E+03];
            
            % H2O coefficients from Aleksandrov2005
            aH2O = [ 3.10409601236035E+01, -3.91422080460869E+01, 3.79695277233575E+01, -2.18374910952284E+01, 7.42251494566339E+00, -1.38178929609470E+00, ...
                1.08807067571454E-01, -1.20771176848589E+01, 3.39105078851732E+00, -5.84520979955060E-01, 5.89930846488082E-02, -3.12970001415882E-03, 6.57460740981757E-05, -1.249994285621061E+04];

            H0_CO2_mol = R * T_star * ( aCO2(1)  * tau.^1   / 1  + aCO2(2)  * tau.^2   / 2  + aCO2(3)  * tau.^3   / 3  + ...
                aCO2(4)  * tau.^4   / 4  + aCO2(5)  * tau.^5   / 5  + aCO2(6)  * tau.^6   / 6  + aCO2(7)  * tau.^7   / 7  + ...
                aCO2(8)  * log(tau)      + aCO2(9)  * tau.^(-1)/(-1)+ aCO2(10) * tau.^(-2)/(-2)+ aCO2(11) * tau.^(-3)/(-3)+ ...
                aCO2(12) * tau.^(-4)/(-4)+ aCO2(13) * tau.^(-5)/(-5) ) + R * aCO2(14) ;%+ 1.2897e+04;   % J/mol (NIST shift)

            H0_CO2_mol_ref = R * T_star * ( aCO2(1)  * tau0^1    / 1  + aCO2(2)  * tau0^2    / 2  + aCO2(3)  * tau0^3    / 3  + aCO2(4)  * tau0^4    / 4  + ...
                aCO2(5)  * tau0^5    / 5  + aCO2(6)  * tau0^6    / 6  + aCO2(7)  * tau0^7    / 7  + aCO2(8)  * log(tau0)      + ...
                aCO2(9)  * tau0^(-1) /(-1)+ aCO2(10) * tau0^(-2)/(-2)+ aCO2(11) * tau0^(-3) /(-3)+ aCO2(12) * tau0^(-4)/(-4)+ ...
                aCO2(13) * tau0^(-5) /(-5)) + R * aCO2(14);

            H0_H2O_mol = R * T_star * ( aH2O(1)  * tau.^1   / 1  + aH2O(2)  * tau.^2   / 2  + aH2O(3)  * tau.^3   / 3  + ...
                aH2O(4)  * tau.^4   / 4  + aH2O(5)  * tau.^5   / 5  + aH2O(6)  * tau.^6   / 6  + aH2O(7)  * tau.^7   / 7  + ...
                aH2O(8)  * log(tau)      + aH2O(9)  * tau.^(-1)/(-1)+ aH2O(10) * tau.^(-2)/(-2)+ aH2O(11) * tau.^(-3)/(-3)+ ...
                aH2O(12) * tau.^(-4)/(-4)+ aH2O(13) * tau.^(-5)/(-5) ) + R * aH2O(14) ;%- 8.0178e+03;     % J/mol (NIST shift)        

            H0_H2O_mol_ref = R * T_star * ( aH2O(1)  * tau0^1    / 1  + aH2O(2)  * tau0^2    / 2  + aH2O(3)  * tau0^3    / 3  + aH2O(4)  * tau0^4    / 4  + ...
                aH2O(5)  * tau0^5    / 5  + aH2O(6)  * tau0^6    / 6  + aH2O(7)  * tau0^7    / 7  + aH2O(8)  * log(tau0)      + ...
                aH2O(9)  * tau0^(-1) /(-1)+ aH2O(10) * tau0^(-2)/(-2)+ aH2O(11) * tau0^(-3) /(-3)+ aH2O(12) * tau0^(-4)/(-4)  + ...
                aH2O(13) * tau0^(-5) /(-5)) + R * aH2O(14);

            H_dep = R * T * (Z - 1) + (T * da_mix_dT - a_mix) / (2*sqrt(2)*b_mix) * ln_term;   % J/mol

            H_mol = y1 * (H0_CO2_mol - H0_CO2_mol_ref) + y2 * (H0_H2O_mol - H0_H2O_mol_ref) + H_dep;   % J/mol
            h = H_mol / Mmix;            % J/kg
        end


        %% Garcia2001 + Spivey2004 LIQUID DENSITY
        function rho_b = garcia2001(obj, P, T, m_NaCl)
            % [x_CO2, x_H2O, x_NaCl]

            aw = [-0.127213, 0.645486, 1.03265,  -0.070291, 0.639589; ...
                4.221,   -3.478,    6.221,     0.5182,  -0.4405; ...
                -11.403,   29.932,   27.952,     0.20684,  0.3768];

            ab = [-7.925e-5, -1.93e-6,  -3.4254e-4,  0.0,      0.0; ...
                1.0998e-3,-2.8755e-3,-3.5819e-3,  -0.72877, 1.92016; ...
                -7.6402e-3, 3.6963e-2, 4.36083e-2, -0.333661,1.185685; ...
                3.746e-4, -3.328e-4, -3.346e-4,   0.0,      0.0; ...
                0.0,       0.0,       0.1353,     0.0,      0.0; ...
                -1.409,    -0.361,    -0.2532,     0.0,      9.216; ...
                0.0,       5.614,     4.6782,    -0.307,    2.6069; ...
                -0.1127,    0.2047,   -0.0452,     0.0,      0.0];

            Pbar = P / 1e5;           % bar
            tc = T - 273.15;          % degC
            tc_100 = tc / 100.0;
            p0 = 700.0;               % bar (reference pressure of 70 MPa)

            % Pure water parameters a_w
            a_w = zeros(3,1);
            for i = 1:3
                num = aw(i,1)*tc_100^2 + aw(i,2)*tc_100 + aw(i,3);
                den = aw(i,4)*tc_100^2 + aw(i,5)*tc_100 + 1.0;
                a_w(i) = num / den;
            end

            rho_w0 = a_w(1);
            Ew     = a_w(2);
            Fw     = a_w(3);

            % Spivey2004 base density rho_b
            if m_NaCl <= 0
                % Pure water
                Iw  = (1.0/Ew) * log(abs(Ew*(Pbar/p0) + Fw));
                Iw0 = (1.0/Ew) * log(abs(Ew*(p0/p0) + Fw));
                rho_b = 1000.0 * rho_w0 * exp(Iw - Iw0);
            else
                % brine
                Cm = m_NaCl;

                a_b = zeros(8,1);
                for i = 1:8
                    num = ab(i,1)*tc_100^2 + ab(i,2)*tc_100 + ab(i,3);
                    den = ab(i,4)*tc_100^2 + ab(i,5)*tc_100 + 1.0;
                    a_b(i) = num / den;
                end

                rho_b0 = rho_w0 + a_b(1) * Cm^2 + a_b(2) * Cm^1.5 + a_b(3) * Cm + a_b(4) * sqrt(Cm);
                Eb = Ew + a_b(5)*Cm;
                Fb = Fw + a_b(6)*Cm^1.5 + a_b(7)*Cm + a_b(8)*sqrt(Cm);

                Ib  = (1.0/Eb) * log(abs(Eb*(Pbar/p0) + Fb));
                Ib0 = (1.0/Eb) * log(abs(Eb*(p0/p0) + Fb));

                rho_b = 1000.0 * rho_b0 * exp(Ib - Ib0);
            end

            % V_app   = (37.51 - 9.585e-2*tc + 8.740e-4*tc^2 - 5.044e-7*tc^3) * 1e-6; % m^3/mol
            % MW_CO2  = 44.01;                                                        % g/mol
            % rho_CO2 = MW_CO2 * 1e-3 / V_app;
            %
            % rho = 1.0 / ( (1.0 - X_CO2) / rho_b + X_CO2 / rho_CO2 );
        end

        %% Fenghour 1998 CO2 VISCOSITY
        function mu = mu_Fenghour1998(obj, T, rho)
            % Fenghour, Wakeham & Vesovic (1998) CO2 viscosity correlation
            a = [0.235156, -0.491266, 5.211155e-2, 5.347906e-2, -1.537102e-2];
            d = [0.4071119e-2, 0.7198037e-4, 0.2411967e-16, 0.2971072e-22, -0.1627888e-22];
            eps_k = 251.196;

            T_ = T / eps_k;
            lnT_ = log(T_);

            lnG = a(1);
            for i = 2:5
                lnG = lnG + a(i) * lnT_^(i-1);
            end
            G = exp(lnG);

            n0 = 1.00697 * sqrt(T) / G; % (microPa.s)

            dn = d(1)*rho + d(2)*rho^2 + d(3)*rho^6/(T_^3) + d(4)*rho^8 + d(5)*rho^8/T_;

            n = (n0 + dn);          % microPa.s
            mu = n * 1e-6;          % microPa.s -> Pa.s
        end

        %% Islam2012 LIQUID VISCOSITY
        function mu_brine = mu_Islam2012(obj, P, T, mNaCl)

            P_bar = P / 1e5; % bar

            % Islam&Carlson (2012): pure water viscosity in microPa.s
            mu_a = 9.03591045e1;
            mu_b = [3.40285740e4, 8.23556123e8, -9.28022905e8];
            mu_c = [1.40090092e-2, 4.86126399e-2, 5.26696663e-2];
            mu_d = [-1.22757462e-1, 2.15995021e-2, -3.65253919e-4, 1.97270835e-6];

            muH2O = mu_a;
            for i = 1:3
                muH2O = muH2O + mu_b(i) * exp(-mu_c(i)*T);
            end
            for i = 0:3
                muH2O = muH2O + (P_bar*0.1) * mu_d(i+1) * (T - 293.15)^i; % 0.1 converts bar->MPa in that fit
            end

            % Mao & Duan (2009): relative viscosity for NaCl solution
            a = [-0.21319213, 0.13651589e-2, -0.12191756e-5];
            b = [0.69161945e-1, -0.27292263e-3, 0.20852448e-6];
            c = [-0.25988855e-2, 0.77989227e-5];

            A = a(1) + a(2)*T + a(3)*T^2;
            B = b(1) + b(2)*T + b(3)*T^2;
            C = c(1) + c(2)*T;

            m = mNaCl;
            mu_r = exp(A*m + B*m^2 + C*m^3);

            mu_brine = (mu_r * muH2O) * 1e-6; % microPa.s -> Pa.s

            % Islam & Carlson CO2 correction
            % mu = mu_brine * (1.0 + 4.65 * x1^1.0134);
        end

        %% Ziabakhsh2012 fugacity coefficient (aq) + Guo2019 ENTHALPY
        function [phi1, phi2, hb, hCO2] = zia_phi_aq_CO2_H2O(obj, P, T, mNa, mCl, varargin)
            have_x = (numel(varargin) >= 2);

            P_bar = P / 1e5;  % bar
            tc = T - 273.15;  % degC

            % constants
            R = 83.145;       % bar*cm3/mol/K
            Vw = 18.1;        % cm3/mol
            Mw = 18.0152;     % g/mol H2O
            Tc = 647.14;      % K
            Pc = 220.64;      % bar

            PswC = [-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502];

            % lambda/kappa coefficients for CO2
            lab = [-0.0652869, 1.6790636e-4, 40.838951, 0, 0, -3.9266518e-2, 0, 2.1157167e-2, 6.5486487e-6, 0, 0, 0];
            ksi = [-1.144624e-2, 2.8274958e-5, 1.3980876e-2, -1.4349005e-2];

            % CO2 Henry term lnk_H
            eta_CO2   = -0.114535;
            tau_CO2   = -5.279063;
            beta_CO2  = 6.187967;
            Gamma_CO2 = 0.0;

            % Eq.A2
            denom = 1./(0.9998396 + 18.224944e-3 * tc - 7.922210e-6  * tc.^2 - 55.44846e-9  * tc.^3 + 149.7562e-12 * tc.^4 - 393.2952e-15 * tc.^5);
            V0 = (1.0 + 18.159725e-3 * tc) .* denom;   % cm3/g

            ddenom_dtc = - (18.224944e-3 - 2*7.922210e-6*tc - 3*55.44846e-9*tc^2 + 4*149.7562e-12*tc^3 - 5*393.2952e-15*tc^4) / (0.9998396 + 18.224944e-3*tc - 7.922210e-6*tc^2 - 55.44846e-9*tc^3 ...
                + 149.7562e-12*tc^4 - 393.2952e-15*tc^5 )^2;
            dV0_dT = 18.159725e-3 * denom + ( 1 + 18.159725e-3*tc) * ddenom_dtc;

            B  = 19654.320 + 147.037 * tc - 2.21554 * tc.^2 + 1.0478e-2 * tc.^3 - 2.2789e-5 * tc.^4;
            A1 = 3.2891   - 2.3910e-3 * tc + 2.8446e-4 * tc.^2 - 2.8200e-6 * tc.^3 + 8.477e-9 * tc.^4;
            A2 = 6.245e-5 - 3.913e-6 * tc - 3.499e-8 * tc.^2 + 7.942e-10 * tc.^3 - 3.299e-12 * tc.^4;

            dB_dT  = 147.037 - 2*2.21554*tc + 3*1.0478e-2*tc^2 - 4*2.2789e-5*tc^3;
            dA1_dT = -2.3910e-3 + 2*2.8446e-4*tc - 3*2.8200e-6*tc^2 + 4*8.477e-9*tc^3;
            dA2_dT = -3.913e-6 - 2*3.499e-8*tc + 3*7.942e-10*tc^2 - 4*3.299e-12*tc^3;

            denomP = 1 ./ (B + A1 * (P_bar - 1.0) + A2 * (P_bar-1.0) .^ 2);
            V = V0 - V0 .* (P_bar - 1.0)  .* denomP;     % cm3/g
            rho0_H2O = 1 ./ V;                           % g/cm3

            ddenomP_dT = - (dB_dT + dA1_dT * (P_bar - 1.0) + dA2_dT * (P_bar - 1.0) ^ 2) / (B + A1 * (P_bar - 1.0) + A2 * (P_bar - 1.0) ^ 2) ^ 2;
            dV_dT = dV0_dT - (dV0_dT * (P_bar - 1.0) * denomP + V0 * (P_bar - 1.0) * ddenomP_dT);
            drho0_dT = -(1/V^2)*dV_dT;

            tau = 1 - T / Tc;
            Psw = Pc * exp( Tc / T * (PswC(1)*tau + PswC(2)*tau^(1.5) + PswC(3)*tau^3 + PswC(4)*tau^(3.5) + PswC(5)*tau^4 + PswC(6)*tau^(7.5)));
            dPsw_dT = Psw * (-(Tc/T^2) * (PswC(1)*tau + PswC(2)*tau^(1.5) + PswC(3)*tau^3 + PswC(4)*tau^(3.5) + PswC(5)*tau^4 + PswC(6)*tau^(7.5)) - (1/T)*(PswC(1) + 1.5*PswC(2)*tau^(0.5) + 3*PswC(3)*tau^2 + 3.5*PswC(4)*tau^(2.5) + 4*PswC(5)*tau^3 + 7.5*PswC(6)*tau^(6.5)));

            % f0_H2O (bar)
            f0_H2O = Psw * exp( (P_bar - Psw) * Mw * V / (R * T) );  % bar
            df0_dT = f0_H2O * (dPsw_dT / Psw + (-(dPsw_dT) * Mw * V/(R*T) + (P_bar - Psw)*Mw*dV_dT/(R*T) - (P_bar - Psw)*Mw*V/(R*T^2)));  

            dB  = tau_CO2 + Gamma_CO2 * P_bar + beta_CO2 * sqrt(1000.0 / T);
            dBdT = beta_CO2 * (-0.5) * sqrt(1000.0) * T^(-3/2);
            lnkH = (1 - eta_CO2) * log(f0_H2O) + eta_CO2 * log( (R*T/Mw) * rho0_H2O ) + 2.0 * rho0_H2O * dB;
            dlnkH_dT = (1 - eta_CO2) * (df0_dT / f0_H2O) ...
                + eta_CO2 * ( ((R/Mw)*rho0_H2O + (R*T/Mw)*drho0_dT) / ((R*T/Mw)*rho0_H2O) ) ...
                + 2 * (drho0_dT * dB + rho0_H2O * dBdT);
            mc  = mNa;
            mac = mNa * mCl;

            Tinv = 1/T; pinv = 1/P_bar;

            lambda = lab(1) + lab(2)*T + lab(3)*Tinv + lab(4)*P_bar + lab(5)*pinv + lab(6)*P_bar*Tinv + ...
                lab(7)*T*(P_bar^(-2)) + lab(8)*P_bar/(630-T) + lab(9)*T*log(P_bar) + ...
                lab(10)*P_bar*(T^(-2)) + lab(11)*(P_bar^2)*T + lab(12)*P_bar*T;

            kappa = ksi(1) + ksi(2)*T + ksi(3)*P_bar*Tinv + ksi(4)*P_bar/(630-T);

            lnji = 2.0 * mc * lambda + mac * kappa;

            dlambda_dT = lab(2) - lab(3)/T^2 - lab(6)*P_bar/T^2 + lab(7)/P_bar^2 + lab(8)*P_bar/(630 - T)^2 ...
                + lab(9)*log(P_bar) - 2*lab(10)*P_bar/T^3 + lab(11)*P_bar^2 + lab(12)*P_bar;
            dksi_dT = ksi(2) - ksi(3)*P_bar/T^2 + ksi(4)*P_bar/(630 - T)^2;
            dlnji_dT = 2*mc*dlambda_dT + mac*dksi_dT;

            % lnphi_aq_CO2 = lnji + lnkH - ln(p)
            lnphiCO2 = lnji + lnkH - log(P_bar);
            

            % --- water term: lnphi_aq_H2O = lnKw - ln(p) ---
            log10K0w = -2.209 + 3.097e-2*tc - 1.098e-4*tc^2 + 2.048e-7*tc^3;
            K0w = 10^(log10K0w);
            Kw = K0w * exp((P_bar-1.0)*Vw/(R*T));

            dlog10K0w_dT = 3.097e-2 - 2*1.098e-4*tc + 3*2.048e-7*tc^2; % d/dT since dtc/dT=1
            dlnphiH2O_dT = log(10)*dlog10K0w_dT - (P_bar-1.0)*Vw/(R*T^2);

            lnphiH2O = log(Kw) - log(P_bar);
            phi1 = exp(lnphiCO2);
            phi2 = exp(lnphiH2O);
            hb = NaN;
            hCO2 = NaN;

            if have_x
                T_star = 1000; % K
                T_refe = 298.15;
                tau = T / T_star;
                tau0 = T_refe / T_star;

                Ru = 8.314462618;   % J/mol/K
                Mw_CO2 = 44.0e-3;   % kg/mol
                Mw_H2O = 18.0152e-3;   % kg/mol

                aCO2 = [ -1.8188731, 12.903022, -9.6634864, 4.2251879, -1.042164, 0.12683515, -0.49939675E-02, 2.4950242, -0.8272375, 0.15372481, ...
                    -0.015861243, 8.601715E-04, -0.19222165E-04, 2.10820700252E+03];

                % H2O coefficients from Aleksandrov2005
                aH2O = [ 3.10409601236035E+01, -3.91422080460869E+01, 3.79695277233575E+01, -2.18374910952284E+01, 7.42251494566339E+00, -1.38178929609470E+00, ...
                    1.08807067571454E-01, -1.20771176848589E+01, 3.39105078851732E+00, -5.84520979955060E-01, 5.89930846488082E-02, -3.12970001415882E-03, 6.57460740981757E-05, -1.249994285621061E+04];

                H0_CO2_mol = Ru * T_star * ( aCO2(1)  * tau.^1   / 1  + aCO2(2)  * tau.^2   / 2  + aCO2(3)  * tau.^3   / 3  + ...
                    aCO2(4)  * tau.^4   / 4  + aCO2(5)  * tau.^5   / 5  + aCO2(6)  * tau.^6   / 6  + aCO2(7)  * tau.^7   / 7  + ...
                    aCO2(8)  * log(tau)      + aCO2(9)  * tau.^(-1)/(-1)+ aCO2(10) * tau.^(-2)/(-2)+ aCO2(11) * tau.^(-3)/(-3)+ ...
                    aCO2(12) * tau.^(-4)/(-4)+ aCO2(13) * tau.^(-5)/(-5) ) + Ru * aCO2(14) ;%+ 1.2897e+04;   % J/mol (NIST shift)

                H0_CO2_mol_ref = Ru * T_star * ( aCO2(1)  * tau0^1    / 1  + aCO2(2)  * tau0^2    / 2  + aCO2(3)  * tau0^3    / 3  + aCO2(4)  * tau0^4    / 4  + ...
                    aCO2(5)  * tau0^5    / 5  + aCO2(6)  * tau0^6    / 6  + aCO2(7)  * tau0^7    / 7  + aCO2(8)  * log(tau0)      + ...
                    aCO2(9)  * tau0^(-1) /(-1)+ aCO2(10) * tau0^(-2)/(-2)+ aCO2(11) * tau0^(-3) /(-3)+ aCO2(12) * tau0^(-4)/(-4)+ ...
                    aCO2(13) * tau0^(-5) /(-5)) + Ru * aCO2(14);

                H0_H2O_mol = Ru * T_star * ( aH2O(1)  * tau.^1   / 1  + aH2O(2)  * tau.^2   / 2  + aH2O(3)  * tau.^3   / 3  + ...
                    aH2O(4)  * tau.^4   / 4  + aH2O(5)  * tau.^5   / 5  + aH2O(6)  * tau.^6   / 6  + aH2O(7)  * tau.^7   / 7  + ...
                    aH2O(8)  * log(tau)      + aH2O(9)  * tau.^(-1)/(-1)+ aH2O(10) * tau.^(-2)/(-2)+ aH2O(11) * tau.^(-3)/(-3)+ ...
                    aH2O(12) * tau.^(-4)/(-4)+ aH2O(13) * tau.^(-5)/(-5) ) + Ru * aH2O(14) ;%- 8.0178e+03;     % J/mol (NIST shift)

                H0_H2O_mol_ref = Ru * T_star * ( aH2O(1)  * tau0^1    / 1  + aH2O(2)  * tau0^2    / 2  + aH2O(3)  * tau0^3    / 3  + aH2O(4)  * tau0^4    / 4  + ...
                    aH2O(5)  * tau0^5    / 5  + aH2O(6)  * tau0^6    / 6  + aH2O(7)  * tau0^7    / 7  + aH2O(8)  * log(tau0)      + ...
                    aH2O(9)  * tau0^(-1) /(-1)+ aH2O(10) * tau0^(-2)/(-2)+ aH2O(11) * tau0^(-3) /(-3)+ aH2O(12) * tau0^(-4)/(-4)  + ...
                    aH2O(13) * tau0^(-5) /(-5)) + Ru * aH2O(14);

                HCO2 = (H0_CO2_mol - H0_CO2_mol_ref) ./Mw_CO2;
                Hw   = (H0_H2O_mol - H0_H2O_mol_ref) ./Mw_H2O;

                hdep_CO2 = -Ru * T^2 * (dlnkH_dT + dlnji_dT) ./ Mw_CO2;
                hdep_H2O = -Ru * T^2 * dlnphiH2O_dT / Mw_H2O;
                
                hb = Hw + hdep_H2O;          % J/kg brine
                hCO2 = HCO2 + hdep_CO2;

            end

        end
        

        function dA = fd3D_axis(obj, A, grid, axis_id)
            dA = nan(size(A));

            if axis_id == 1
                dA(1,:,:) = (A(2,:,:) - A(1,:,:)) ./ (grid(2) - grid(1));
                dA(end,:,:) = (A(end,:,:) - A(end-1,:,:)) ./ (grid(end) - grid(end-1));
                dg = (grid(3:end) - grid(1:end-2)); % (n-2)x1
                dA(2:end-1,:,:) = (A(3:end,:,:) - A(1:end-2,:,:)) ./ reshape(dg, [], 1, 1);

            elseif axis_id == 2
                dA(:,1,:) = (A(:,2,:) - A(:,1,:)) ./ (grid(2) - grid(1));
                dA(:,end,:) = (A(:,end,:) - A(:,end-1,:)) ./ (grid(end) - grid(end-1));
                dg = (grid(3:end) - grid(1:end-2)).'; % 1x(n-2)
                dA(:,2:end-1,:) = (A(:,3:end,:) - A(:,1:end-2,:)) ./ reshape(dg, 1, [], 1);

            elseif axis_id == 3
                dA(:,:,1) = (A(:,:,2) - A(:,:,1)) ./ (grid(2) - grid(1));
                dA(:,:,end) = (A(:,:,end) - A(:,:,end-1)) ./ (grid(end) - grid(end-1));
                dg = (grid(3:end) - grid(1:end-2)).'; % 1x(n-2)
                dA(:,:,2:end-1) = (A(:,:,3:end) - A(:,:,1:end-2)) ./ reshape(dg, 1, 1, []);
            else
                error('axis_id must be 1,2,3');
            end
        end

    end
end
