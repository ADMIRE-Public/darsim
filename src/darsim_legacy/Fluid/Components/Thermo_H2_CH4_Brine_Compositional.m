% Thermodynamics of H2-CH4-Brine for compositional Fluid Model class
%--------------------------------------------------------------------------
% DARSim 2 Reservoir Simulator
% Author: Amin MB
% TU Delft
% Created: 22 July 2025
%--------------------------------------------------------------------------
classdef Thermo_H2_CH4_Brine_Compositional < matlab.mixin.Heterogeneous & handle
    %{
This MATLAB class calculates thermodynamic properties and phase 
equilibria for a system containing hydrogen (H2), methane (CH4), and brine
(H2O with salinity). The code is based on the Peng-Robinson Equation of 
State (PR EOS) and follows the methodology described by 
Ziabakhsh-Ganji and Kooi (2012).

References:

1. Ziabakhsh-Ganji, Z., & Kooi, H. (2012). An Equation of State for 
thermodynamic equilibrium of gas mixtures and brines to allow simulation 
of the effects of impurities in subsurface CO2 storage. International 
Journal of Greenhouse Gas Control, 11, S21-S34.

2. Shabani, B., & Vilcáez, J. (2019). TOUGHREACT-CO2Bio–A new module to 
simulate geological carbon storage under biotic conditions (Part 1): 
The multiphase flow of CO2-CH4-H2-H2S gas mixtures. Journal of Natural 
Gas Science and Engineering, 63, 85-94.

3. Nielsen, M. H., Arekhov, V., Whitson, C. H., Clemens, T., Zhainakov, 
T., & Wegner, J. (2023, June). Fluid Modeling of Underground Hydrogen 
Storage in a Depleted Natural Gas Field. In SPE Europec featured at EAGE 
Conference and Exhibition (p. D021S001R002). SPE.

INPUTS:
    T: Temperature
    Salinity: Salinity of brine

OUTPUTS:
    gridded interpolant object for K-values and compressibility factor


    %}
    properties
        comp % a cell array containing components' data

        zcf
        dzcfdp
        dzcfdy1

        mu % Viscosity
        dmudp
        dmudy1
    end
    methods
        function obj = Thermo_H2_CH4_Brine_Compositional(T, Salinity)
            obj.ComputeLookupTable(T, Salinity);
        end
        function ComputeLookupTable(obj, T, Salinity)
            %% Thermodynamics of Hydrogen-Methan-brine
            fprintf('##### Calculating the thermodynamics of H2-CH4-brine...\n');

            %% Initial inputs and functions

            R = 83.1447; % Gas constant [bar cm^3 mol^(-1) K^(-1)]

            % Defining
            P_vec = 10:10:500; % Pressure range [bar]
            y1_vec = [1e-8, 1e-4:1e-4:1-1e-4, 1-1e-8]; % Mole fraction of H2 in vapor phase
            y2_vec = 1 - y1_vec;   % Mole fraction of CH4 in vapor phase
            fprintf('##### Applicable pressure range: %d - %d [bar]\n', min(P_vec), max(P_vec));
            fprintf('##### You can change pressure range in Fluid/Components. \n');

            P = zeros(length(P_vec), length(y1_vec));
            y1 = zeros(length(P_vec), length(y1_vec));

            for i = 1:length(y1_vec)
                y1(:, i) = y1_vec(i);
            end
            for i = 1:length(P_vec)
                P(i, :) = P_vec(i);
            end
            y2 = 1 - y1;
            T = ones(length(P_vec), length(y1_vec)) * T;

            % Temperature
            T_C = T - 273.15; % Temperature in [C]
            P_app = P - 1.01325 ; % Applied pressure at P_app = 0 bars, the absolute pressure is 1 atm or 1.01325 bars

            % brine
            % Salinity is in [ppm] (1000 [ppm] means 1 [g] solute in 1 [kg] water)
            m_Na = Salinity / 1e3 / 58.44; % [mol/kg]
            m_Cl = m_Na; % [mol/kg]

            % H2
            H2.Pc = 12.964; % Critical pressure [bar]
            H2.Tc = 33.145; % Critical temperature [K]
            H2.w = -0.219;  % Acentric factor
            H2.Mw = 2.0157; % Molecular weight [g/mol]
            H2.c = -2.959;  % Volume shift factor (from Shabani (2019))
            H2.X = 4.8070; % Constants fitting parameter for correct density
            H2.Vc = 64.4828; % the molar volume in the critical point. cm3/mol

            % CH4
            CH4.Pc = 45.992;  % Critical pressure [bar]
            CH4.Tc = 190.564; % Critical temperature [K]
            CH4.w = 0.01142;  % Acentric factor
            CH4.Mw = 16.0425; % Molecular weight [g/mol]
            CH4.c = -0.1557;  % Volume shift factor (from Shabani (2019))
            CH4.X = 6.404378; % Constants fitting parameter for correct density calculation for PR EOS.
            CH4.Vc = 99.3; % the molar volume in the critical point. cm3/mol

            % H2O
            H2O.Pc = 220.640; % Critical pressure [bar]
            H2O.Tc = 647.096; % Critical temperature [K]
            H2O.w = 0.3443;   % Acentric factor
            H2O.Mw = 18.0152; % Molecular weight [g/mol]
            H2O.Vc = 55.9478; % the molar volume in the critical point. cm3/mol

            %% Thermodynamics
            % Peng-Robinson EOS
            alpha = @(T, Tc, w) (1 + (0.37464 + 1.54226*w - 0.26992*w^2) * (1 - sqrt(T / Tc))).^2;
            a = @(Tc, Pc, alpha) (0.45724 * R^2 * Tc^2 * alpha) / Pc;
            b = @(Tc, Pc) 0.07780 * R * Tc / Pc;

            H2.alpha = alpha(T, H2.Tc, H2.w);
            H2.a = a(H2.Tc, H2.Pc, H2.alpha);
            H2.b = b(H2.Tc, H2.Pc);
            H2.B = H2.b .* P ./ (R * T);

            CH4.alpha = alpha(T, CH4.Tc, CH4.w);
            CH4.a = a(CH4.Tc, CH4.Pc, CH4.alpha);
            CH4.b = b(CH4.Tc, CH4.Pc);
            CH4.B = CH4.b .* P ./ (R * T);


            H2O.alpha = alpha(T, H2O.Tc, H2O.w);
            H2O.a = a(H2O.Tc, H2O.Pc, H2O.alpha);
            H2O.b = b(H2O.Tc, H2O.Pc);
            H2O.B = H2O.b .* P ./ (R * T);

            % Mixture rule
            mixture.a11 = H2.a;
            mixture.a22 = CH4.a;
            mixture.a12 = sqrt(H2.a .* CH4.a) * (1 - 0.362);
            mixture.b1 = H2.b;
            mixture.b2 = CH4.b;

            mixture.a = y1.^2 .* mixture.a11 + 2 * y1 .* y2 .* mixture.a12 + y2.^2 .* mixture.a22;
            mixture.b = y1 * mixture.b1 + y2 * mixture.b2;

            mixture.A = mixture.a .* P ./ (R*T).^2;
            mixture.B = mixture.b .* P ./ (R*T);

            % Compressibility factor
            Z_guess = 1; % Initial guess
            Z = arrayfun(@(A, B)  fzero(@(Z) Z^3 - (1 - B) * Z^2 + (A - 2*B - 3*B^2) * Z - (A*B - B^2 - B^3), Z_guess),  mixture.A, mixture.B);
            if any(Z <= 0, 'all') || any(Z <= mixture.B, 'all')
                error('There is a problem in calculation of compressibility factor (Z).')
            end

            % Fugacity coefficient
            H2.phi = exp( ...
                H2.B ./ mixture.B .* (Z - 1) ...
                - log(Z - mixture.B) ...
                + (mixture.A ./ (2.828 * mixture.B)) .* ...
                (H2.B ./ mixture.B - (2 * (y1 .* mixture.a11 + y2 .* mixture.a12)) ./ mixture.a) .* ...
                (log((Z + 2.414 * mixture.B) ./ (Z - 0.414 * mixture.B))) ...
                );


            CH4.phi = exp( ...
                CH4.B ./ mixture.B .* (Z - 1) ...
                - log(Z - mixture.B) ...
                + (mixture.A ./ (2.828 * mixture.B)) .* ...
                (CH4.B ./ mixture.B - (2 * (y1 .* mixture.a12 + y2 .* mixture.a22)) ./ mixture.a) .* ...
                (log((Z + 2.414 * mixture.B) ./ (Z - 0.414 * mixture.B))) ...
                );

            % Properties of pure water
            V0 = (1 + 18.159725E-3 * T_C) ./ (0.9998396 + 18.224944E-3 * T_C - 7.922210E-6 * T_C.^2 - 55.44846E-9 * T_C.^3 + 149.7562E-12 * T_C.^4 - 393.2952E-15 * T_C.^5);
            B  = 19654.320 + 147.037   * T_C - 2.21554   * T_C.^2 + 1.0478E-2  * T_C.^3 - 2.2789E-5 * T_C.^4;
            A1 = 3.2891    - 2.3910E-3 * T_C + 2.8446E-4 * T_C.^2 - 2.8200E-6  * T_C.^3 + 8.477E-9  * T_C.^4;
            A2 = 6.2450E-5 - 3.9130E-6 * T_C + 3.4990E-8 * T_C.^2 - 7.9420E-10 * T_C.^3 + 3.299E-12 * T_C.^4;

            V = V0 - (V0 .* P_app) ./ (B + A1 .* P_app + A2 .* P_app.^2); % [cm3/g], typo in Ganji2012
            H2O.rho = 1 ./ V; % [g/cm3]

            % Fugucity of pure water (Appendix A)
            V_molar = 18.0152 * V;
            tau = 1 - T/H2O.Tc;

            % vapor pressure of pure water, only applicable up to critical point
            P_s = H2O.Pc .* exp(H2O.Tc ./ T .*(-7.85951783 * tau + 1.84408259 * tau.^1.5 - 11.7866497 * tau.^3 + 22.6807411 * tau.^3.5 - 15.9618719 * tau.^4 + 1.80122502 * tau.^7.5));
            H2O.f0 = P_s .* exp((P - P_s) .* V_molar ./ (R*T)); % Poynting equation
            % H2O.phi = H2O.f ./ P;
            H2O.phi = exp( ...
                H2O.B ./ mixture.B .* (Z - 1) ...
                - log(Z - mixture.B) ...
                + (mixture.A ./ (2.828 * mixture.B)) .* ...
                (H2O.B ./ mixture.B - (2 * (y1 .* sqrt(H2.a .* H2O.a) .* (1 - 0.1) + y2 .* sqrt(CH4.a .* H2O.a) .* (1 - 0.47893))) ./ mixture.a) .* ...
                (log((Z + 2.414 * mixture.B) ./ (Z - 0.414 * mixture.B))) ...
                );


            % Henry's constant
            eta = 0.294994;
            tau = -4.23407;
            beta = 5.41551;
            H2.deltaB = tau + beta * sqrt(1E3 ./ T);
            H2.KH = exp((1 - eta) * log(H2O.f0) + eta * log(R * T / H2O.Mw .* H2O.rho) + 2 * H2O.rho .* H2.deltaB);

            eta = -0.092248;
            tau = -5.779280;
            beta = 7.262730;
            CH4.deltaB = tau + beta * sqrt(1E3 ./ T);
            CH4.KH = exp((1 - eta) * log(H2O.f0) + eta * log(R * T / H2O.Mw .* H2O.rho) + 2 * H2O.rho .* CH4.deltaB);

            % Activity coefficient
            Par = @(c, T, P) c(1) + c(2) * T + c(3) ./ T + c(4) * P + c(5) ./ P + c(6) * P ./ T + c(7) * T ./ P.^2 + c(8) * P ./ (630 - T) + c(9) * T .* log(P) + c(10) * P ./ T.^2;

            H2.c_Na =     [-2.1432831E+00 +3.1411257E-03 +3.9220546E+02 -2.8601200E-05 +0.0000000E+00 +2.3527160E-03 +0.0000000E+00 -2.4422000E-03 +2.9806000E-06 -3.8900000E-01];
            H2.c_Na_Cl =  [-4.0631000E-03 +0.0000000E+00 +0.0000000E+00 -3.6650000E-06 +1.7004000E-01 -4.1800000E-04 +5.8800000E-04 -3.9100000E-04 +0.0000000E+00 +1.9913000E-01];
            H2.gamma = exp(2 * m_Na * Par(H2.c_Na, T, P) + m_Na * m_Cl * Par(H2.c_Na_Cl, T, P));

            CH4.c_Na =    [-5.7066455E-01 +7.2997588E-04 +1.5176903E+02 +3.1927112E-05 +0.0000000E+00 -1.6426510E-05 +0.0000000E+00 +0.0000000E+00 +0.0000000E+00 +0.0000000E+00];
            CH4.c_Na_Cl = [-2.9990084E-03 +0.0000000E+00 +0.0000000E+00 +0.0000000E+00 +0.0000000E+00 +0.0000000E+00 +0.0000000E+00 +0.0000000E+00 +0.0000000E+00 +0.0000000E+00];
            CH4.gamma = exp(2 * m_Na * Par(CH4.c_Na, T, P) + m_Na * m_Cl * Par(CH4.c_Na_Cl, T, P));

            % K-value
            H2.K_mole = H2.KH .* H2.gamma ./ (H2.phi .* P);
            CH4.K_mole = CH4.KH .* CH4.gamma ./ (CH4.phi.* P);

            H2O.K0 = 10 .^ (-2.209 + 3.097E-2 * T_C - 1.098E-4 * T_C.^2 + 2.048E-7 * T_C.^3);
            p_ref = 1; % [bar];
            H2O.Vm_avg = 18.1; % [cm^3 mol^(-1)];
            H2O.K_mole = H2O.K0 .* exp((P - p_ref) * H2O.Vm_avg ./ (R * T))./ (H2O.phi .* P); %(eq 11)


            %% Viscosity
            
            % ==== H2
            V_H2 = Z(:, 1) .* R .* T ./ P;
            d1 = -1.55199;
            d2 =  2.92788;
            d3 =  0.645731;
            eta0_H2 = d1 * (T .^ 0.5) + d2 * (T .^ d3);
            Pa_H2 = -H2.a ./( V_H2.^2 + 2*H2.b*V_H2 - H2.b.^2 );
            Pr_H2 = R * T ./(V_H2 - H2.b);
            kr_H2 = -0.00185308;
            ka_H2 = -0.332575;
            krr_H2 = (1.35146e-4)./((T./H2.Tc).^2);        
            % etaf_H2 = ka_H2 .* Pa_H2 + kr_H2 .* Pr_H2 + krr_H2 .* Pr_H2 .^ 2;
            % mu_H2 = (eta0_H2 + etaf_H2); % uP

            % ==== CH4
            V_CH4 = Z(:, end) .* R .* T ./ P;
            d1 = 13.3919;
            d2 = -47.9429 ;
            d3 = 0.160913;
            eta0_CH4 = d1 * (T .^ 0.5) + d2 * (T .^ d3);
            Pa_CH4 = -CH4.a ./( V_CH4.^2 + 2*CH4.b*V_CH4 - CH4.b.^2 );
            Pr_CH4 = R * T ./(V_CH4 - CH4.b);
            kr_CH4 = 0.0731796; 
            ka_CH4 = -0.382909;
            krr_CH4 = (6.63615e-5)./((T./CH4.Tc).^2);
            % etaf_CH4 = ka_CH4 .* Pa_CH4 + kr_CH4 .* Pr_CH4 + krr_CH4 .* Pr_CH4 .^ 2;
            % mu_CH4 = (eta0_CH4 + etaf_CH4);

            % eta0 Zeberg Mikkelsen 2002 Wilke
            a1 = (1+((eta0_H2./eta0_CH4).^0.5).*(CH4.Mw./H2.Mw).^0.25).^2;
            b1 = (4/(2.^0.5)).*(1+(H2.Mw/CH4.Mw)).^0.5; 
            phi_H2_CH4 = a1./b1;
            
            a2 = (1+((eta0_CH4./eta0_H2).^0.5).*(H2.Mw./CH4.Mw).^0.25).^2;
            b2 = (4/(2.^0.5)).*(1+(CH4.Mw/H2.Mw)).^0.5; 
            phi_CH4_H2 = a2 ./ b2;
            phi_H2_H2 = 1;
            phi_CH4_CH4 = 1;
            
            % Mixture
            V_mixture = Z .* R .* T ./ P;
            eta0_mix= (y1.*eta0_H2) ./ (y1.*phi_H2_H2 + y2.*phi_H2_CH4) + (y2.*eta0_CH4) ./ (y1.*phi_CH4_H2 + y2.*phi_CH4_CH4);
            
            eps=0.3;
            MM = y1./(H2.Mw.^eps) + y2./(CH4.Mw.^eps);
            Z2 = y1 ./((H2.Mw.^eps).*MM);
            
            krI = Z2.*kr_H2 + (1-Z2).*kr_CH4 ;
            kaI = Z2.*ka_H2 + (1-Z2).*ka_CH4;
            krrI = Z2.*krr_H2 + (1-Z2).*krr_CH4;
            
            Pa_mix = -mixture.a ./ (V_mixture .^2 + 2*mixture.b .* V_mixture - mixture.b.^2);
            Pr_mix = R * T ./ (V_mixture - mixture.b);
            
            etaf_mix = krI.*Pr_mix + kaI.*Pa_mix + krrI.*Pr_mix.^2;
            
            mu_mix = eta0_mix + etaf_mix;
            mu_mix = mu_mix * 1e-7; % Converting unti to [Pa.s]
            clear d1 d2 d3 eta0_H2 Pa_H2 Pr_H2 kr_H2 ka_H2 krr_H2 etaf_H2 d1 d2 d3 eta0_CH4 Pa_CH4 Pr_CH4 kr_CH4 ka_CH4 krr_CH4 etaf_CH4 V_mixture eta0_mix eps MM Z2 krI kaI krrI Pa_mix Pr_mix etaf_mix

            %% creating interpolations
            fprintf('##### Final step: creating interpolations!\n');

            % H2.dKdp  = ( H2.K(2:end, :) -  H2.K(1:end-1, :)) ./( P(2:end, :) * 1E5 - P(1:end-1, :) * 1E5);
            % CH4.dKdp = (CH4.K(2:end, :) - CH4.K(1:end-1, :)) ./( P(2:end, :) * 1E5 - P(1:end-1, :) * 1E5);
            % H2O.dKdp = (H2O.K(2:end, :) - H2O.K(1:end-1, :)) ./( P(2:end, :) * 1E5 - P(1:end-1, :) * 1E5);
            % dZdp = (Z(2:end, :) - Z(1:end-1, :)) ./(P(2:end, :) * 1E5 - P(1:end-1, :) * 1E5);
            %
            % H2.dKdy1  = ( H2.K(:, 2:end) -  H2.K(:, 1:end-1)) ./( y1(:, 2:end) - y1(:, 1:end-1));
            % CH4.dKdy1 = (CH4.K(:, 2:end) - CH4.K(:, 1:end-1)) ./( y1(:, 2:end) - y1(:, 1:end-1));
            % H2O.dKdy1 = (H2O.K(:, 2:end) - H2O.K(:, 1:end-1)) ./( y1(:, 2:end) - y1(:, 1:end-1));
            % dZdy1 = (Z(:, 2:end) - Z(:, 1:end-1)) ./(y1(:, 2:end) - y1(:, 1:end-1));
            %
            % H2.K  = griddedInterpolant(P * 1E5, y1, H2.K);
            % CH4.K = griddedInterpolant(P * 1E5, y1, CH4.K);
            % H2O.K = griddedInterpolant(P * 1E5, y1, H2O.K);
            % obj.zcf = griddedInterpolant(P * 1E5, y1, Z);
            %
            % P_mid = (P(1:end-1, :) + (P(2:end, :))) / 2;
            % H2.dKdp  = griddedInterpolant(P_mid * 1E5, y1(2:end, :), H2.dKdp);
            % CH4.dKdp = griddedInterpolant(P_mid * 1E5, y1(2:end, :), CH4.dKdp);
            % H2O.dKdp = griddedInterpolant(P_mid * 1E5, y1(2:end, :), H2O.dKdp);
            % obj.dzcfdp = griddedInterpolant(P_mid, y1(2:end, :), dZdp);
            %
            % y1_mid = (y1(:, 1:end-1) + (y1(:, 2:end))) / 2;
            % H2.dKdy1  = griddedInterpolant(P(:, 2:end) * 1E5, y1_mid, H2.dKdy1);
            % CH4.dKdy1 = griddedInterpolant(P(:, 2:end) * 1E5, y1_mid, CH4.dKdy1);
            % H2O.dKdy1 = griddedInterpolant(P(:, 2:end) * 1E5, y1_mid, H2O.dKdy1);
            % obj.dzcfdy1 = griddedInterpolant(P(:, 2:end) * 1E5, y1_mid, dZdy1);

            % Save y1 as an input for interpolations
            y1_ = y1;

            % Update y_i, x_i based on mole fraction of brine in vapor phase
            % (doing an averaging for accuracy!)
            x1 = y1 ./ H2.K_mole; x2 = y2 ./ CH4.K_mole; x3 = 1 - x1 - x2;
            sum_x = x1 + x2 + x3;
            x1 = x1 ./ sum_x; x2 = x2 ./ sum_x; x3 = x3 ./ sum_x;

            y3 = x3 .* H2O.K_mole;
            sum_y = y1 + y2 + y3;
            y1 = y1 ./ sum_y; y2 = y2 ./ sum_y; y3 = y3 ./ sum_y;


            % Creating interpolations of K_mole-values (these will be used
            % in flash calculations)
            input_1 = P * 1E5;
            input_2 = y1_;
             H2.K_mole = griddedInterpolant(input_1, input_2,  H2.K_mole);
            CH4.K_mole = griddedInterpolant(input_1, input_2, CH4.K_mole);
            H2O.K_mole = griddedInterpolant(input_1, input_2, H2O.K_mole);

            % Calculating mass fractions
            Y1 = y1 .*  H2.Mw ./ (y1 .* H2.Mw + y2 .* CH4.Mw + y3 .* H2O.Mw);
            Y2 = y2 .* CH4.Mw ./ (y1 .* H2.Mw + y2 .* CH4.Mw + y3 .* H2O.Mw);
            Y3 = y3 .* H2O.Mw ./ (y1 .* H2.Mw + y2 .* CH4.Mw + y3 .* H2O.Mw);
            X1 = x1 .*  H2.Mw ./ (x1 .* H2.Mw + x2 .* CH4.Mw + x3 .* H2O.Mw);
            X2 = x2 .* CH4.Mw ./ (x1 .* H2.Mw + x2 .* CH4.Mw + x3 .* H2O.Mw);
            X3 = x3 .* H2O.Mw ./ (x1 .* H2.Mw + x2 .* CH4.Mw + x3 .* H2O.Mw); 

            % calculating K_mass-values and its derivatives(these will be used
            % in calculation of mass fraction derivatives and compressibility) 
             H2.K_mass = Y1 ./ X1;
            CH4.K_mass = Y2 ./ X2;
            H2O.K_mass = Y3 ./ X3;

            % Derivaties of K_mass-values
             H2.dKdp = ( H2.K_mass(2:end, :) -  H2.K_mass(1:end-1, :)) ./( P(2:end, :) * 1E5 - P(1:end-1, :) * 1E5);
            CH4.dKdp = (CH4.K_mass(2:end, :) - CH4.K_mass(1:end-1, :)) ./( P(2:end, :) * 1E5 - P(1:end-1, :) * 1E5);
            H2O.dKdp = (H2O.K_mass(2:end, :) - H2O.K_mass(1:end-1, :)) ./( P(2:end, :) * 1E5 - P(1:end-1, :) * 1E5);
            dZdp = (Z(2:end, :) - Z(1:end-1, :)) ./(P(2:end, :) * 1E5 - P(1:end-1, :) * 1E5);

            H2.dKdy1  = ( H2.K_mass(:, 2:end) -  H2.K_mass(:, 1:end-1)) ./( Y1(:, 2:end) - Y1(:, 1:end-1));
            CH4.dKdy1 = (CH4.K_mass(:, 2:end) - CH4.K_mass(:, 1:end-1)) ./( Y1(:, 2:end) - Y1(:, 1:end-1));
            H2O.dKdy1 = (H2O.K_mass(:, 2:end) - H2O.K_mass(:, 1:end-1)) ./( Y1(:, 2:end) - Y1(:, 1:end-1));
            dZdy1 = (Z(:, 2:end) - Z(:, 1:end-1)) ./(Y1(:, 2:end) - Y1(:, 1:end-1));

            % Creating interpolations
            input_1 = P * 1E5;
            input_2 = y1_;
             H2.K_mass = griddedInterpolant(input_1, input_2,  H2.K_mass);
            CH4.K_mass = griddedInterpolant(input_1, input_2, CH4.K_mass);
            H2O.K_mass = griddedInterpolant(input_1, input_2, H2O.K_mass);
            obj.zcf = griddedInterpolant(input_1, input_2, Z);

            input_1 = (P(1:end-1, :) + (P(2:end, :))) / 2 * 1E5; % mid-point of P
            input_2 = y1_(2:end, :);
            H2.dKdp  = griddedInterpolant(input_1, input_2,  H2.dKdp);
            CH4.dKdp = griddedInterpolant(input_1, input_2, CH4.dKdp);
            H2O.dKdp = griddedInterpolant(input_1, input_2, H2O.dKdp);
            obj.dzcfdp = griddedInterpolant(input_1, input_2, dZdp);

            input_1 = P(:, 2:end)* 1E5;
            input_2 = (y1_(:, 1:end-1) + (y1_(:, 2:end))) / 2; % mid-point y1
            H2.dKdy1  = griddedInterpolant(input_1, input_2, H2.dKdy1); 
            CH4.dKdy1 = griddedInterpolant(input_1, input_2, CH4.dKdy1);
            H2O.dKdy1 = griddedInterpolant(input_1, input_2, H2O.dKdy1);
            obj.dzcfdy1 = griddedInterpolant(input_1, input_2, dZdy1);


            obj.comp = {H2, CH4, H2O};


            input_1 = P * 1E5;
            input_2 = y1_;
            obj.mu = griddedInterpolant(input_1, input_2, mu_mix);

            input_1 = (P(1:end-1, :) + (P(2:end, :))) / 2 * 1E5; % mid-point of P
            input_2 = y1_(2:end, :);
            dmudp_mix = ( mu_mix(2:end, :) -  mu_mix(1:end-1, :)) ./( P(2:end, :) * 1E5 - P(1:end-1, :) * 1E5);
            obj.dmudp  = griddedInterpolant(input_1, input_2,  dmudp_mix);


            input_1 = P(:, 2:end)* 1E5;
            input_2 = (y1_(:, 1:end-1) + (y1_(:, 2:end))) / 2; % mid-point y1
            dmudy1_mix = ( mu_mix(:, 2:end) -  mu_mix(:, 1:end-1)) ./( Y1(:, 2:end) - Y1(:, 1:end-1));
            obj.dmudy1  = griddedInterpolant(input_1, input_2,  dmudy1_mix);

            fprintf('##### Thermodynamics of H2-CH4-brine done!\n');

        end

    end
end
