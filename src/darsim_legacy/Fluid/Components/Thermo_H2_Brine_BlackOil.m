% Thermodynamics of H2-Brine for BlackOil Fluid Model class
%--------------------------------------------------------------------------
% DARSim 2 Reservoir Simulator
% Author: Amin MB, Willemijn van Roijen
% TU Delft
% Created: 06 May 2025
%--------------------------------------------------------------------------
classdef Thermo_H2_Brine_BlackOil < matlab.mixin.Heterogeneous & handle
    properties
        % gridded interpolant object for data
        % Rs and zcf are pressure dependent
        Rs_p_table;
        dRs_p_table;
        zcf_p_table; % attn. zcf: compressibility factor; z: mole(mass) fraction of a component;
        dzcf_p_table;
    end
    methods
        function obj = Thermo_H2_Brine_BlackOil(T, Salinity)
            obj.ComputeLookupTable(T, Salinity);
        end
        % Calculate Rs, dRs, zcf, dzcf table;
        function ComputeLookupTable(obj, T, Salinity)
            % Thermodynamics of Hydrogen-brine,
            % Formulation drived by Willemijn
            % Implemented into Artur's version of the code by Amin
            
            syms V;
            R = 83.1447; % bar cm^3 mol^(-1) K^(-1);

            % ===== H2 input =====
            % Hydrogen initial temperature in [C])
            T_C = T - 273.15;
            % Critical pressure of H2 [bar]
            P_C_H2 = 13.130;
            % Critical Temperature H2 [K]
            T_C_H2 = 33.18;
            % Accentric factor of H2
            w_H2 = -0.2150;
            % Molar volume in the critical point [cm3/mol]
            Vc_H2 = 64.4828;
            % Molecular weight [g/mol]
            Mw_H2 = 2.01588;
            % Density at standard conditions H2 [kg/m3]
            Rho_H2 = 0.0898;
            
            k_H2 = 0.37464 + 1.4522 * w_H2 - 0.26992 * w_H2^2; % Pitzer's acentric factor correlation
            rho_molar_H2 = Rho_H2 * (1/Mw_H2) ; % Molar density of H2 at standard conditions [mol/L])
           
            % ===== H2O input =====
            % Critical pressure (Critical for H2 [bar])
            P_C_H2O = 220.640;
            % Critical temperature (Critical Temperature H2 [K])
            T_C_H2O = 647.096;
            % Accentric factor of H2O
            w_H2O = 0.3443;
            % Molecular weight [g/mol]
            Mw_H2O = 18.01528;
            % Density at standard conditions H2 [kg/m3]
            Rho_H2O= 1000;
            % Average molecular molume [cm3/mol]
            Vm_H2O_avg = 18.1;
            
            k_H2O = 0.37464 + 1.5422 * w_H2O - 0.26992 * w_H2O^2; % Part of Eq. 20
            rho_molar_H2O = Rho_H2O * (1/Mw_H2O) ; % Molar density of H2 at standard conditions [mol/L])
            
            % ===== Salinity =====
            m_Na = Salinity/1e3/58.44; % mol kg^(-1);
            m_Cl = Salinity/1e3/58.44; % mol kg^(-1);            
            
            % Pressure and Temperature
            p = (2:2:400)'; % Pressure range for the lookuptable bar

            %% Thermodynamic parameters using PR - EoS see Ganji et al. 2012
            alpha_H2_T = (1 + (k_H2) * (1 - (T / T_C_H2)^0.5))^2; %Eq. 20
            a_H2_T = 0.45724 * R^2 * T_C_H2^2 * alpha_H2_T / P_C_H2; %Eq. 18
            b_H2 = 0.07780 * R * T_C_H2 / P_C_H2; %Eq. 19
            A = a_H2_T * p / (R * T)^2; %Eq. 16
            B = b_H2 * p / (R * T); % Eq. 17

            alpha_H2O_T = (1 + k_H2O * (1 - (T ./ T_C_H2O).^0.5)).^2; %Eq. 20
            a_H2O_T = 0.45724 * R^2 * T_C_H2O^2 * alpha_H2O_T ./ P_C_H2O; %Eq. 18
            b_H2O = 0.07780 * R * T_C_H2O / P_C_H2O; %Eq. 19   
            A_H2O = a_H2O_T .* p ./ (R * T).^2; %Eq. 16
            B_H2O = b_H2O .* p ./ (R * T); % Eq. 17

            delta1 = 1 + sqrt(2); %Defined in PR EoS
            delta2 = 1 - sqrt(2); %Defined in PR EoS
            
            %% Compressibility
            P_app = p - 1.01325 ; % (applied pressure at P_app = 0 bars, the absolute pressure is 1 atm or 1.01325 bars)
            Z = zeros(length(p),1);
            Vm = zeros(length(p),1);

            for i = 1:length(p)
                S = vpasolve(V^3 - (1 - B(i)) * V^2 + (A(i) - 2*B(i) - 3*B(i).^2) * V - (A(i).*B(i) - B(i).^2 - B(i).^3) == 0, V); % Eq. 15
                Z_h = max(double(S(imag(S)==0)));
                Z_l = min(double(S(S>0)));
                deltaG = (Z_h - Z_l) + log((Z_l-B(i))/(Z_h-B(i))) + A(i)/B(i)/(delta2-delta1)*log((Z_l+delta1*B(i))*(Z_h+delta2*B(i))/(Z_l+delta2*B(i))*(Z_h+delta1*B(i))); %Eq. 23
                if deltaG > 0
                    Z(i) = Z_l;
                else
                    Z(i) = Z_h;          
                end
                Vm(i) = Z(i) * R * T / p(i); %Molecular Volume
            end
            
            dzcf_value = [(Z(2)-Z(1))/(p(2)*10^5-p(1)*10^5); (Z(2:end)-Z(1:end-1))./(p(2:end)*10^5-p(1:end-1)*10^5)]; % convert bar to pa
            
            obj.zcf_p_table = griddedInterpolant(p*10^5,Z); % convert bar to pa
            obj.dzcf_p_table = griddedInterpolant(p*10^5,dzcf_value); % convert bar to pa

            %% 2. Fugucity cofficients (Eq22) Prausnitz et al

            phi_H2 = exp((Z-1) - log(Z-B) - A./(2.828.*B).* log((Z+2.414.*B)./(Z-0.414.*B)));  
            phi_H2O = exp(B_H2O ./ B .* (Z-1) - log(Z-B) + A./(2.828.*B).*(B_H2O./B - 2.*(sqrt(a_H2_T.*a_H2O_T).*(1-0.1))./a_H2_T) .* log((Z+2.414.*B)./(Z-0.414.*B)));
            
            %% 3. Henry's constant
            % Density of pure water (Appendix A)
            V0 = (1+18.159725e-3*T_C) ./ (0.9998396 + 18.224944e-3*T_C - 7.922210e-6*T_C.^2 - 55.44846e-9*T_C.^3 + 149.7562e-12*T_C.^4 - 393.2952e-15*T_C.^5);
            B = 19654.320 + 147.037*T_C - 2.21554*T_C.^2 + 1.0478e-2*T_C.^3 - 2.2789e-5*T_C.^4;
            A1 = 3.2891 - 2.3910e-3*T_C + 2.8446e-4*T_C.^2 - 2.8200e-6*T_C.^3 + 8.477e-9*T_C.^4;
            A2 = 6.245e-5 - 3.913e-6*T_C + 3.499e-8*T_C.^2 - 7.942e-10*T_C.^3 + 3.299e-12*T_C.^4;
            V = V0 - (V0.*P_app) ./ (B+A1.*P_app+A2.*P_app.^2); % cm3/g, typo in Ganji2012
            rho_H2O = 1./V; %g/cm3 
            
            % Fugucity of pure water (Appendix A)
            V_molar = 18.0152*V;
            tau = 1 - T/T_C_H2O;
            % vapor pressure of pure water, only applicable up to critical point
            P_s = P_C_H2O .* exp(T_C_H2O./T.*(-7.85951783*tau + 1.84408259*tau.^1.5 - 11.7866497*tau.^3 + 22.6807411*tau.^3.5 - 15.9618719*tau.^4 + 1.80122502*tau.^7.5)); 
            f_H2O = P_s.*exp((p - P_s).*V_molar ./ (R*T)); % Poynting equation
            
            % Henry's constant of H2 %Coefficients from Shabani et al. 
            deltaB_H2 = -4.23407 + 5.41551*sqrt(1000./T); % The difference in interaction between dissimilar solvent molecules
            kH_H2 = exp((1-0.294994)*log(f_H2O) + 0.294994*log(R*T/Mw_H2O.*rho_H2O) + 2*rho_H2O.*deltaB_H2);
    
            
            %% 4. Activity coefficient (Eq26)
            par = @(c_1,c_2,c_3,c_4,c_5,c_6, c_7, c_8, c_9, c_10, T, P) c_1 + c_2*T + c_3./T + c_4*P + c_5./P + c_6*P./T + c_7*T./(P.^2) + c_8*P./(630-T) + c_9*T.*log(P) + c_10*P./(T.^2); 
            
            %Coefficients from Shabani et al. 
            gamma_H2 = exp(2*m_Na*par(-2.1432831, 0.0031411257, 392.20546, -0.0000286012, 0, 0.002352716, 0, -0.0024422, 0.0000029806, -0.389, T, p)...
                        + m_Na*m_Cl*par(-0.0040631, 0, 0, -0.000003665, 0.17004, -0.000418, 0.000588, -0.000391, 0, 0.19913, T, p));   
            
            %% 5. K-value
            K_H2 = kH_H2 .* gamma_H2 ./ (phi_H2.*p);
            K_H2O_0 = 10.^(-2.209 + 3.097e-2*T_C - 1.098e-4*(T_C).^2 + 2.048e-7*(T_C).^3);
            p_ref = 1; % bar;          
            K_H2O = K_H2O_0 .* exp((p-p_ref)*Vm_H2O_avg./(R*T)) ./ (phi_H2O.*p); %eq.11 Ganji
            
            %% Flash calculation
            
            % brine density
            a_1 = 5.916365 - 0.01035794*T + 0.9270048e-5*T^2 - 1127.522/T + 100674.1/(T^2);
            a_2 = 0.5204914e-2 - 0.10482101810e-4*T + 0.8328532e-8*T^2 - 1.1702939/T + 102.2783/(T^2);
            a_3 = 0.118547e-7 - 0.6599143e-10*T;
            a_4 = -2.5166 + 0.0111766*T - 0.170552e-4*T^2;
            a_5 = 2.84851 - 0.0154305*T + 0.223982e-4*T^2;
            a_6 = -0.0014814 + 0.82969e-5*T - 0.12469e-7*T^2;
            a_7 = 0.0027141 - 0.15391e-4 + 0.22655e-7*T^2;
            a_8 = 0.62158e-6 - 0.40075e-8*T + 0.65972e-11*T^2;

            S = Salinity*1.0e-3/(Salinity*1.0e-3 + 1000); % NaCl mass fraction in brine;
            % p_f = 1.01972*p; % pressure in kgf/(cm^2)
            p_f = 1.0332; % standard condition;

            rho_brine = ((a_1 - p_f*a_2 - p_f.^2*a_3 + a_4*S + a_5*S^2 - p_f*a_6*S - p_f*a_7*S^2 - 0.5*p_f.^2*a_8*S)*1.0e-3).^(-1);%[g/L]

            % brine molar density at standard condition
            n_NaCl = Salinity/1e3/58.44; % mol/kg
            rho_brine_molar_STC = rho_brine/(58.44*(n_NaCl*2/(n_NaCl*2+rho_molar_H2O)) + 18.052*(rho_molar_H2O/(n_NaCl*2+rho_molar_H2O)));%[mol/L]

            %% Eq. 29b in Ganji 2014           
            x_H2 = (p .* phi_H2) ./ (kH_H2 .* gamma_H2) .* ...
                   (1 - ((1 - (p .* phi_H2) ./ (kH_H2 .* gamma_H2)) ./ ...
                   ((1 ./ ((K_H2O_0 ./ (phi_H2O.*p)) .* exp((p - p_ref) .* Vm_H2O_avg ./ (R .* T)))) - ...
                   (p .* phi_H2) ./ (kH_H2 .* gamma_H2))));

            % H2 molality in brines
            m_H2_0 = 55.508.*x_H2./(1-x_H2);
            m_H2 = m_H2_0./gamma_H2;
            x_H2 = m_H2./(m_H2 + 55.508 + 2*m_Na);

            % solution gas-oil ratio;
            Rs_value = rho_brine_molar_STC*x_H2./(rho_molar_H2*(1-x_H2));
            dRs_value = [(Rs_value(2)-Rs_value(1))/(p(2)*10^5-p(1)*10^5); (Rs_value(2:end)-Rs_value(1:end-1))./(p(2:end)*10^5-p(1:end-1)*10^5)]; % convert bar to pa
            
            obj.Rs_p_table = griddedInterpolant(p*10^5,Rs_value); % convert bar to pa
            obj.dRs_p_table = griddedInterpolant(p*10^5,dRs_value); % convert bar to pa
        end

    end
end
