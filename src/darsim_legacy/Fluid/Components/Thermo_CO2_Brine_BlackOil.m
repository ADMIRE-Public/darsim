% Thermodynamics of H2-Brine for BlackOil Fluid Model class
%--------------------------------------------------------------------------
% DARSim 2 Reservoir Simulator
% Author: Yuhang Wang
% TU Delft
% Modified: 06 May 2025
%--------------------------------------------------------------------------
classdef Thermo_CO2_Brine_BlackOil < matlab.mixin.Heterogeneous & handle
    properties
        % gridded interpolant object for data
        % Rs and zcf are pressure dependent
        Rs_p_table;
        dRs_p_table;
        zcf_p_table; % attn. zcf: compressibility factor; z: mole(mass) fraction of a component;
        dzcf_p_table;
    end
    methods
        function obj = Thermo_CO2_Brine_BlackOil(T, Salinity)
            obj.ComputeLookupTable(T, Salinity);
        end
        % Calculate Rs, dRs, zcf, dzcf table;
        function ComputeLookupTable(obj, T, Salinity)
            %% compute lookp table for compressibility
            % compressibility
            syms V;
            R = 83.1447; % bar cm^3 mol^(-1) K^(-1);
            T_C = T - 273.15;
            a_CO2 = 7.54e7 - 4.13e4*T; % bar cm^6 K^(0.5) mol^(-2);
            b_CO2 = 27.8; % cm^3 mol^(-1);

            m_Na = Salinity/1e3/58.44; % mol kg^(-1);
            m_Cl = Salinity/1e3/58.44; % mol kg^(-1);
            % m_Na = 1; % mol kg^(-1);
            % m_Cl = 1; % mol kg^(-1);

            p = (2:2:400)'; % bar;
            Z = zeros(length(p),1);
            Vm = zeros(length(p),1);

            for i = 1:length(p)
                S = vpasolve(V^3-V^2*(R*T/p(i))-V*(R*T*b_CO2/p(i) - a_CO2/(p(i)*T^0.5) + b_CO2^2) - (a_CO2*b_CO2/(p(i)*T^0.5)) == 0, V);
                Vm(i) = max(double(S(imag(S)==0)));

                Z(i) = p(i)*Vm(i)/(R*T);
            end
            dzcf_value = [(Z(2)-Z(1))/(p(2)*10^5-p(1)*10^5); (Z(2:end)-Z(1:end-1))./(p(2:end)*10^5-p(1:end-1)*10^5)]; % convert bar to pa

            obj.zcf_p_table = griddedInterpolant(p*10^5,Z); % convert bar to pa
            obj.dzcf_p_table = griddedInterpolant(p*10^5,dzcf_value); % convert bar to pa

            % fugacity
            a_H2O_CO2 = 7.89e7; % bar cm^6 K^(0.5) mol^(-2);
            b_H2O = 18.18; % cm^3 mol^(-1);

            phi_CO2 = exp(log(Vm./(Vm-b_CO2)) + (b_CO2./(Vm-b_CO2)) - 2*(a_CO2./(R*T^(1.5)*b_CO2))*log((Vm+b_CO2)./Vm) + (a_CO2*b_CO2./(R*T^(1.5)*b_CO2^2))*(log((Vm+b_CO2)./Vm)-(b_CO2./(Vm+b_CO2))) - log((p.*Vm)./(R*T)));
            phi_H2O = exp(log(Vm./(Vm-b_CO2)) + (b_H2O./(Vm-b_CO2)) - 2*(a_H2O_CO2./(R*T^(1.5)*b_CO2))*log((Vm+b_CO2)./Vm) + (a_CO2*b_H2O./(R*T^(1.5)*b_CO2^2))*(log((Vm+b_CO2)./Vm)-(b_CO2./(Vm+b_CO2))) - log((p.*Vm)./(R*T)));

            % equilibrium constants
            K_CO2_0_g = 10.^(1.189 + 1.304e-2*T_C - 5.446e-5*(T_C)^2);
            K_H2O_0 = 10.^(-2.209 + 3.097e-2*T_C - 1.098e-4*(T_C)^2 + 2.048e-7*(T_C)^3);

            % CO2 molality in pure water
            p_ref = 1; % bar;

            Vm_CO2_g_avg = 32.6; % cm^3 mol^(-1);
            Vm_H2O_avg = 18.1; % cm^3 mol^(-1);

            A = K_H2O_0./(phi_H2O.*p).*exp((p-p_ref)*Vm_H2O_avg./(R*T));
            B = phi_CO2.*p./(55.508*K_CO2_0_g).*exp(-(p-p_ref)*Vm_CO2_g_avg./(R*T));

            y_H2O = (1-B)./(1./A-B);
            x_CO2 = B.*(1-y_H2O);

            m_CO2_0 = 55.508.*x_CO2./(1-x_CO2);

            % CO2 activity coefficient
            par = @(c_1,c_2,c_3,c_4,c_5,c_6,T,p) c_1 + c_2*T + c_3/T + c_4*p./T + c_5*p./(630-T) + c_6*T.*log(p); 

            gamma_CO2 = exp(2*m_Na*par(-0.411370585,6.07632013e-4,97.5347708,-0.0237622469,0.0170656236,1.41335834e-5,T,p)...
                        + m_Na*m_Cl*par(3.36389723e-4,-1.98298980e-5,0,2.122220830e-3,-5.24873303e-3,0,T,p));   

            % CO2 molality in brines
            m_CO2 = m_CO2_0./gamma_CO2;
            x_CO2 = m_CO2./(m_CO2 + 55.508 + 2*m_Na);

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

            rho_brine = ((a_1 - p_f*a_2 - p_f.^2*a_3 + a_4*S + a_5*S^2 - p_f*a_6*S - p_f*a_7*S^2 - 0.5*p_f.^2*a_8*S)*1.0e-3).^(-1);

            % brine molar density at standard condition
            n_NaCl = Salinity/1e3/58.44;
            rho_brine_molar_STC = rho_brine/((58.44*(n_NaCl/(n_NaCl*2+55.508)) + 18.052*(55.508/(n_NaCl*2+55.508)))*1.0e-3);

            % solution gas-oil ratio;
            Rs_value = rho_brine_molar_STC*x_CO2./(44.990*(1-x_CO2));
            dRs_value = [(Rs_value(2)-Rs_value(1))/(p(2)*10^5-p(1)*10^5); (Rs_value(2:end)-Rs_value(1:end-1))./(p(2:end)*10^5-p(1:end-1)*10^5)]; % convert bar to pa

            obj.Rs_p_table = griddedInterpolant(p*10^5,Rs_value); % convert bar to pa
            obj.dRs_p_table = griddedInterpolant(p*10^5,dRs_value); % convert bar to pa
        end
    end
end
