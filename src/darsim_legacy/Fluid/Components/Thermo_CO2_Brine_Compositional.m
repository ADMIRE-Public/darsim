% Thermodynamics of CO2-Brine for Compositional Fluid Model class
%--------------------------------------------------------------------------
% DARSim 2 Reservoir Simulator
% Author: Mengjie Zhao
% TU Delft
% Modified: 06 May 2025
%--------------------------------------------------------------------------
classdef Thermo_CO2_Brine_Compositional < matlab.mixin.Heterogeneous & handle
    properties
        % gridded interpolant object for data
        X_comp1_p_table
        Y_comp1_p_table
        X_comp2_p_table
        Y_comp2_p_table

        dX_comp1_p_table
        dY_comp1_p_table
        dX_comp2_p_table
        dY_comp2_p_table

        zcf_p_table
        dzcf_p_table
    end
    methods
        function obj = Thermo_CO2_Brine_Compositional(T, Salinity)
            obj.ComputeLookupTable(T, Salinity);
        end
        function ComputeLookupTable(obj, T, Salinity)

            R = 83.1446261815324; % bar cm^3 mol^(-1) K^(-1);
            T_C = T - 273.15; % temperature in C;

            m_Na = Salinity/1e3/58.44; % mol kg^(-1);
            m_Cl = Salinity/1e3/58.44; % mol kg^(-1);
            % m_Na = 1; % mol kg^(-1);
            % m_Cl = 1; % mol kg^(-1);

            % Pressure range in bar
            P = linspace(1e5, 4e7, 400)'; % Pressure range from 1e5 Pa to 3e7 Pa in Pascal
            p = P / 1e5; % Convert pressure to bar for some calculations

            % Redlich-Kwong parameters
            a_CO2 = 7.54e7 - 4.13e4 * T; 
            a_CO2_H2O = 7.89e7; % bar cm^6 K^(0.5) mol^(-2);
            b_CO2 = 27.8;       % cm^3 mol^(-1);
            b_H2O = 18.18;

            % CO2 equilibrium constant
            logk0 = 1.189 + 1.304e-2 * T_C - 5.446e-5 * T_C^2;
            k0_CO2 = 10.^logk0;

            % H2O equilibrium constant
            logk0 = -2.209 + 3.097e-2 * T_C - 1.098e-4 * T_C^2 + 2.048e-7 * T_C^3;
            k0_H2O = 10.^logk0;

            max_iter = 200;
            tol = 1e-5;

            y_CO2 = 0.5;
            y_H2O = 1 - y_CO2;
            for i = 1 : length(p)
                P_bar = p(i);
                converged = false;
                for iter = 1:max_iter
                    % Mixture parameters
                    a_mix = y_CO2.^2 * a_CO2 + 2 * y_CO2 .* y_H2O .* a_CO2_H2O;
                    b_mix = y_CO2 .* b_CO2 + y_H2O .* b_H2O;
                    a_y_CO2 = y_H2O * a_CO2_H2O + y_CO2 * a_CO2;
                    a_y_H2O = y_CO2 * a_CO2_H2O;

                    % Solve cubic EOS for Z
                    RTp = R * T ./ P_bar;
                    A = a_mix ./ (P_bar * sqrt(T));
                    B = b_mix;

                    e2 = -RTp;
                    e1 = -(RTp .* B - A + B.^2);
                    e0 = - A .* B;

                    coeffs = [1, e2, e1, e0];
                    rootsZ = roots(coeffs);
                    realRoots = rootsZ(imag(rootsZ) == 0); % Keep only real roots
                    V = max(realRoots); % Assuming vapor phase dominant

                    % Fugacity coefficients
                    ln_phi_CO2 = log(V ./ (V - b_mix)) ...
                        + b_CO2 ./ (V - b_mix) ...
                        - 2 * a_y_CO2 ./ (R * T^1.5 * b_mix) .* log((V + b_mix) ./ V) ...
                        + a_mix .* b_CO2 ./ (R * T^1.5 .* b_mix.^2) ...
                        .* (log((V + b_mix) ./ V) - b_mix ./ (V + b_mix)) ...
                        - log(P_bar .* V ./ (R*T));
                    phi_CO2 = exp(ln_phi_CO2);

                    ln_phi_H2O = log(V ./ (V - b_mix)) ...
                        + b_H2O ./ (V - b_mix) ...
                        - 2 * a_y_H2O ./ (R * T^1.5 .* b_mix) .* log((V + b_mix) ./ V) ...
                        + a_mix .* b_H2O ./ (R * T^1.5 .* b_mix.^2) ...
                        .* (log((V + b_mix) ./ V) - b_mix ./ (V + b_mix)) ...
                        - log(P_bar .* V ./ (R*T));
                    phi_H2O = exp(ln_phi_H2O);

                    deltaP = P_bar - 1; % Pressure difference in bar
                    v_av_H2O = 18.1; % Average molar volume of H2O [cm3/mol]
                    v_av_CO2 = 32.6; % Average molar volume of CO2 [cm3/mol]

                    A = k0_H2O ./ (phi_H2O .* P_bar) .* exp(deltaP .* v_av_H2O ./ (R * T));
                    B = phi_CO2 .* P_bar ./ (55.508 * k0_CO2) .* exp(-deltaP .* v_av_CO2 ./ (R * T));

                    % CO2 activity coefficient
                    par = @(c_1,c_2,c_3,c_4,c_5,c_6,T,p) c_1 + c_2*T + c_3/T + c_4*p./T + c_5*p./(630-T) + c_6*T.*log(p);

                    gamma_CO2 = exp(2*m_Na*par(-0.411370585,6.07632013e-4,97.5347708,-0.0237622469,0.0170656236,1.41335834e-5,T,P_bar)...
                        + m_Na*m_Cl*par(3.36389723e-4,-1.98298980e-5,0,2.122220830e-3,-5.24873303e-3,0,T,P_bar));

                    % salting out
                    B_eq      = B ./ gamma_CO2;
                    y_H2O_new = (1 - B_eq) ./ (1 ./ A - B_eq);
                    y_CO2_new = 1 - y_H2O_new;

                    x_CO2     = B_eq .* (1 - y_H2O_new);
                    x_H2O     = 1 - x_CO2;

                    diff_CO2 = y_CO2_new - y_CO2;

                    if abs(diff_CO2) < tol
                        converged = true;
                        y_H2O = y_H2O_new;
                        y_CO2 = y_CO2_new;
                        break;
                    else
                        y_H2O = y_H2O_new;
                        y_CO2 = y_CO2_new;
                    end
                end

                if converged
                    a_mix = y_CO2.^2 * a_CO2 + 2 * y_CO2 .* y_H2O .* a_CO2_H2O;
                    b_mix = y_CO2 .* b_CO2 + y_H2O .* b_H2O;

                    RTp = R * T ./ P_bar;
                    A = a_mix ./ (P_bar * sqrt(T));
                    B = b_mix;

                    e2 = -RTp;
                    e1 = -(RTp .* B - A + B.^2);
                    e0 = - A .* B;

                    coeffs = [1, e2, e1, e0];
                    rootsZ = roots(coeffs);
                    realRoots = rootsZ(imag(rootsZ) == 0); 
                    V = max(realRoots); 
                    Z(i,1) = P_bar * V / (R * T);

                    Y_CO2(i,1) = y_CO2 * 44.01 ./ (y_CO2 * 44.01 + y_H2O * 18.016);
                    X_CO2(i,1) = x_CO2 * 44.01 ./ (x_CO2 * 44.01 + x_H2O * 18.016);

                    Y_H2O(i,1) = y_H2O * 18.016 ./ (y_CO2 * 44.01 + y_H2O * 18.016);
                    X_H2O(i,1) = x_H2O * 18.016 ./ (x_CO2 * 44.01 + x_H2O * 18.016);

                    K_CO2(i,1) = Y_CO2(i,1) ./ X_CO2(i,1);
                    K_H2O(i,1) = Y_H2O(i,1) ./ X_H2O(i,1);
                end
            end

            obj.X_comp1_p_table = griddedInterpolant(p*10^5,X_CO2);
            obj.Y_comp1_p_table = griddedInterpolant(p*10^5,Y_CO2);
            obj.X_comp2_p_table = griddedInterpolant(p*10^5,X_H2O);
            obj.Y_comp2_p_table = griddedInterpolant(p*10^5,Y_H2O);

            dXCO2_value = [(X_CO2(2)-X_CO2(1))/(p(2)*10^5-p(1)*10^5); (X_CO2(2:end)-X_CO2(1:end-1))./(p(2:end)*10^5-p(1:end-1)*10^5)]; % convert bar to pa
            dXH2O_value = [(X_H2O(2)-X_H2O(1))/(p(2)*10^5-p(1)*10^5); (X_H2O(2:end)-X_H2O(1:end-1))./(p(2:end)*10^5-p(1:end-1)*10^5)]; % convert bar to pa
            obj.dX_comp1_p_table = griddedInterpolant(p*10^5, dXCO2_value);
            obj.dX_comp2_p_table = griddedInterpolant(p*10^5, dXH2O_value);
            dYCO2_value = [(Y_CO2(2)-Y_CO2(1))/(p(2)*10^5-p(1)*10^5); (Y_CO2(2:end)-Y_CO2(1:end-1))./(p(2:end)*10^5-p(1:end-1)*10^5)]; % convert bar to pa
            dYH2O_value = [(Y_H2O(2)-Y_H2O(1))/(p(2)*10^5-p(1)*10^5); (Y_H2O(2:end)-Y_H2O(1:end-1))./(p(2:end)*10^5-p(1:end-1)*10^5)]; % convert bar to pa
            obj.dY_comp1_p_table = griddedInterpolant(p*10^5, dYCO2_value);
            obj.dY_comp2_p_table = griddedInterpolant(p*10^5, dYH2O_value);

            dzcf_value = [(Z(2)-Z(1))/(p(2)*10^5-p(1)*10^5); (Z(2:end)-Z(1:end-1))./(p(2:end)*10^5-p(1:end-1)*10^5)]; % convert bar to pa

            obj.zcf_p_table = griddedInterpolant(p*10^5,Z); % convert bar to pa
            obj.dzcf_p_table = griddedInterpolant(p*10^5,dzcf_value); % convert bar to pa
        end
    end
end
