% Thermodynamics of CH4-Brine for Compositional Fluid Model class
%--------------------------------------------------------------------------
% DARSim 2 Reservoir Simulator
% Author: Amin MB
% TU Delft
% Created: 09 Feb 2026
%--------------------------------------------------------------------------
classdef Thermo_CH4_Brine_Compositional < matlab.mixin.Heterogeneous & handle
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
        function obj = Thermo_CH4_Brine_Compositional(T, Salinity)
            obj.ComputeLookupTable(T, Salinity);
        end
        % Calculate Rs, dRs, zcf, dzcf table;
        function ComputeLookupTable(obj, T, Salinity)
            % Thermodynamics of Methane-brine

            R = 83.1447; % bar cm^3 mol^(-1) K^(-1);
            T_C = T - 273.15;

            %% ----------- pure‑component data -----------------------------------
            Tc = [190.564; 647.096 ];  % Critical Temperature K   (CH4 , H2O)
            Pc = [45.992; 220.64  ];  % Critical pressure bar
            w  = [0.01142; 0.3443 ];
            kij= [0 0.47893; 0.47893 0]; % Ziabakhsh-Ganji & Kooi 2012

            kap = zeros(2,1);
            kap(1)= 0.37464 + 1.4522*w(1) - 0.26992*w(1)^2;   % CH4
            kap(2)= 0.37464 + 1.5422*w(2) - 0.26992*w(2)^2;   % H2O

            alpha= (1 + kap.*(1 - sqrt(T./Tc))).^2;
            ai    = 0.45724*R^2 .* Tc.^2 .* alpha ./ Pc;
            bi    = 0.07780*R .* Tc ./ Pc;

            % Average molecular molume [cm3/mol]
            Vm_H2O_avg = 18.1;
            K_H2O_0 = 10.^(-2.209 + 0.03097*T_C - 1.098e-4*T_C.^2 + 2.048e-7*T_C.^3);   % Eq.(10)
            p_ref = 1;

            % ===== Salinity =====
            m_Na = Salinity/1e3/58.44; % mol kg^(-1);
            m_Cl = Salinity/1e3/58.44; % mol kg^(-1);

            % Pressure and Temperature
            p = (1:2:400)'; % Pressure range for the lookuptable bar

            % Activity coefficient
            par = @(c_1,c_2,c_3,c_4,c_5,c_6, c_7, c_8, c_9, c_10, T, P) c_1 + c_2*T + c_3./T + c_4*P + c_5./P + c_6*P./T + c_7*T./(P.^2) + c_8*P./(630-T) + c_9*T.*log(P) + c_10*P./(T.^2);

            %Coefficients from Shabani et al.
            gamma_CH4 = exp(2*m_Na*par(-5.7066455E-01, +7.2997588E-04, +1.5176903E+02, +3.1927112E-05, +0.0000000E+00, -1.6426510E-05, +0.0000000E+00, +0.0000000E+00, +0.0000000E+00, +0.0000000E+00, T, p)...
                       + m_Na*m_Cl*par(-2.9990084E-03, +0.0000000E+00, +0.0000000E+00, +0.0000000E+00, +0.0000000E+00, +0.0000000E+00, +0.0000000E+00, +0.0000000E+00, +0.0000000E+00, +0.0000000E+00, T, p));

            Np   = numel(p);
            x_CH4= zeros(Np,1);       % x_CH4 in liquid
            y_CH4= zeros(Np,1);       % y_CH4 in vapour
            Z = zeros(Np,1);
            converged = false;

            for ip = 1:Np
                P   = p(ip);
                % Wilson initial K
                Ki  = Pc./P .* exp(5.373.*(1+w).*(1-Tc./T));
                K_CH4  = Ki(1);   K_H2O = Ki(2);

                for iter = 1:200
                    x = (1-Ki(2)) / (Ki(1)-Ki(2));
                    xvec = [x ; 1-x];
                    yvec = [K_CH4 ; K_H2O] .* xvec;

                    [phiV,~]  = obj.PRphi(T,P,yvec,ai,bi,kij);

                    kH     = obj.Henry_CH4(P,T);
                    K_CH4   = (kH*gamma_CH4(ip))/(phiV(1)*P);
                    K_H2O  = K_H2O_0 .* exp((P - p_ref).*Vm_H2O_avg./(R*T)) ./ (phiV(2) .* P);
                    Ki_new = [K_CH4 ; K_H2O];
                    
                    if max(abs(log(Ki_new./Ki))) < 1e-8
                        converged = true;
                        Ki = Ki_new;
                        x = (1-Ki(2)) / (Ki(1)-Ki(2));
                        xvec = [x;1-x];
                        yvec = Ki .* xvec;
                        break
                    end
                    Ki = Ki_new;
                end

                if converged
                    x_CH4(ip) = x;
                    y_CH4(ip) = yvec(1);
                    [~,Z(ip)]  = obj.PRphi(T,P,yvec,ai,bi,kij);
                end
            end

            dzcf_value = [(Z(2)-Z(1))/(p(2)*10^5-p(1)*10^5); (Z(2:end)-Z(1:end-1))./(p(2:end)*10^5-p(1:end-1)*10^5)]; % convert bar to pa

            obj.zcf_p_table = griddedInterpolant(p*10^5,Z); % convert bar to pa
            obj.dzcf_p_table = griddedInterpolant(p*10^5,dzcf_value); % convert bar to pa

            m_CH4_0 = 55.508.*x_CH4./(1-x_CH4);

            % CH4 molality in brines
            m_CH4 = m_CH4_0./gamma_CH4;
            x_CH4 = m_CH4./(m_CH4 + 55.508 + 2*m_Na);

            Y_CH4 = y_CH4 * 16.0425 ./ (y_CH4 * 16.0425 + (1-y_CH4) * 18.016);
            X_CH4 = m_CH4 * 16.0425 ./ (1000 + m_Na * 58.448 + m_CH4 * 16.0425);

            Y_H2O = (1-y_CH4) * 18.016 ./ (y_CH4 * 16.0425 + (1-y_CH4) * 18.016);
            X_H2O = 1000 ./ (1000 + m_Na * 58.448 + m_CH4 * 16.0425);

            K_CH4 = Y_CH4 ./ X_CH4;
            K_H2O = Y_H2O ./ X_H2O;

            obj.X_comp1_p_table = griddedInterpolant(p*10^5,X_CH4);
            obj.Y_comp1_p_table = griddedInterpolant(p*10^5,Y_CH4);
            obj.X_comp2_p_table = griddedInterpolant(p*10^5,X_H2O);
            obj.Y_comp2_p_table = griddedInterpolant(p*10^5,Y_H2O);

            dXCH4_value = [(X_CH4(2)-X_CH4(1))/(p(2)*10^5-p(1)*10^5); (X_CH4(2:end)-X_CH4(1:end-1))./(p(2:end)*10^5-p(1:end-1)*10^5)]; % convert bar to pa
            dXH2O_value = [(X_H2O(2)-X_H2O(1))/(p(2)*10^5-p(1)*10^5); (X_H2O(2:end)-X_H2O(1:end-1))./(p(2:end)*10^5-p(1:end-1)*10^5)]; % convert bar to pa
            obj.dX_comp1_p_table = griddedInterpolant(p*10^5, dXCH4_value);
            obj.dX_comp2_p_table = griddedInterpolant(p*10^5, dXH2O_value);
            dYCH4_value = [(Y_CH4(2)-Y_CH4(1))/(p(2)*10^5-p(1)*10^5); (Y_CH4(2:end)-Y_CH4(1:end-1))./(p(2:end)*10^5-p(1:end-1)*10^5)]; % convert bar to pa
            dYH2O_value = [(Y_H2O(2)-Y_H2O(1))/(p(2)*10^5-p(1)*10^5); (Y_H2O(2:end)-Y_H2O(1:end-1))./(p(2:end)*10^5-p(1:end-1)*10^5)]; % convert bar to pa
            obj.dY_comp1_p_table = griddedInterpolant(p*10^5, dYCH4_value);
            obj.dY_comp2_p_table = griddedInterpolant(p*10^5, dYH2O_value);

        end

        function [phi,Z] = PRphi(obj,T,P,y,ai0,bi,kij)
            R = 83.1447;
            Aij = sqrt(ai0.*ai0.').*(1-kij);
            A   = y.'*Aij*y;
            Bm = y.'*bi;
            Afac= P*A /(R^2*T^2);
            Bfac= P*Bm/(R*T);
            coef= [1 -(1-Bfac) (Afac-3*Bfac^2-2*Bfac) -(Afac*Bfac-Bfac^2-Bfac^3)];
            % Zr  = roots(coef); Zr = Zr(abs(imag(Zr))<1e-9);
            % Zl=min(Zr(Zr>0)); Zh=max(Zr);
            % DG=(Zh-Zl)+log((Zl-Bfac)/(Zh-Bfac))+Afac/Bfac/2.828*log((Zl+2.414*Bfac)*(Zh-0.414*Bfac)/((Zl-0.414*Bfac)*(Zh+2.414*Bfac)));
            % Z=Zl;
            % if DG<0
            %     Z=Zh;
            % end
            % 1. Get the roots
            Zr = roots(coef);
            % 2. Keep only real and physical roots (must be > Bfac)
            % We use a small buffer (1e-6) to avoid log(0) issues
            Zr = Zr(imag(Zr) == 0 & Zr > Bfac + 1e-6);
            Zr = sort(Zr);
            
            if isempty(Zr)
                % Fallback if no roots are found > Bfac
                Z = Bfac + 1e-6; 
            elseif isscalar(Zr)
                % Single phase region
                Z = Zr;
            else
                % Multiple roots: Calculate DG to find the stable one
                Zl = Zr(1);
                Zh = Zr(end);
                % Gibbs comparison formula
                DG = (Zh - Zl) + ...
                     log((Zl - Bfac) / (Zh - Bfac)) + ...
                     (Afac / (Bfac * 2.828)) * log((Zl + 2.414 * Bfac) * (Zh - 0.414 * Bfac) / ((Zl - 0.414 * Bfac) * (Zh + 2.414 * Bfac)));
                if DG > 0  % Note: check your sign convention here
                    Z = Zl;
                else
                    Z = Zh;
                end
            end
            sqrtai = sqrt(ai0);
            Gamma = (sqrtai.*sqrtai.' .* (1-kij))/A;
            lnphi = (bi/Bm).*(Z-1) - log(Z-Bfac) - Afac/(2.828*Bfac).*(2*(Gamma*y)-bi/Bm).*log((Z+(1+sqrt(2))*Bfac)/(Z+(1-sqrt(2))*Bfac));
            phi = exp(lnphi);
        end

        % ----- Henry constant (Shabani 2017) ---------------------------------
        function kH = Henry_CH4(obj,p,T)
            R = 83.1447;                       % bar·cm3 mol-1 K-1
            Mw_H2O = 18.01528;                 % g/mol
            T_C = T-273.15;                    % °C
            P_app = p - 1.01325;               % bar
            T_c_H2O=647.096; P_c_H2O=220.640;
            % density of water
            V0 = (1+18.159725e-3*T_C)./(0.9998396+18.224944e-3*T_C-7.922210e-6*T_C.^2-55.44846e-9*T_C.^3+149.7562e-12*T_C.^4-393.2952e-15*T_C.^5);
            B  = 19654.320+147.037*T_C-2.21554*T_C.^2+1.0478e-2*T_C.^3-2.2789e-5*T_C.^4;
            A1 = 3.2891-2.3910e-3*T_C+2.8446e-4*T_C.^2-2.8200e-6*T_C.^3+8.477e-9*T_C.^4;
            A2 = 6.245e-5-3.913e-6*T_C+3.499e-8*T_C.^2-7.942e-10*T_C.^3+3.299e-12*T_C.^4;
            V  = V0 - (V0.*P_app)./(B+A1.*P_app+A2.*P_app.^2);   % cm3/g
            rho = 1./V;                                         % g/cm3
            % pure‑water fugacity
            Vmol = 18.01528*V;
            tau  = 1 - T/T_c_H2O;
            Ps   = P_c_H2O .* exp(T_c_H2O./T.*(-7.85951783*tau + 1.84408259*tau.^1.5 - 11.7866497*tau.^3 + 22.6807411*tau.^3.5 - 15.9618719*tau.^4 + 1.80122502*tau.^7.5));
            fw   = Ps .* exp((p-Ps).*Vmol./(R*T));
            % Henry constant
            eta = -0.092248;
            tau = -5.779280;
            beta = 7.262730;
            deltaB = tau + beta * sqrt(1E3 ./ T);
            kH = exp((1 - eta) * log(fw) + eta * log(R * T / Mw_H2O .* rho) + 2 * rho .* deltaB);
        end

    end
end
