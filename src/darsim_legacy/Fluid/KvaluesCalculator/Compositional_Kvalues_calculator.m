% Co Oil K-values 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: 
%TU Delft
%Created: 15 Nov 2024
%Last modified: 15 Nov 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Compositional_Kvalues_calculator < Kvalues_calculator
    properties

    end
    methods
        function k = Compute(obj, Status, Thermo)
            p = Status.Properties('P_2').Value;
            if ismember('T', Status.Properties.keys)    % Thermal
                T = Status.Properties('T').Value;  
                k(:,1) = Thermo.Y_comp1_pT_table(p,T) ./ Thermo.X_comp1_pT_table(p,T);
                k(:,2) = Thermo.Y_comp2_pT_table(p,T) ./ Thermo.X_comp2_pT_table(p,T);  
            else                                        % non-thermal
                k(:,1) = Thermo.Y_comp1_p_table(p) ./Thermo.X_comp1_p_table(p);
                k(:,2) = Thermo.Y_comp2_p_table(p) ./Thermo.X_comp2_p_table(p);
            end
        end
        function dk = DKvalDp(obj, Status, Thermo)
            p = Status.Properties('P_2').Value;
            if ismember('T', Status.Properties.keys)    % Thermal
                T = Status.Properties('T').Value;  
                dk(:,1) = (Thermo.dY_comp1_dP_pT_table(p,T).*Thermo.X_comp1_pT_table(p,T)-Thermo.dX_comp1_dP_pT_table(p,T).*Thermo.Y_comp1_pT_table(p,T))./(Thermo.X_comp1_pT_table(p,T)).^2;
                dk(:,2) = (Thermo.dY_comp2_dP_pT_table(p,T).*Thermo.X_comp2_pT_table(p,T)-Thermo.dX_comp2_dP_pT_table(p,T).*Thermo.Y_comp2_pT_table(p,T))./(Thermo.X_comp2_pT_table(p,T)).^2;
            else
                dk(:,1) = (Thermo.dY_comp1_p_table(p).*Thermo.X_comp1_p_table(p)-Thermo.dX_comp1_p_table(p).*Thermo.Y_comp1_p_table(p))./(Thermo.X_comp1_p_table(p)).^2;
                dk(:,2) = (Thermo.dY_comp2_p_table(p).*Thermo.X_comp2_p_table(p)-Thermo.dX_comp2_p_table(p).*Thermo.Y_comp2_p_table(p))./(Thermo.X_comp2_p_table(p)).^2;
            end
        end
        function dk = DKvalDT(obj, Status, Thermo)
            p = Status.Properties('P_2').Value;
            T = Status.Properties('T').Value;
            dk(:,1) = (Thermo.dY_comp1_dT_pT_table(p,T).*Thermo.X_comp1_pT_table(p,T)-Thermo.dX_comp1_dT_pT_table(p,T).*Thermo.Y_comp1_pT_table(p,T))./(Thermo.X_comp1_pT_table(p,T)).^2;
            dk(:,2) = (Thermo.dY_comp2_dT_pT_table(p,T).*Thermo.X_comp2_pT_table(p,T)-Thermo.dX_comp2_dT_pT_table(p,T).*Thermo.Y_comp2_pT_table(p,T))./(Thermo.X_comp2_pT_table(p,T)).^2;         
        end
    end
end