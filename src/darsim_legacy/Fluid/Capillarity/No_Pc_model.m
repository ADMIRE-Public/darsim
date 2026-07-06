% No Pc model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 8 September 2016
%Last modified: 8 September 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef  No_Pc_model < capillary_pressure_model
    properties
    end
    methods 
        function Initialise(obj, ProductionSystem)
        end
        function Pc = ComputePc(obj, S)
            %Compute Pc and dPc analytically
            Pc = zeros(length(S), 1); 
        end
        function dPc = dPcdS(obj, S)
            dPc = zeros(length(S), 1); 
        end
        function d2Pc = dPcdSdS(obj, S)
            d2Pc = zeros(length(S), 1); 
        end
        
        function dPc = dPcdSnoHysteresis(obj, RelPermModel, State)
            S = State.Properties('S_1').Value; 
            dPc = zeros(length(S), 1); 
        end 
        
        function ComputePcnoHysteresis(obj, RelPermModel, State)
            status_P_1 = State.Properties('P_1');
            status_P_2 = State.Properties('P_2');
            
            % update pressure;
            status_P_1.Value = status_P_2.Value; 
        end  
    end
end
