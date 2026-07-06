 % Table relative permeability model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: 
%TU Delft
%Created: 29 Jun 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef relperm_model_fractures < relperm_model
    properties
        srnw
        srw
        kr
        dkr
    end
    methods
        function obj = relperm_model_fractures(srw)
            % obj.srnw = obj.srnw;
            obj.srw  = srw;
        end
        function kr = ComputeRelPerm(obj, Phases, Medium, RelPermMatrix)
            status_S_1 = Medium.State.Properties('S_1');
            S = status_S_1.Value;
            % safe-guard
            % S = max(S, obj.srnw);
            S = min(S, 1-obj.srw);

            % phase 1 (gas phase) 
            kr(:,1) = S/(1-obj.srw);
            % phase 2 (oil phase) 
            kr(:,2) = 1.0 - S/(1-obj.srw);           
  
            status_S_1.Value = S; 
            status_S_2 = Medium.State.Properties('S_2');
            status_S_2.Value = 1.0-S;
            
        end
        
        function dkr = ComputeDerivative(obj, Phases, Medium, RelPermMatrix)
            status_S_1 = Medium.State.Properties('S_1');
            S = status_S_1.Value;            % safe-guard
            % S = max(S, obj.srnw);
            S = min(S, 1-obj.srw);
            
            % phase 1 (gas phase)
            dkr(:,1) = zeros(length(S),1) + 1/(1-obj.srw);
            % phase 2 (oil phase) 
            dkr(:,2) = zeros(length(S),1) - 1/(1-obj.srw);
        end 
    end

end