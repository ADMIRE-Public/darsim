% Reservoir 
%--------------------------------------------------------------------------
% DARSim 2 Reservoir Simulator
% Author: Matteo Cusini
% TU Delft
% Modified by Túlio de M. Cavalcante (UFPE) in 08 August 2022
%--------------------------------------------------------------------------
classdef reservoir < handle
    properties
        Length
        Width
        Thickness
        Temp
        K
        KRotAngle
        K_coarse
        Por0 % initial porosity
        Por % current porosity
        DPor
        TotalPV
        State
        level
        State_old
        State_old_old
        MaxLevel
        P0
        Cr % rock compressibility
        K_Cond_rock % rock conductivity
        Diff        % diffusion constants
        Disp        % dispersivity
        Cpr % rock specific heat
        Rho % density of reservoir rock
        Dp = 1e-3;% grain diameter
        FaciesIndex
    end
    methods
        function obj = reservoir(length, width, thickness, temperature)
            obj.Length = length;
            obj.Width = width;
            obj.Thickness = thickness;
            obj.Temp = temperature;
            obj.State = status();
            obj.State_old = status();
            obj.State_old_old = status();
        end
        function ComputePorosity(obj, P)
            obj.Por = obj.Por0.*exp(obj.Cr.*(P-obj.P0));
        end
        function ComputeDerPorosity(obj, P)
            obj.DPor = obj.Cr.*obj.Por0.*exp(obj.Cr.*(P-obj.P0));
        end
        function AddPermeabilityPorosity(obj, k, por0, theta)
            obj.K = k;
            obj.KRotAngle = theta;
            obj.Por0 = por0;
            obj.Por = por0;
            obj.TotalPV = obj.Length * obj.Width * obj.Thickness * obj.Por;
        end
        function AddCoarsePermeability(obj, k_coarse)
            obj.K_coarse = k_coarse;
        end
    end
end