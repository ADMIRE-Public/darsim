% Compositional K-values 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Amin MB
%TU Delft
%Created: 24 July 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Kvalues_calculator_3component < Kvalues_calculator
    properties

    end
    methods
        function k = Compute(obj, Status, Thermo)
            p = Status.Properties('P_2').Value;
            y1 = Status.Properties('x_1ph1').Value; 
            N = length(p);
            nc = length(Thermo.comp);
            k = zeros(N, nc);
            for i = 1:nc
                k(:,i) = Thermo.comp{i}.K(p, y1);
            end
        end
        function [dkdp, dkdy1] = DKvalDp(obj, Status, Thermo)
            p = Status.Properties('P_2').Value;
            y1 = Status.Properties('x_1ph1').Value; 
            N = length(p);
            nc = length(Thermo.comp);
            dkdp = zeros(N, nc);
            dkdy1 = zeros(N, nc);
            for i = 1:nc
                dkdp(:,i)  = Thermo.comp{i}.dKdp(p, y1);
                dkdy1(:,i) = Thermo.comp{i}.dKdy1(p, y1);
            end
        end
    end
end