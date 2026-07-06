% Component class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Mengjie Zhao
%TU Delft
%Created: 15 Nov 2024
%Last modified: 22 July 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef component < matlab.mixin.Heterogeneous & handle
    properties
        rho
        Pc
        Tc
        w
        Mw
        Vc
    end
    methods
        function obj = component()
        end
    end
end