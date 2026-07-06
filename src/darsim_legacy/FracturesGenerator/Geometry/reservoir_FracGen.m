% class reservoir for DARSim2FracGen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini & Mousa HosseiniMehr
%TU Delft
%Created: 2016-07-04
%Last modified: 2017-03-01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef reservoir_FracGen < handle
    properties
        Lx
        Ly
        Lz
        Nx
        Ny
        Nz
        N
        dx
        dy
        dz
        ADM
        Epsilon
        OvelapList = [];
        CellVolume
        Stress
    end
    methods
        function obj = reservoir_FracGen(ReservoirPropertiesInput)
            obj.Lx = ReservoirPropertiesInput.Size(1);
            obj.Ly = ReservoirPropertiesInput.Size(2);
            obj.Lz = ReservoirPropertiesInput.Size(3);
            obj.Nx = ReservoirPropertiesInput.Grid.N(1);
            obj.Ny = ReservoirPropertiesInput.Grid.N(2);
            obj.Nz = ReservoirPropertiesInput.Grid.N(3);
            obj.N  = ReservoirPropertiesInput.Grid.N_ActiveCells;
            obj.ADM.Active      = ReservoirPropertiesInput.ADM(1);
            obj.ADM.CoarseLevel = ReservoirPropertiesInput.ADM(2);
            obj.ADM.CoarseRatio = ReservoirPropertiesInput.ADM(3:5);
            obj.Stress = ReservoirPropertiesInput.Stress;
        end
    end
end