% Formulation base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef formulation < handle
    properties
        MatrixAssembler
        NofPhases
        UpWind
        U
        U_T
        U_D
        Tph
        Gph
        Dph % Diffusion
        Mob
        S
        dMob
        dPc
        drhodp
        GravityModel
        SinglePhase
    end
    methods (Abstract)
        %obj = BuildResidual(obj)
        %obj = BuildJacobian(obj)
        obj = UpdateState(obj)
    end
    methods
        function obj = formulation()
            obj.UpWind = cell(1);  %struct('x',{},'y',{},'z',{});
            obj.U =  cell(1);%struct('x',{},'y',{},'z',{});
        end
        function SavePhaseState(obj)
            % virtual call
        end
        function Reset(obj)
            % virtual call
        end
        function UpWindAndPhaseRockFluxes(obj, DiscretizationModel, Phases, ProductionSystem)
            P_Global = zeros(DiscretizationModel.N, obj.NofPhases);
            Depth_Global = zeros(DiscretizationModel.N,1);
            %% Reservoir upwind
            Grid = DiscretizationModel.ReservoirGrid;
            obj.GravityModel.ComputeInterfaceDensities(Grid, ProductionSystem.Reservoir.State, 0);
            % Compute phase rock velocities and Upwind operators
            for i=1:obj.NofPhases
                P = ProductionSystem.Reservoir.State.Properties(['P_', num2str(i)]).Value;
                P_Global(1:Grid.N, i) = P;
                [obj.UpWind{i, 1}, obj.U{i, 1}] = Phases(i).UpWindAndRockFluxes(Grid, P, obj.GravityModel.RhoInt{i, 1});
            end
            for j = 1: obj.NofComponents
                for i=1:obj.NofPhases
                    xph = ProductionSystem.Reservoir.State.Properties(['x_', num2str(j),'ph', num2str(i)]).Value;
                    obj.U_D{2*(j-1)+i, 1} = Phases(i).UpWind_diffusion(Grid, xph);   % add diffusion flux
                end
            end
            if isa(obj,'Overall_Composition_Energy_formulation')
                obj.U_T{1,1} = Phases(1).UpWind_T(Grid, ProductionSystem.Reservoir.State.Properties('T').Value);
            end
            Depth_Global(1:Grid.N) = DiscretizationModel.ReservoirGrid.Depth(:);    % matrix part
            % Fractures upwind operators
            if ProductionSystem.FracturesNetwork.Active
                End = Grid.N;
                Grids = DiscretizationModel.FracturesGrid.Grids;
                for f=1:ProductionSystem.FracturesNetwork.NumOfFrac
                    obj.GravityModel.ComputeInterfaceDensities(Grids(f), ProductionSystem.FracturesNetwork.Fractures(f).State, f);
                    % Compute phase rock velocities and Upwind operators
                    Start = End+1;
                    End = Start + Grids(f).N - 1;
                    for i=1:obj.NofPhases
                        P = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['P_', num2str(i)]).Value;
                        P_Global(Start:End, i) = P;
                        [obj.UpWind{i, f+1}, obj.U{i, f+1}] = Phases(i).UpWindAndRockFluxes(Grids(f), P, obj.GravityModel.RhoInt{i, f+1});
                    end
                    for j = 1: obj.NofComponents
                        for i=1:obj.NofPhases
                            xph = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['x_', num2str(j),'ph', num2str(i)]).Value;
                            obj.U_D{2*(j-1)+i, 1+f} = Phases(i).UpWind_diffusion(Grids(f), xph);   % add diffusion flux
                        end
                    end
                    if isa(obj,'Overall_Composition_Energy_formulation')
                        obj.U_T{1,1+f} = Phases(1).UpWind_T(Grids(f), ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('T').Value);
                    end
                    Depth_Global(Start:End) = Grids(f, 1).Depth;
                end
                for CrossConn = 1:length(DiscretizationModel.CrossConnections)
                    %loop over all cross connections
                    I = DiscretizationModel.ReservoirGrid.N + CrossConn;
                    DiscretizationModel.CrossConnections(CrossConn).CrossConnectionsUpWind(P_Global, Depth_Global, obj.rho, obj.GravityModel.g,Grid.N, I);
                end
            end
        end
    end
end