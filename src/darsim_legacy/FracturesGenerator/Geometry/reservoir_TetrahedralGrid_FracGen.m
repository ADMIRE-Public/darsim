% class reservoir for DARSim2FracGen
%--------------------------------------------------------------------------
% DARSim 2 Reservoir Simulator
% Author: Túlio de M. Cavalcante
% UFPE
% Created: 14 July 2022
%--------------------------------------------------------------------------
classdef reservoir_TetrahedralGrid_FracGen < reservoir_FracGen
    properties
        TetrahedralGridData
        pEDFM_alpha_Transmissibility
        % alpha is a the fraction of blockage (from 0 to 1) between
        % neighboring matrix cells due to existence of fracture planes when
        % considering pEDFM method. The alpha=0 means no blockage and
        % alpha=1 means total blockage.
    end
    methods
        function Discretize(obj, CompareAccuracy, Type)
            obj.CellVolume = obj.TetrahedralGridData.Cell.Volume;
            obj.Epsilon = CompareAccuracy;
            if strcmp(Type,'pEDFM')
                obj.pEDFM_alpha_Transmissibility = zeros( size(obj.TetrahedralGridData.Face.Inner.Vertices,1) , 1 );
            end
        end
        function PrintInfo(obj)
            disp( char(5) );
            if obj.Nz > 1
                disp( '3D Reservoir geometry:');
            else
                disp( '2D Reservoir geometry:');
            end
            disp( ['Length: ', num2str(obj.Lx), ' [m]'] );
            disp( ['Width : ', num2str(obj.Ly), ' [m]'] );
            disp( ['Depth : ', num2str(obj.Lz), ' [m]'] );
            disp( ['Total Grids  : ', num2str(size(obj.TetrahedralGridData.Cell.Vertices,1))] );
            disp( ['Active Grids  : ', num2str(size(obj.TetrahedralGridData.Cell.Vertices,1))] );
            disp( '---------------------------------------------------------' );
        end
    end
end