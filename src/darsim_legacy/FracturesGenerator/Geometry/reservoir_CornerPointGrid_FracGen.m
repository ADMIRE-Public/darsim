% class reservoir for DARSim2FracGen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini & Mousa HosseiniMehr
%TU Delft
%Created: 2016-07-04
%Last modified: 2017-03-01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef reservoir_CornerPointGrid_FracGen < reservoir_FracGen
    properties
        CornerPointGridData
        pEDFM_alpha_Transmissibility
        % alpha is a the fraction of blockage (from 0 to 1) between
        % neighboring matrix cells due to existence of fracture planes when
        % considering pEDFM method. The alpha=0 means no blockage and
        % alpha=1 means total blockage.
    end
    methods
        function Discretize(obj, CompareAccuracy, Type)
            obj.dx = obj.CornerPointGridData.Cells.dx;
            obj.dy = obj.CornerPointGridData.Cells.dy;
            obj.dz = obj.CornerPointGridData.Cells.dz;
            obj.CellVolume = obj.CornerPointGridData.Cells.Volume;
            obj.Epsilon = CompareAccuracy * min([obj.dx;obj.dy;obj.dz]);
            if strcmp(Type,'pEDFM')
                obj.pEDFM_alpha_Transmissibility = zeros( obj.CornerPointGridData.N_InternalFaces , 1 );
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
            disp( ['Total Grids  : ', num2str(obj.Nx), ' x ',  num2str(obj.Ny), ' x ', num2str(obj.Nz),...
                               ' = ', num2str(obj.Nx * obj.Ny * obj.Nz)] );
            disp( ['Active Grids  : ', num2str(obj.CornerPointGridData.N_ActiveCells)] );
            disp( '---------------------------------------------------------' );
        end
    end
end