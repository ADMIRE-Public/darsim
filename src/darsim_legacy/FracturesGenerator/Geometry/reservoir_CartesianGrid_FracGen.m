% class reservoir for DARSim2FracGen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini & Mousa HosseiniMehr
%TU Delft
%Created: 2016-07-04
%Last modified: 2017-03-01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef reservoir_CartesianGrid_FracGen < reservoir_FracGen
    properties
        Xcm
        Ycm
        Zcm
        Xim
        Yim
        Zim
        Neighbours
        pEDFM_alpha_Tx
        pEDFM_alpha_Ty
        pEDFM_alpha_Tz
        % alpha is a the fraction of blockage (from 0 to 1) between
        % neighboring matrix cells due to existence of fracture planes when
        % considering pEDFM method. The alpha=0 means no blockage and
        % alpha=1 means total blockage.
    end
    methods
        function Discretize(obj, CompareAccuracy, Type)
            obj.dx  = obj.Lx/obj.Nx;
            obj.dy  = obj.Ly/obj.Ny;
            obj.dz  = obj.Lz/obj.Nz;
            obj.CellVolume = obj.dx .* obj.dy .* obj.dz .* ones(obj.N,1);
            obj.Xcm = linspace(obj.dx/2, obj.Lx - obj.dx/2, obj.Nx)';              % Grid centers locations in x-direction
            obj.Ycm = linspace(obj.dy/2, obj.Ly - obj.dy/2, obj.Ny)';              % Grid centers locations in y-direction
            obj.Zcm = linspace(obj.dz/2, obj.Lz - obj.dz/2, obj.Nz)';              % Grid centers locations in z-direction
            obj.Xim = linspace(0, obj.Lx, obj.Nx+1)';                              % Interface locations in x-direction
            obj.Yim = linspace(0, obj.Ly, obj.Ny+1)';                              % Interface locations in y-direction
            obj.Zim = linspace(0, obj.Lz, obj.Nz+1)';                              % Interface locations in z-direction
            obj.Epsilon = CompareAccuracy * min(obj.CellVolume)^(1/3);
            obj.AssignNeighbours();
            if strcmp(Type,'pEDFM')
                obj.pEDFM_alpha_Tx = zeros(obj.Nx+1, obj.Ny  , obj.Nz  );
                obj.pEDFM_alpha_Ty = zeros(obj.Nx  , obj.Ny+1, obj.Nz  );
                obj.pEDFM_alpha_Tz = zeros(obj.Nx  , obj.Ny  , obj.Nz+1);
            end
        end
        function AssignNeighbours(obj)
            % West Neighbors
            [i,j,k] = meshgrid( 0:obj.Nx-1 , 1:obj.Ny , 1:obj.Nz );  i(i<1) = nan;
            IW = (k-1).*obj.Nx.*obj.Ny + (j-1).*obj.Nx + i;
            IW = reshape(permute(IW,[2 1 3]),obj.N,1);
            
            % East Neighbors
            [i,j,k] = meshgrid( 2:obj.Nx+1 , 1:obj.Ny , 1:obj.Nz );  i(i>obj.Nx) = nan;
            IE = (k-1).*obj.Nx.*obj.Ny + (j-1).*obj.Nx + i;
            IE = reshape(permute(IE,[2 1 3]),obj.N,1);
            
            % South Neighbors
            [i,j,k] = meshgrid( 1:obj.Nx , 0:obj.Ny-1 , 1:obj.Nz );  j(j<1) = nan;
            IS = (k-1).*obj.Nx.*obj.Ny + (j-1).*obj.Nx + i;
            IS = reshape(permute(IS,[2 1 3]),obj.N,1);
            
            % North Neighbors
            [i,j,k] = meshgrid( 1:obj.Nx , 2:obj.Ny+1 , 1:obj.Nz );  j(j>obj.Ny) = nan;
            IN = (k-1).*obj.Nx.*obj.Ny + (j-1).*obj.Nx + i;
            IN = reshape(permute(IN,[2 1 3]),obj.N,1);
            
            % Bottom Neighbors
            [i,j,k] = meshgrid( 1:obj.Nx , 1:obj.Ny , 0:obj.Nz-1 );  k(k<1) = nan;
            IB = (k-1).*obj.Nx.*obj.Ny + (j-1).*obj.Nx + i;
            IB = reshape(permute(IB,[2 1 3]),obj.N,1);
            
            % Top Neighbors
            [i,j,k] = meshgrid( 1:obj.Nx , 1:obj.Ny , 2:obj.Nz+1 );  k(k>obj.Nz) = nan;
            IT = (k-1).*obj.Nx.*obj.Ny + (j-1).*obj.Nx + i;
            IT = reshape(permute(IT,[2 1 3]),obj.N,1);
            
            % Assembling the array for all the neighbors
            obj.Neighbours = [IW,IE,IS,IN,IB,IT];
            obj.Neighbours = num2cell(obj.Neighbours,2);
            obj.Neighbours = cellfun(@(x) x(~isnan(x)),obj.Neighbours ,'UniformOutput' ,false);
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
            disp( ['Grid  : ', num2str(obj.Nx), ' x ',  num2str(obj.Ny), ' x ', num2str(obj.Nz),...
                        ' = ', num2str(obj.Nx * obj.Ny * obj.Nz)] );
            disp( '---------------------------------------------------------' );
        end
    end
end