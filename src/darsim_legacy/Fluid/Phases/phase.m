% Phase class
%--------------------------------------------------------------------------
% DARSim 2 Reservoir Simulator
% Author: Matteo Cusini, Yuhang Wang
% TU Delft
% Created: 14 July 2016
% Modified by Túlio de M. Cavalcante (UFPE) in 11 June 2022
%--------------------------------------------------------------------------
classdef phase < matlab.mixin.Heterogeneous & handle
    properties
        mu % Reference Viscosity
        sr % Irriducible saturation
    end
    methods
        function obj = phase()
        end
        function [Upwind, U] = UpWindAndRockFluxes(obj, Grid, P, RhoInt)
            switch class(Grid)
                case('corner_point_grid')
                    Neighbor1Index = Grid.CornerPointGridData.Internal_Faces.CellNeighbor1Index;
                    Neighbor2Index = Grid.CornerPointGridData.Internal_Faces.CellNeighbor2Index;
                    
                    % Gravitational velocities
                    Ug =  ( Grid.Depth(Neighbor2Index) - Grid.Depth(Neighbor1Index) ) .* RhoInt;
                    % Gravity test r2101
                    Ug =  -Ug;

                    % Compute 'rock' fluxes ([m^3/s])
                    U = - Grid.Trans .* ( +(P(Neighbor2Index) - P(Neighbor1Index)) + Ug );
                    
                    % Use velocity to build upwind operator
                    nc = Grid.N;
                    nf = length(Grid.Trans);
                    C = (U>=0) .* Grid.CornerPointGridData.Internal_Faces.CellNeighbor1Index + ...
                        (U< 0) .* Grid.CornerPointGridData.Internal_Faces.CellNeighbor2Index;
                    Upwind = sparse( (1:nf)', C , ones(nf,1) , nf , nc );
                    
                case('cartesian_grid')
                    Nx = Grid.Nx;
                    Ny = Grid.Ny;
                    Nz = Grid.Nz;
                    N = Grid.N;

                    % Gravitational velocities
                    depth = reshape(Grid.Depth, Nx, Ny, Nz);
                    Ugx = zeros(Nx+1, Ny, Nz);
                    Ugy = zeros(Nx, Ny+1, Nz);
                    Ugz = zeros(Nx, Ny, Nz+1);
                    Ugx(2:Nx,:,:) = (depth(1:Nx-1,:,:) - depth(2:Nx,:,:)) .* RhoInt.x(2:Nx,:,:);
                    Ugy(:,2:Ny,:) = (depth(:,1:Ny-1,:) - depth(:,2:Ny,:)) .* RhoInt.y(:,2:Ny,:);
                    Ugz(:,:,2:Nz) = (depth(:,:,1:Nz-1) - depth(:,:,2:Nz)) .* RhoInt.z(:,:,2:Nz);
                    
                    P = reshape(P, Nx, Ny, Nz);
                    %% Compute 'rock' fluxes ([m^3/s])
                    U.x = zeros(Nx+1,Ny,Nz);
                    U.y = zeros(Nx,Ny+1,Nz);
                    U.z = zeros(Nx,Ny,Nz+1);
                    U.x(2:Nx,:,:) = Grid.Tx(2:Nx,:,:) .* ( (P(1:Nx-1,:,:)-P(2:Nx,:,:)) - Ugx(2:Nx,:,:) );
                    U.y(:,2:Ny,:) = Grid.Ty(:,2:Ny,:) .* ( (P(:,1:Ny-1,:)-P(:,2:Ny,:)) - Ugy(:,2:Ny,:) );
                    U.z(:,:,2:Nz) = Grid.Tz(:,:,2:Nz) .* ( (P(:,:,1:Nz-1)-P(:,:,2:Nz)) - Ugz(:,:,2:Nz) );
                    
                    %convert_to_unstructured(U,Nx,Ny,Nz)

                    %% Use velocity to build upwind operator
                    L = reshape((U.x(2:Nx+1,:,:) >= 0), N, 1);
                    R = reshape((U.x(1:Nx,:,:) < 0), N, 1);
                    B = reshape((U.y(:,2:Ny+1,:) >= 0), N, 1);
                    T = reshape((U.y(:,1:Ny,:) < 0), N, 1);
                    Down = reshape((U.z(:,:,2:Nz+1) >= 0), N, 1);
                    Up = reshape((U.z(:,:,1:Nz) < 0), N, 1);
                    
                    DiagVecs = [L, R];
                    DiagIndx = [0, 1];
                    Upwind.x = spdiags(DiagVecs, DiagIndx, N, N);
                    DiagVecs = [B, T];
                    DiagIndx = [0, Nx];
                    Upwind.y = spdiags(DiagVecs, DiagIndx, N, N);
                    DiagVecs = [Down, Up];
                    DiagIndx = [0, Nx*Ny];
                    Upwind.z = spdiags(DiagVecs, DiagIndx, N, N);
                    
                case('tetrahedral_grid')
                    D = Grid.Depth;
                    
                    pR = P(Grid.TetrahedralGridData.Face.Inner.JusCell);
                    pL = P(Grid.TetrahedralGridData.Face.Inner.MontCell);
                    dR = D(Grid.TetrahedralGridData.Face.Inner.JusCell);
                    dL = D(Grid.TetrahedralGridData.Face.Inner.MontCell);
                    
                    I = Grid.TetrahedralGridData.Face.Inner.Vertices(:,1);
                    J = Grid.TetrahedralGridData.Face.Inner.Vertices(:,2);
                    K = Grid.TetrahedralGridData.Face.Inner.Vertices(:,3);
                    
                    p = Grid.Weights * P; pI = p(I); pJ = p(J); pK = p(K);
                    
                    % Gravitational velocities
                    Ug = ( Grid.GTrans .* dR - Grid.GTrans .* dL ) .* RhoInt;
                    
                    % Compute 'rock' fluxes ([m^3/s])
                    U = - Grid.Trans .* pR + Grid.Trans .* pL + Grid.CDfactors(:,1) .* pI + Grid.CDfactors(:,2) .* pJ + Grid.CDfactors(:,3) .* pK + Ug;
                                                            
                    % Use velocity to build upwind operator
                    nc = Grid.N;
                    nf = length(Grid.Trans);
                    C = (U>=0) .* Grid.TetrahedralGridData.Face.Inner.MontCell + ...
                        (U< 0) .* Grid.TetrahedralGridData.Face.Inner.JusCell;
                    Upwind = sparse( (1:nf)', C , ones(nf,1) , nf , nc );
                 case('unstructured_grid')
                    Neighbor1Index = Grid.mesh.face.inner.upstream;
                    Neighbor2Index = Grid.mesh.face.inner.downstream;
                    
                    % Gravitational velocities
                    Ug = ( Grid.Depth(Neighbor2Index) - Grid.Depth(Neighbor1Index) ) .* RhoInt;
                    
                    % Compute 'rock' fluxes ([m^3/s])
                    U = - Grid.T .* ( P(Neighbor2Index) - P(Neighbor1Index) - Ug );
                    
                    % Use velocity to build upwind operator
                    nc = Grid.N;
                    nf = length(Grid.T);
                    C = (U>=0) .* Grid.mesh.face.inner.upstream + ...
                        (U< 0) .* Grid.mesh.face.inner.downstream;
                    Upwind = sparse( (1:nf)', C , ones(nf,1) , nf , nc );
                 case('unstructured_fracture_grid')
                    Neighbor1Index = Grid.fracmesh.faces.inner.upstream;
                    Neighbor2Index = Grid.fracmesh.faces.inner.downstream;
                    
                    % Gravitational velocities
                    Ug = ( Grid.Depth(Neighbor2Index) - Grid.Depth(Neighbor1Index) ) .* RhoInt;
                    
                    % Compute 'rock' fluxes ([m^3/s])
                    U = - Grid.T .* ( P(Neighbor2Index) - P(Neighbor1Index) - Ug );
                    
                    % Use velocity to build upwind operator
                    nc = Grid.N;
                    nf = length(Grid.T);
                    C = (U>=0) .* Grid.fracmesh.faces.inner.upstream + ...
                        (U< 0) .* Grid.fracmesh.faces.inner.downstream;
                    Upwind = sparse( (1:nf)', C , ones(nf,1) , nf , nc );
            end
        end

        function U_D = UpWind_diffusion(obj, Grid, xph)
            switch class(Grid)
                case('unstructured_grid')

                    Neighbor1Index = Grid.mesh.face.inner.upstream;
                    Neighbor2Index = Grid.mesh.face.inner.downstream;

                    % diffusion driving force on each internal face
                    U_D = xph(Neighbor1Index) - xph(Neighbor2Index);
                
                case('unstructured_fracture_grid')

                    Neighbor1Index = Grid.fracmesh.faces.inner.upstream;
                    Neighbor2Index = Grid.fracmesh.faces.inner.downstream;

                    % diffusion driving force on each internal face
                    U_D = xph(Neighbor1Index) - xph(Neighbor2Index);

                case('corner_point_grid')

                    Neighbor1Index = Grid.CornerPointGridData.Internal_Faces.CellNeighbor1Index;
                    Neighbor2Index = Grid.CornerPointGridData.Internal_Faces.CellNeighbor2Index;

                    % diffusion driving force on each internal face
                    U_D = xph(Neighbor1Index) - xph(Neighbor2Index);
                    
                case('cartesian_grid')
                    Nx = Grid.Nx;
                    Ny = Grid.Ny;
                    Nz = Grid.Nz;
  
                    xph = reshape(xph, Nx, Ny, Nz);
                    %% Determine the direction of diffusion
                    U_D.x = zeros(Nx+1,Ny,Nz);
                    U_D.y = zeros(Nx,Ny+1,Nz);
                    U_D.z = zeros(Nx,Ny,Nz+1);
                    U_D.x(2:Nx,:,:) = (xph(1:Nx-1,:,:) - xph(2:Nx,:,:));
                    U_D.y(:,2:Ny,:) = (xph(:,1:Ny-1,:) - xph(:,2:Ny,:));
                    U_D.z(:,:,2:Nz) = (xph(:,:,1:Nz-1) - xph(:,:,2:Nz));

            end
        end

        function U_T = UpWind_T(obj, Grid, T)
            switch class(Grid)
                case('corner_point_grid')

                    Neighbor1Index = Grid.CornerPointGridData.Internal_Faces.CellNeighbor1Index;
                    Neighbor2Index = Grid.CornerPointGridData.Internal_Faces.CellNeighbor2Index;

                    % diffusion driving force on each internal face
                    U_T = T(Neighbor1Index) - T(Neighbor2Index);
                    
                case('cartesian_grid')
                    Nx = Grid.Nx;
                    Ny = Grid.Ny;
                    Nz = Grid.Nz;
  
                    T = reshape(T, Nx, Ny, Nz);
                    %% Determine the direction of diffusion
                    U_T.x = zeros(Nx+1,Ny,Nz);
                    U_T.y = zeros(Nx,Ny+1,Nz);
                    U_T.z = zeros(Nx,Ny,Nz+1);
                    U_T.x(2:Nx,:,:) = (T(1:Nx-1,:,:) - T(2:Nx,:,:));
                    U_T.y(:,2:Ny,:) = (T(:,1:Ny-1,:) - T(:,2:Ny,:));
                    U_T.z(:,:,2:Nz) = (T(:,:,1:Nz-1) - T(:,:,2:Nz));
            end
        end

    end
end