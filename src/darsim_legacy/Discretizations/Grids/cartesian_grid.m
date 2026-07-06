%  Cartesian grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%Last modified: 14 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef cartesian_grid < grid_darsim
    properties
        Nx
        Ny
        Nz
        dx
        dy
        dz
        Ax
        Ay
        Az
        Por
        Volume
        Tx
        Ty
        Tz
        Diffx
        Diffy
        Diffz
        DispT
        pEDFM_alpha_Tx = 0; % Correction for pEDFM Connectivities
        pEDFM_alpha_Ty = 0; % Correction for pEDFM Connectivities
        pEDFM_alpha_Tz = 0; % Correction for pEDFM Connectivities
        condTrans_x % transmisibility of heat conductivity
        condTrans_y % transmisibility of heat conductivity
        condTrans_z % transmisibility of heat conductivity
        I
        J
        K
        Centroids
    end
    methods
        function obj = cartesian_grid(n)
            obj.Nx = n(1);
            obj.Ny = n(2);
            obj.Nz = n(3);
            obj.N = prod(n);
            obj.Active = ones(obj.N, 1);
            obj.ActiveTime = ones(obj.N, 1);
            obj.Tx = zeros(obj.Nx+1, obj.Ny, obj.Nz);
            obj.Ty = zeros(obj.Nx, obj.Ny+1, obj.Nz);
            obj.Tz = zeros(obj.Nx, obj.Ny, obj.Nz+1);
            obj.pEDFM_alpha_Tx = zeros(obj.Nx+1, obj.Ny, obj.Nz);
            obj.pEDFM_alpha_Ty = zeros(obj.Nx, obj.Ny+1, obj.Nz);
            obj.pEDFM_alpha_Tz = zeros(obj.Nx, obj.Ny, obj.Nz+1);
            obj.condTrans_x = zeros(obj.Nx+1, obj.Ny, obj.Nz);
            obj.condTrans_y = zeros(obj.Nx, obj.Ny+1, obj.Nz);
            obj.condTrans_z = zeros(obj.Nx, obj.Ny, obj.Nz+1);
        end
        function Initialize(obj, Reservoir)
            obj.dx = Reservoir.Length/obj.Nx;
            obj.dy = Reservoir.Width/obj.Ny;
            obj.dz = Reservoir.Thickness/obj.Nz;
            obj.Ax = obj.dy * obj.dz;
            obj.Ay = obj.dx * obj.dz;
            obj.Az = obj.dx * obj.dy;
            obj.Volume = obj.dx * obj.dy * obj.dz * ones(obj.N,1);
            obj.AssignNeighbours();
            obj.ComputeRockTransmissibilities(Reservoir.K);
            obj.CoarseFactor = [1, 1, 1];
            obj.DualCoarseType = zeros(obj.N, 1);
            obj.Children = cell(obj.N, 1);
            obj.Por = Reservoir.Por;
            obj.AddCoordinates();
            obj.Depth = zeros(obj.N, 1);
        end
        function Initialize_fracture(obj, Reservoir)
            obj.dx = Reservoir.Length/obj.Nx;
            obj.dy = Reservoir.Width/obj.Ny;
            obj.dz = Reservoir.Thickness/obj.Nz;
            obj.Ax = obj.dy * obj.dz;
            obj.Ay = obj.dx * obj.dz;
            obj.Az = obj.dx * obj.dy;
            obj.Volume = obj.dx * obj.dy * obj.dz * ones(obj.N,1);
            obj.AssignNeighbours();
            obj.Por = Reservoir.Por;
            obj.ComputeRockTransmissibilities(Reservoir.K);
            obj.CoarseFactor = [1, 1, 1];
            obj.DualCoarseType = zeros(obj.N, 1);
            obj.Children = cell(obj.N, 1);
            obj.AddFractureCoordinates();
            obj.Depth = zeros(obj.N, 1);
        end
        function AddFractureCoordinates(obj)
            % Add I, J coordinates to Fracture Grid
            % obj.I = ones(obj.N, 1);
            % obj.J = ones(obj.N, 1);
            % obj.K = ones(obj.N, 1);

            X = reshape(obj.GridCoords(:,1), obj.Nx+1, obj.Ny+1);
            Y = reshape(obj.GridCoords(:,2), obj.Nx+1, obj.Ny+1);
            Z = reshape(obj.GridCoords(:,3), obj.Nx+1, obj.Ny+1);

            P11x = X(1:obj.Nx,   1:obj.Ny);   P11y = Y(1:obj.Nx,   1:obj.Ny);   P11z = Z(1:obj.Nx,   1:obj.Ny);
            P21x = X(2:obj.Nx+1, 1:obj.Ny);   P21y = Y(2:obj.Nx+1, 1:obj.Ny);   P21z = Z(2:obj.Nx+1, 1:obj.Ny);
            P22x = X(2:obj.Nx+1, 2:obj.Ny+1); P22y = Y(2:obj.Nx+1, 2:obj.Ny+1); P22z = Z(2:obj.Nx+1, 2:obj.Ny+1);
            P12x = X(1:obj.Nx,   2:obj.Ny+1); P12y = Y(1:obj.Nx,   2:obj.Ny+1); P12z = Z(1:obj.Nx,   2:obj.Ny+1);

            Cx = (P11x + P21x + P22x + P12x)/4;
            Cy = (P11y + P21y + P22y + P12y)/4;
            Cz = (P11z + P21z + P22z + P12z)/4;

            obj.Centroids = [Cx(:) Cy(:) Cz(:)];
        end
        function ComputeRockTransmissibilities(obj, K)
            %Harmonic average of permeability.
            Kx = reshape(K(:,1), obj.Nx, obj.Ny, obj.Nz);
            Ky = reshape(K(:,2), obj.Nx, obj.Ny, obj.Nz);
            Kz = reshape(K(:,3), obj.Nx, obj.Ny, obj.Nz);

            KHx = zeros(obj.Nx+1, obj.Ny, obj.Nz);
            KHy = zeros(obj.Nx, obj.Ny+1, obj.Nz);
            KHz = zeros(obj.Nx, obj.Ny, obj.Nz+1);
            KHx(2:obj.Nx,:,:) = 2*Kx(1:obj.Nx-1,:,:).*Kx(2:obj.Nx,:,:)./(Kx(1:obj.Nx-1,:,:)+Kx(2:obj.Nx,:,:));
            KHy(:,2:obj.Ny,:) = 2*Ky(:,1:obj.Ny-1,:).*Ky(:,2:obj.Ny,:)./(Ky(:,1:obj.Ny-1,:)+Ky(:,2:obj.Ny,:));
            KHz(:,:,2:obj.Nz) = 2*Kz(:,:,1:obj.Nz-1).*Kz(:,:,2:obj.Nz)./(Kz(:,:,1:obj.Nz-1)+Kz(:,:,2:obj.Nz));

            %Transmissibility
            obj.Tx(2:obj.Nx,:,:) = obj.Ax./obj.dx.*KHx(2:obj.Nx,:,:);
            obj.Ty(:,2:obj.Ny,:) = obj.Ay./obj.dy.*KHy(:,2:obj.Ny,:);
            obj.Tz(:,:,2:obj.Nz) = obj.Az./obj.dz.*KHz(:,:,2:obj.Nz);
        end
        function AddpEDFMCorrections(obj,pEDFM_alpha_Tx,pEDFM_alpha_Ty,pEDFM_alpha_Tz)
            obj.pEDFM_alpha_Tx = pEDFM_alpha_Tx;
            obj.pEDFM_alpha_Ty = pEDFM_alpha_Ty;
            obj.pEDFM_alpha_Tz = pEDFM_alpha_Tz;
        end
        function CorrectTransmissibilitiesForpEDFM(obj)
            obj.Tx(2:obj.Nx,:,:) = obj.Tx(2:obj.Nx,:,:) .* ( 1 - obj.pEDFM_alpha_Tx(2:obj.Nx,:,:) );
            obj.Ty(:,2:obj.Ny,:) = obj.Ty(:,2:obj.Ny,:) .* ( 1 - obj.pEDFM_alpha_Ty(:,2:obj.Ny,:) );
            obj.Tz(:,:,2:obj.Nz) = obj.Tz(:,:,2:obj.Nz) .* ( 1 - obj.pEDFM_alpha_Tz(:,:,2:obj.Nz) );
            %             obj.condTrans_x(2:obj.Nx,:,:) = obj.condTrans_x(2:obj.Nx,:,:) .* ( 1 - obj.pEDFM_alpha_Tx(2:obj.Nx,:,:) );
            %             obj.condTrans_y(:,2:obj.Ny,:) = obj.condTrans_y(:,2:obj.Ny,:) .* ( 1 - obj.pEDFM_alpha_Ty(:,2:obj.Ny,:) );
            %             obj.condTrans_z(:,:,2:obj.Nz) = obj.condTrans_z(:,:,2:obj.Nz) .* ( 1 - obj.pEDFM_alpha_Tz(:,:,2:obj.Nz) );
        end
        function ComputeHeatConductivitiyTransmissibilities(obj, CondEff)
            % Harmonic average of effective heat conductivity.
            % For now, we assume that heat conductivity is isotropic.
            % This means that CondEff is only one-column vector.
            Cx = reshape(CondEff(:,1), obj.Nx, obj.Ny, obj.Nz);
            Cy = reshape(CondEff(:,1), obj.Nx, obj.Ny, obj.Nz);
            Cz = reshape(CondEff(:,1), obj.Nx, obj.Ny, obj.Nz);

            CHx = zeros(obj.Nx+1, obj.Ny, obj.Nz);
            CHy = zeros(obj.Nx, obj.Ny+1, obj.Nz);
            CHz = zeros(obj.Nx, obj.Ny, obj.Nz+1);
            CHx(2:obj.Nx,:,:) = 2*Cx(1:obj.Nx-1,:,:).*Cx(2:obj.Nx,:,:)./(Cx(1:obj.Nx-1,:,:)+Cx(2:obj.Nx,:,:));
            CHy(:,2:obj.Ny,:) = 2*Cy(:,1:obj.Ny-1,:).*Cy(:,2:obj.Ny,:)./(Cy(:,1:obj.Ny-1,:)+Cy(:,2:obj.Ny,:));
            CHz(:,:,2:obj.Nz) = 2*Cz(:,:,1:obj.Nz-1).*Cz(:,:,2:obj.Nz)./(Cz(:,:,1:obj.Nz-1)+Cz(:,:,2:obj.Nz));

            %Transmissibility
            obj.condTrans_x(2:obj.Nx,:,:) = obj.Ax./obj.dx.*CHx(2:obj.Nx,:,:);
            obj.condTrans_y(:,2:obj.Ny,:) = obj.Ay./obj.dy.*CHy(:,2:obj.Ny,:);
            obj.condTrans_z(:,:,2:obj.Nz) = obj.Az./obj.dz.*CHz(:,:,2:obj.Nz);
        end
        function ComputeDiffusionTransmissibilities(obj, D)
            %Harmonic average
            nPh = size(D,2);
            obj.Diffx = zeros(obj.Nx+1, obj.Ny,   obj.Nz,   nPh);
            obj.Diffy = zeros(obj.Nx,   obj.Ny+1, obj.Nz,   nPh);
            obj.Diffz = zeros(obj.Nx,   obj.Ny,   obj.Nz+1, nPh);

            Dx = zeros(obj.Nx+1, obj.Ny, obj.Nz);
            Dy = zeros(obj.Nx, obj.Ny+1, obj.Nz);
            Dz = zeros(obj.Nx, obj.Ny, obj.Nz+1);

            for iph = 1:nPh

                Dph = reshape(D(:,iph), obj.Nx, obj.Ny, obj.Nz);

                % Harmonic average
                Dx(2:obj.Nx,:,:) = 2 * Dph(1:obj.Nx-1,:,:) .* Dph(2:obj.Nx,:,:)./ (Dph(1:obj.Nx-1,:,:) + Dph(2:obj.Nx,:,:) + 1e-10);
                Dy(:,2:obj.Ny,:) = 2 * Dph(:,1:obj.Ny-1,:) .* Dph(:,2:obj.Ny,:)./ (Dph(:,1:obj.Ny-1,:) + Dph(:,2:obj.Ny,:) + 1e-10);
                Dz(:,:,2:obj.Nz) = 2 * Dph(:,:,1:obj.Nz-1) .* Dph(:,:,2:obj.Nz)./ (Dph(:,:,1:obj.Nz-1) + Dph(:,:,2:obj.Nz) + 1e-10);

                % Transmissibility
                obj.Diffx(2:obj.Nx,:,:,iph) = obj.Ax ./ obj.dx .* Dx(2:obj.Nx,:,:);
                obj.Diffy(:,2:obj.Ny,:,iph) = obj.Ay ./ obj.dy .* Dy(:,2:obj.Ny,:);
                obj.Diffz(:,:,2:obj.Nz,iph) = obj.Az ./ obj.dz .* Dz(:,:,2:obj.Nz);
            end
        end
        function ComputeDispersionTransmissibilities(obj, E, S, UpWind, U, Mob, f)
            nph = size(U,1);
            
            E = reshape(E, obj.Nx, obj.Ny, obj.Nz );
            ExT = zeros(obj.Nx+1, obj.Ny, obj.Nz);
            EyT = zeros(obj.Nx, obj.Ny+1, obj.Nz);
            EzT = zeros(obj.Nx, obj.Ny, obj.Nz+1);

            ExT(2:obj.Nx,:,:) = (2 .* E(1:obj.Nx-1,:,:) .* E(2:obj.Nx,:,:) ./ (E(1:obj.Nx-1,:,:) + E(2:obj.Nx,:,:) + 1e-10)) .* (obj.Ax ./ obj.dx);
            EyT(:,2:obj.Ny,:) = (2 .* E(:,1:obj.Ny-1,:) .* E(:,2:obj.Ny,:) ./ (E(:,1:obj.Ny-1,:) + E(:,2:obj.Ny,:) + 1e-10)) .* (obj.Ay ./ obj.dy);
            EzT(:,:,2:obj.Nz) = (2 .* E(:,:,1:obj.Nz-1) .* E(:,:,2:obj.Nz) ./ (E(:,:,1:obj.Nz-1) + E(:,:,2:obj.Nz) + 1e-10)) .* (obj.Az ./ obj.dz);

            phase_tol = 1e-12;

            phase_presence_mult_x = zeros(obj.Nx+1, obj.Ny, obj.Nz);
            phase_presence_mult_y = zeros(obj.Nx, obj.Ny+1, obj.Nz);
            phase_presence_mult_z = zeros(obj.Nx, obj.Ny, obj.Nz+1);

            vx = zeros(obj.Nx+1, obj.Ny, obj.Nz);
            vy = zeros(obj.Nx, obj.Ny+1, obj.Nz);
            vz = zeros(obj.Nx, obj.Ny, obj.Nz+1);

            Umagx = zeros(obj.Nx+1, obj.Ny, obj.Nz);
            Umagy = zeros(obj.Nx, obj.Ny+1, obj.Nz);
            Umagz = zeros(obj.Nx, obj.Ny, obj.Nz+1);

            ux_on_x = zeros(obj.Nx+1, obj.Ny, obj.Nz);ux_on_y = zeros(obj.Nx, obj.Ny+1, obj.Nz);ux_on_z = zeros(obj.Nx, obj.Ny, obj.Nz+1);
            uy_on_x = zeros(obj.Nx+1, obj.Ny, obj.Nz);uy_on_y = zeros(obj.Nx, obj.Ny+1, obj.Nz);uy_on_z = zeros(obj.Nx, obj.Ny, obj.Nz+1);
            uz_on_x = zeros(obj.Nx+1, obj.Ny, obj.Nz);uz_on_y = zeros(obj.Nx, obj.Ny+1, obj.Nz);uz_on_z = zeros(obj.Nx, obj.Ny, obj.Nz+1);

            for iph = 1 : nph

                obj.DispT{iph,1+f}.x = zeros(obj.Nx+1, obj.Ny,   obj.Nz);
                obj.DispT{iph,1+f}.y = zeros(obj.Nx,   obj.Ny+1, obj.Nz);
                obj.DispT{iph,1+f}.z = zeros(obj.Nx,   obj.Ny,   obj.Nz+1);

                Sph = reshape(S(:, iph), obj.Nx, obj.Ny, obj.Nz);
                phase_presence_mult_x(2:obj.Nx,:,:) = double(Sph(1:obj.Nx-1,:,:) .* Sph(2:obj.Nx,:,:) > phase_tol );
                phase_presence_mult_y(:,2:obj.Ny,:) = double(Sph(:,1:obj.Ny-1,:) .* Sph(:,2:obj.Ny,:) > phase_tol );    
                phase_presence_mult_z(:,:,2:obj.Nz) = double(Sph(:,:,1:obj.Nz-1) .* Sph(:,:,2:obj.Nz) > phase_tol );

                Mx = reshape(UpWind{iph,1+f}.x*(Mob(:,iph)), obj.Nx, obj.Ny, obj.Nz);
                My = reshape(UpWind{iph,1+f}.y*(Mob(:,iph)), obj.Nx, obj.Ny, obj.Nz);
                Mz = reshape(UpWind{iph,1+f}.z*(Mob(:,iph)), obj.Nx, obj.Ny, obj.Nz);

                vx(2:obj.Nx,:,:) = U{iph,1+f}.x(2:obj.Nx,:,:) .* Mx(1:obj.Nx-1,:,:) / obj.Ax;
                vy(:,2:obj.Ny,:) = U{iph,1+f}.y(:,2:obj.Ny,:) .* My(:,1:obj.Ny-1,:) / obj.Ay;
                vz(:,:,2:obj.Nz) = U{iph,1+f}.z(:,:,2:obj.Nz) .* Mz(:,:,1:obj.Nz-1) / obj.Az;

                Uxc = 0.5 .* (vx(1:obj.Nx,:,:) + vx(2:obj.Nx+1,:,:));
                Uyc = 0.5 .* (vy(:,1:obj.Ny,:) + vy(:,2:obj.Ny+1,:));
                Uzc = 0.5 .* (vz(:,:,1:obj.Nz) + vz(:,:,2:obj.Nz+1));
 
                ux_on_x(2:obj.Nx,:,:) = vx(2:obj.Nx,:,:);
                uy_on_x(2:obj.Nx,:,:) = 0.5 .* (Uyc(1:obj.Nx-1,:,:) + Uyc(2:obj.Nx,:,:));
                uz_on_x(2:obj.Nx,:,:) = 0.5 .* (Uzc(1:obj.Nx-1,:,:) + Uzc(2:obj.Nx,:,:));

                % y-faces: between cell j-1 and j
                ux_on_y(:,2:obj.Ny,:) = 0.5 .* (Uxc(:,1:obj.Ny-1,:) + Uxc(:,2:obj.Ny,:));
                uy_on_y(:,2:obj.Ny,:) = vy(:,2:obj.Ny,:);
                uz_on_y(:,2:obj.Ny,:) = 0.5 .* (Uzc(:,1:obj.Ny-1,:) + Uzc(:,2:obj.Ny,:));

                % z-faces: between cell k-1 and k
                ux_on_z(:,:,2:obj.Nz) = 0.5 .* (Uxc(:,:,1:obj.Nz-1) + Uxc(:,:,2:obj.Nz));
                uy_on_z(:,:,2:obj.Nz) = 0.5 .* (Uyc(:,:,1:obj.Nz-1) + Uyc(:,:,2:obj.Nz));
                uz_on_z(:,:,2:obj.Nz) = vz(:,:,2:obj.Nz);

                Umagx(2:obj.Nx,:,:) = sqrt(ux_on_x(2:obj.Nx,:,:).^2 + uy_on_x(2:obj.Nx,:,:).^2 + uz_on_x(2:obj.Nx,:,:).^2 );
                Umagy(:,2:obj.Ny,:) = sqrt(ux_on_y(:,2:obj.Ny,:).^2 + uy_on_y(:,2:obj.Ny,:).^2 + uz_on_y(:,2:obj.Ny,:).^2 );
                Umagz(:,:,2:obj.Nz) = sqrt(ux_on_z(:,:,2:obj.Nz).^2 + uy_on_z(:,:,2:obj.Nz).^2 + uz_on_z(:,:,2:obj.Nz).^2 );

                obj.DispT{iph,1+f}.x(2:obj.Nx,:,:) =  phase_presence_mult_x(2:obj.Nx,:,:) .* ExT(2:obj.Nx,:,:) .* Umagx(2:obj.Nx,:,:);
                obj.DispT{iph,1+f}.y(:,2:obj.Ny,:) =  phase_presence_mult_y(:,2:obj.Ny,:) .* EyT(:,2:obj.Ny,:) .* Umagy(:,2:obj.Ny,:);
                obj.DispT{iph,1+f}.z(:,:,2:obj.Nz) =  phase_presence_mult_z(:,:,2:obj.Nz) .* EzT(:,:,2:obj.Nz) .* Umagz(:,:,2:obj.Nz);
            end
        end
        function AddCoordinates(obj)
            % Add I, J coordinates to Grid
            obj.I = ones(obj.N, 1);
            obj.J = ones(obj.N, 1);
            obj.K = ones(obj.N, 1);
            Jindexes = 1:1:obj.Ny;
            for k=1:obj.Nz
                for i=1:obj.Ny
                    a = obj.Nx*(i-1)+ (k-1)*obj.Nx*obj.Ny + 1;
                    obj.I(a:a+obj.Nx-1) = 1:1:obj.Nx;
                    obj.J(a:a+obj.Nx-1) = Jindexes(i)*ones(obj.Nx, 1);
                end
            end
            for i=1:obj.Nz
                obj.K((i-1)*obj.Nx*obj.Ny+1:i*obj.Nx*obj.Ny) = i*ones(obj.Nx*obj.Ny, 1);
            end
            XCM = linspace(obj.dx/2 , obj.Nx*obj.dx - obj.dx/2 , obj.Nx);
            YCM = linspace(obj.dy/2 , obj.Ny*obj.dy - obj.dy/2 , obj.Ny);
            ZCM = linspace(obj.dz/2 , obj.Nz*obj.dz - obj.dz/2 , obj.Nz);
            [x,y,z] = meshgrid(XCM,YCM,ZCM);
            x = permute(x,[2 1 3]);
            y = permute(y,[2 1 3]);

            obj.Centroids = [ x(:) , y(:) , z(:) ];
        end
        function AddGridCoordinates(obj)
            Xi = linspace(0,obj.dx*obj.Nx,obj.Nx+1)';
            Yi = linspace(0,obj.dy*obj.Ny,obj.Ny+1)';
            Zi = linspace(0,obj.dz*obj.Nz,obj.Nz+1)';
            [x,y,z] = meshgrid(Xi,Yi,Zi);
            x = permute(x,[2 1 3]);
            y = permute(y,[2 1 3]);
            z = permute(z,[2 1 3]);
            obj.GridCoords = [x(:),y(:),z(:)];
        end
        function ComputeDepth(obj, alpha, Thickness, index)
            if (index) % Reservoir
                % x_centres = (obj.I - 1/2) * obj.dx;
                % y_centres = (obj.J - 1/2) * obj.dy;
                z_centres = (obj.K - 1/2) * obj.dz;
                obj.Depth = Thickness - z_centres;
            else % Fracture

                obj.Depth = Thickness - obj.Centroids(:,3);
                % GridCoords_z = obj.GridCoords(:,3);
                % GridCoords_z = GridCoords_z(1:length(GridCoords_z)/2);
                % obj.Depth = Thickness - (GridCoords_z(1:end-1)+GridCoords_z(2:end))/2;
                % one column
                %                 GridCoords_z = obj.GridCoords(:,3);
                %                 GridCoords_z = GridCoords_z(1:length(GridCoords_z)/2);
                %                 obj.Depth = Thickness - (GridCoords_z(1:end-1)+GridCoords_z(2:end))/2;
                % multiple columns (square)
                %               obj.Depth = zeros(length(obj.I),1);
                %                 edge_length = sqrt(length(obj.GridCoords));
                %                 for i = 1:edge_length-1
                %                     start_index = 1+(i-1)*(edge_length-1);
                %                     end_index = i*(edge_length-1);
                %                     obj.Depth(start_index:end_index) = Thickness - (GridCoords_z(1+(i-1)*edge_length:i*edge_length-1)+GridCoords_z(2+(i-1)*edge_length:i*edge_length))/2;
                %                 end

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
    end
end