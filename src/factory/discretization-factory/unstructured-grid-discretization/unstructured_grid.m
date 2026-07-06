%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel;
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef unstructured_grid < grid_darsim
    properties
        mesh
        dx
        A
        T
        Por
        Diff
        pEDFM_alpha_T; % Correction for pEDFM Connectivities
        condTrans % transmisibility of heat conductivity
        I
        Centroids
        Volume
        ConnectivityMatrix
    end
    methods
        function obj = unstructured_grid(mesh)
            % since DARSim uses only null flux condition faces transm are
            % exclusively dependent of the inner faces
            obj.mesh = mesh;
            obj.GridCoords = mesh.vertex.coord;
            obj.Centroids = mesh.element.centroid;
            obj.Volume = mesh.element.volume;
            obj.N = mesh.element.num;
            obj.Active = ones(obj.N, 1);
            obj.ActiveTime = ones(obj.N, 1);
            obj.T = zeros(mesh.face.inner.num,1);
            obj.pEDFM_alpha_T = obj.T;
            obj.condTrans = obj.T;
            ind = [1:size(mesh.face.inner.upstream,1)]';
            auxmat = [ind mesh.face.inner.upstream  ones(size(mesh.face.inner.upstream)); ...
                ind mesh.face.inner.downstream  -1*ones(size(mesh.face.inner.upstream))];
            obj.ConnectivityMatrix = sparse(auxmat(:,1),auxmat(:,2),auxmat(:,3),max(auxmat(:,1)),max(auxmat(:,2)))';

        end
        function Initialize(obj, Reservoir)
            obj.AssignNeighbours();
            obj.ComputeRockTransmissibilities(Reservoir.K);
            obj.CoarseFactor = [1, 1, 1];
            obj.DualCoarseType = zeros(obj.N, 1);
            obj.Children = cell(obj.N, 1);
            obj.Depth = zeros(obj.N, 1);
            obj.Por = Reservoir.Por;
        end
        function AddGridCoordinates(obj)
            %uncessary function
        end
        function ComputeRockTransmissibilities(obj, K)
            % Artur's Harmonic average of permeability.
            % K_left = K(obj.Neighbours(:,1));
            % K_right = K(obj.Neighbours(:,2));
            % H_left = vecnorm(obj.mesh.face.inner.centroid - obj.mesh.element.centroid(obj.Neighbours(:,1),:),2,2);
            % H_right = vecnorm(obj.mesh.face.inner.centroid - obj.mesh.element.centroid(obj.Neighbours(:,2),:),2,2);
            % KH = K_left.*K_right./((K_left.*H_right)+(K_right.*H_left));
            % obj.T = obj.mesh.face.inner.area.*KH;

            %%% diagonal anisotropic and non-orthogonal TPFA
            %%% Keff=K*e*n  e: vector of face center-->cell center
            L = obj.Neighbours(:,1);
            R = obj.Neighbours(:,2);

            fc = obj.mesh.face.inner.centroid;    % face center [nF x dim]
            cL = obj.mesh.element.centroid(L,:);  % left cell center [nF x dim]
            cR = obj.mesh.element.centroid(R,:);  % right cell center [nF x dim]
            n   = obj.mesh.face.inner.normal;     % face normal [nF x dim]
            n = n ./ max(vecnorm(n,2,2),1e-12);   % normalize

            vL = fc - cL;
            vR = cR - fc;

            halfL = abs( sum(n .* (K(L,:) .* vL), 2) ./ sum(vL .^ 2, 2));
            halfR = abs( sum(n .* (K(R,:) .* vR), 2) ./ sum(vR .^ 2, 2));

            obj.T = obj.mesh.face.inner.area .* halfL .* halfR ./ (halfL + halfR);

        end
        function ComputeDiffusionTransmissibilities(obj, D)

            nPh = size(D, 2);

            L = obj.Neighbours(:,1);
            R = obj.Neighbours(:,2);

            fc = obj.mesh.face.inner.centroid;      % [nF x dim]
            cL = obj.mesh.element.centroid(L,:);    % [nF x dim]
            cR = obj.mesh.element.centroid(R,:);    % [nF x dim]

            n = obj.mesh.face.inner.normal;         % face normal
            n = n ./ max(vecnorm(n,2,2), 1e-12);    % unit normal

            vL = fc - cL;
            vR = cR - fc;

            area = obj.mesh.face.inner.area;
            nf = size(fc,1);

            obj.Diff = zeros(nf, nPh);

            geomL = abs(sum(n .* vL, 2) ./ (sum(vL.^2, 2) + 1e-30));   % ~ 1/dL
            geomR = abs(sum(n .* vR, 2) ./ (sum(vR.^2, 2) + 1e-30));   % ~ 1/dR

            for iph = 1:nPh
                halfL = D(L, iph) .* geomL;
                halfR = D(R, iph) .* geomR;

                obj.Diff(:, iph) = area .* halfL .* halfR ./ (halfL + halfR + 1e-30);
            end
        end
        function AddpEDFMCorrections(obj,pEDFM_alpha_T)
            obj.pEDFM_alpha_T(:) = pEDFM_alpha_T;
        end
        function CorrectTransmissibilitiesForpEDFM(obj)
            obj.T = obj.T .* (1 - obj.pEDFM_alpha_T);

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
            CHz(:,:,2:obj.Nz) = 2*Cz(:,:,1:obj.Nz-1).*Cz(:,:,2:obj.Nz)./(Cy(:,:,1:obj.Nz-1)+Cz(:,:,2:obj.Nz));

            %Transmissibility
            obj.condTrans_x(2:obj.Nx,:,:) = obj.Ax./obj.dx.*CHx(2:obj.Nx,:,:);
            obj.condTrans_y(:,2:obj.Ny,:) = obj.Ay./obj.dy.*CHy(:,2:obj.Ny,:);
            obj.condTrans_z(:,:,2:obj.Nz) = obj.Az./obj.dz.*CHz(:,:,2:obj.Nz);
        end
        function ComputeDepth(obj, alpha, Thickness, index)
            if (index) % Reservoir
                obj.Depth = Thickness - obj.Centroids(:,3);
            end
        end
        function AssignNeighbours(obj)
            obj.Neighbours = [obj.mesh.face.inner.upstream obj.mesh.face.inner.downstream];
        end
    end
end
