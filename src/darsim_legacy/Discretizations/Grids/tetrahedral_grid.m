%  Create the tetrahedral grid class
%--------------------------------------------------------------------------
% DARSim 2 Reservoir Simulator
% Author: Túlio de M. Cavalcante
% UFPE
%--------------------------------------------------------------------------
classdef tetrahedral_grid < grid_darsim
    properties
        Trans
        GTrans
        CDfactors
        ConnectivityMatrix
        Weights
        TetrahedralGridData
        pEDFM_alpha_Trans
        HeatTrans
        Volume
    end
    methods
        function obj = tetrahedral_grid(ReservoirProperties)
            obj.TetrahedralGridData = ReservoirProperties.TetrahedralGridData;
            obj.N = size(ReservoirProperties.TetrahedralGridData.Cell.Vertices,1);
            obj.Active = ones(obj.N, 1);
            obj.ActiveTime = ones(obj.N, 1);
            obj.Trans = zeros(size(ReservoirProperties.TetrahedralGridData.Face.Inner.Vertices,1),1);
            obj.pEDFM_alpha_Trans = zeros(size(ReservoirProperties.TetrahedralGridData.Face.Inner.Vertices,1),1);
            obj.HeatTrans = zeros(size(ReservoirProperties.TetrahedralGridData.Face.Inner.Vertices,1),1);
        end
        function Initialize(obj, Reservoir)
            obj.ComputeRockTransmissibilities(Reservoir);
            obj.ComputeWeights(Reservoir);
        end
        function ComputeRockTransmissibilities(obj, Reservoir)
            %--------------------------------------------------------------
            nface = size(obj.TetrahedralGridData.Face.Inner.Vertices,1);
            L = obj.TetrahedralGridData.Face.Inner.MontCell; R = obj.TetrahedralGridData.Face.Inner.JusCell;
            kL = obj.TetrahedralGridData.Cell.Centroid(L,:); kR = obj.TetrahedralGridData.Cell.Centroid(R,:);
            LR = kR - kL; distLR = sum(LR.*LR,2).^0.5;
            %--------------------------------------------------------------
%             permmatrixL = Reservoir.K(L,:); permmatrixR = Reservoir.K(R,:);
            permqueueL = zeros(3,3,nface); permqueueR = zeros(3,3,nface); 
%             permqueueL(1,1,:) = reshape(permmatrixL(:,1),[1 1 nface]); permqueueR(1,1,:) = reshape(permmatrixR(:,1),[1 1 nface]);
%             permqueueL(2,2,:) = reshape(permmatrixL(:,2),[1 1 nface]); permqueueR(2,2,:) = reshape(permmatrixR(:,2),[1 1 nface]);
%             permqueueL(3,3,:) = reshape(permmatrixL(:,3),[1 1 nface]); permqueueR(3,3,:) = reshape(permmatrixR(:,3),[1 1 nface]);
            % Tensor cheio ------------------------------------------------
            thetax = Reservoir.KRotAngle(1); thetay = Reservoir.KRotAngle(2); thetaz = Reservoir.KRotAngle(3);
            Rx = [1 0 0; 0 cos(thetax) -sin(thetax); 0 sin(thetax) cos(thetax)];
            Ry = [cos(thetay) 0 sin(thetay); 0 1 0; -sin(thetay) 0 cos(thetay)];
            Rz = [cos(thetaz) -sin(thetaz) 0; sin(thetaz) cos(thetaz) 0; 0 0 1];
            for icell = 1:size(L,1)
                permqueueL(:,:,icell) = Rz'*Ry'*Rx'*diag(Reservoir.K(L(icell),:))*Rx*Ry*Rz;
                permqueueR(:,:,icell) = Rz'*Ry'*Rx'*diag(Reservoir.K(R(icell),:))*Rx*Ry*Rz;
            end
            %--------------------------------------------------------------
            normals = obj.TetrahedralGridData.Face.Inner.Normals;
            KnL = tridot( obj, normals, permqueueL, normals );
            KnR = tridot( obj, normals, permqueueR, normals );
            %--------------------------------------------------------------
            veck = repmat([0 0 1],[size(normals,1) 1]);
            KgL = tridot( obj, veck, permqueueL, normals );
            KgR = tridot( obj, veck, permqueueR, normals );
            %--------------------------------------------------------------
            I = obj.TetrahedralGridData.Vertex.Coordinates(obj.TetrahedralGridData.Face.Inner.Vertices(:,1),:);
            J = obj.TetrahedralGridData.Vertex.Coordinates(obj.TetrahedralGridData.Face.Inner.Vertices(:,2),:);
            K = obj.TetrahedralGridData.Vertex.Coordinates(obj.TetrahedralGridData.Face.Inner.Vertices(:,3),:);
            JI = I - J; JK = K - J; tJI = cross( normals, JI ); tJK = cross( normals, JK );
            KtJIL = tridot( obj, normals, permqueueL, tJI ); KtJIR = tridot( obj, normals, permqueueR, tJI );
            KtJKL = tridot( obj, normals, permqueueL, tJK ); KtJKR = tridot( obj, normals, permqueueR, tJK );
            %--------------------------------------------------------------
            kL_I = I - kL; hL = abs(sum(kL_I.*normals,2))./sqrt(sum(normals.*normals,2)); 
            kR_I = I - kR; hR = abs(sum(kR_I.*normals,2))./sqrt(sum(normals.*normals,2));
            %--------------------------------------------------------------
            obj.Trans = (KnL.*KnR)./(KnL.*hR+KnR.*hL).*(sqrt(sum(normals.*normals,2)));
            obj.GTrans = abs((2./distLR).*(KgL.*KgR)./(KgL+KgR+1*(KgL+KgR==0)));
            obj.pEDFM_Trans_Correction();    
            DJI = (sum(tJI.*LR,2)./sum(normals.*normals,2))-(sqrt(sum(normals.*normals,2)).^(-1)).*(hL.*(KtJIL./KnL)+hR.*(KtJIR./KnR));
            DJK = (sum(tJK.*LR,2)./sum(normals.*normals,2))-(sqrt(sum(normals.*normals,2)).^(-1)).*(hL.*(KtJKL./KnL)+hR.*(KtJKR./KnR));
            % VertContrib = 0.5 *( -DJK * ( pJ - pI ) + DJI * ( pJ - pK ) )
            obj.CDfactors = [ -0.5*obj.Trans.*DJK -obj.Trans.*( 0.5*DJI - 0.5*DJK ) 0.5*obj.Trans.*DJI ];
            obj.ConnectivityMatrix = obj.TetrahedralGridData.Face.Inner.InnerFaceCellMatrix;
            obj.Volume = obj.TetrahedralGridData.Cell.Volume;
        end
        function ComputeWeights(obj, Reservoir)
            fprintf('---> Calculating weights for vertex unknowns in MPFA-D...')            
            % Global Least Square (Dong & Kang, 2021; 2022) ---------------
            obj.Weights = obj.TetrahedralGridData.Vertex.VertCellMatrix;
            ninface = sum(obj.TetrahedralGridData.Vertex.VertInnerFaceMatrix,2);
            nboface = sum(obj.TetrahedralGridData.Vertex.VertBoundFaceMatrix,2);
            % Tentativa de vetorizar --------------------------------------
            nelem = sum(obj.TetrahedralGridData.Vertex.VertCellMatrix,2);
            Q = obj.TetrahedralGridData.Vertex.Coordinates;
            k = obj.TetrahedralGridData.Cell.Centroid;
            [vertinface, infacevert] = find(obj.TetrahedralGridData.Vertex.VertInnerFaceMatrix==1);
            entinface = sortrows([vertinface infacevert],[1 2]);
            [vertboface, bofacevert] = find(obj.TetrahedralGridData.Vertex.VertBoundFaceMatrix==1);
            entboface = sortrows([vertboface bofacevert],[1 2]);
            [vertcell, cellvert] = find(obj.TetrahedralGridData.Vertex.VertCellMatrix==1); 
            entcell = sortrows([vertcell cellvert],[1 2]);
            %--------------------------------------------------------------
            thetax = Reservoir.KRotAngle(1); thetay = Reservoir.KRotAngle(2); thetaz = Reservoir.KRotAngle(3);
            Rx = [1 0 0; 0 cos(thetax) -sin(thetax); 0 sin(thetax) cos(thetax)];
            Ry = [cos(thetay) 0 sin(thetay); 0 1 0; -sin(thetay) 0 cos(thetay)];
            Rz = [cos(thetaz) -sin(thetaz) 0; sin(thetaz) cos(thetaz) 0; 0 0 1];
            %--------------------------------------------------------------
            cont = 0; inicel = 1; ininface = 1; inboface = 1;
            for i=1:size(Q,1)
                %----------------------------------------------------------
                finalel = inicel + nelem(i) - 1;
                velem = entcell(inicel:finalel,2)';
                M1 = zeros(size(velem,2),3*size(velem,2)+1);
                for j = 1:size(velem,2)
                    M1(j,3*(j-1)+1:3*j) = k(velem(j),:) - Q(i,:); M1(j,end) = 1; 
                end
                %----------------------------------------------------------
                M2 = zeros(ninface(i),3*size(velem,2)+1);
                M3 = zeros(ninface(i),3*size(velem,2)+1); 
                M4 = zeros(ninface(i)+nboface(i),3*size(velem,2)+1);
                N1 = eye(size(velem,2)); N2 = zeros(2*ninface(i),size(velem,2)); 
                N3 = zeros(ninface(i)+nboface(i),size(velem,2)); N = [N1; N2; N3];
                %----------------------------------------------------------
                finalif = ininface + ninface(i) - 1;
                facesconcorrnoi = entinface(ininface:finalif,2)';
                if isempty(facesconcorrnoi)==0
                    for j = 1:size(facesconcorrnoi,2)
                        fi = facesconcorrnoi(j); 
                        vertinface = obj.TetrahedralGridData.Face.Inner.Vertices(fi,:);
                        montelem = obj.TetrahedralGridData.Face.Inner.MontCell(fi); 
                        juselem = obj.TetrahedralGridData.Face.Inner.JusCell(fi);
                        nou = find(vertinface~=i);
                        t1 = Q(vertinface(nou(1)),:) - Q(i,:);
                        t2 = Q(vertinface(nou(2)),:) - Q(i,:);
                        nprov = cross(t1,t2); n = nprov/norm(nprov); 
                        nface = obj.TetrahedralGridData.Face.Inner.Normals(fi,:); nface = nface/norm(nface);
                        if round(norm(n - nface,8))~=0, testelem = juselem; juselem = montelem; montelem = testelem; end
                        cm = k(montelem,:); Km = Reservoir.K(montelem,:); cj = k(juselem,:); Kj = Reservoir.K(juselem,:); 
                        if obj.pEDFM_alpha_Trans(fi)==1, Km = 0*Km; Kj = 0*Kj; end
                        posm = find(velem==montelem); posj = find(velem==juselem);
                        M2(j,3*(posm-1)+1:3*posm) = t1; M2(j,3*(posj-1)+1:3*posj) = -t1; 
                        M3(j,3*(posm-1)+1:3*posm) = t2; M3(j,3*(posj-1)+1:3*posj) = -t2;
%                         M4(j,3*(posm-1)+1:3*posm) = n*diag(Km); M4(j,3*(posj-1)+1:3*posj) = -n*diag(Kj);
                        % Tensor cheio ------------------------------------
                        M4(j,3*(posm-1)+1:3*posm) = n*Rz'*Ry'*Rx'*diag(Km)*Rx*Ry*Rz; M4(j,3*(posj-1)+1:3*posj) = -n*Rz'*Ry'*Rx'*diag(Kj)*Rx*Ry*Rz;
                        %--------------------------------------------------
                    end
                end
                %----------------------------------------------------------
                finalbf = inboface + nboface(i) - 1;
                facesconcorrnoib = entboface(inboface:finalbf,2)';
                if isempty(facesconcorrnoib)==0
                    L3 = zeros(size(facesconcorrnoib,2),size(facesconcorrnoib,2));
                    F = zeros(size(facesconcorrnoib,2),1);
                    for j = 1:size(facesconcorrnoib,2)
                        fb = facesconcorrnoib(j); montelem = obj.TetrahedralGridData.Face.Bound.MontCell(fb);      
                        vertinface = obj.TetrahedralGridData.Face.Bound.Vertices(fb,:);
                        nou = find(vertinface~=i);
                        t1 = Q(vertinface(nou(1)),:) - Q(i,:);
                        t2 = Q(vertinface(nou(2)),:) - Q(i,:);
                        nprov = cross(t1,t2); n = nprov/norm(nprov);
                        nface = obj.TetrahedralGridData.Face.Bound.Normals(fb,:); nface = nface/norm(nface);
                        if round(norm(n - nface,8))~=0, n = nface; end
                        cm = k(montelem,:); Km = Reservoir.K(montelem,:); posm = find(velem==montelem);
%                         M4(j+ninface(i),3*(posm-1)+1:3*posm) = n*diag(Km);
                        % Tensor cheio ------------------------------------
                        M4(j+ninface(i),3*(posm-1)+1:3*posm) = n*Rz'*Ry'*Rx'*diag(Km)*Rx*Ry*Rz;
                        %--------------------------------------------------
                        vflux = 0; % DARSim unique boudary condition
                        F(j,1) = -vflux;
                        L3(j,j) = 1;
                    end
                    L1 = zeros(size(velem,2),size(L3,2)); L2 = zeros(3*ninface(i),size(L3,2)); %factor = 1;
                    factor = 10^(floor(log10(abs(max(max(max(M1)),max(max(max(M2)),max(max(M3))))/min(min(M4)))))); 
                    M4 = factor*M4; L3 = factor*L3; F = factor*F; L = [L1;L2;L3];
                else
                    factor = 10^(floor(log10(abs(max(max(max(M1)),max(max(max(M2)),max(max(M3))))/min(min(M4))))));
                    L = zeros(size(M1,1)+size(M2,1)+size(M3,1)+size(M4,1),1); %factor = 1;
                    M4 = factor*M4; F = 0; 
                end
                %----------------------------------------------------------
                M = [M1;M2;M3;M4]; e = zeros(1,3*size(velem,2)+1); e(end) = 1; 
                A = M'*M; B = M'*N; w = e*(A\B);          
                obj.Weights(i,velem) = w/sum(w);
                %----------------------------------------------------------
                per = floor(100*(i/size(Q,1)));
                if per>cont
                    if cont>0 && per<=10
                        fprintf('\b\b\b');
                    elseif per>10
                        fprintf('\b\b\b\b');
                    end
                    fprintf(' %u%%',floor(100*(i/size(Q,1))))
                    cont = cont + 1;                    
                end
                %----------------------------------------------------------
                inicel = inicel + nelem(i);
                ininface = ininface + ninface(i);
                inboface = inboface + nboface(i);
            end
            fprintf('\b\b\b\b\b');
            fprintf(' Done!\n')
            fprintf('---------------------------------------------------------\n')
        end
        function AddpEDFMCorrections(obj,pEDFM_alpha_Trans)
            obj.pEDFM_alpha_Trans = pEDFM_alpha_Trans;
        end
        function pEDFM_Trans_Correction(obj)
            obj.Trans = obj.Trans .* ( 1 - obj.pEDFM_alpha_Trans );
        end
        function CorrectTransmissibilitiesForpEDFM(obj)
            % Já foi feito;
        end
        function AddGridCoordinates(obj)
            obj.GridCoords = obj.TetrahedralGridData.Cell.Vertices;
        end
        function ComputeDepth(obj, alpha, Thickness, index)
            if (index) % Reservoir
                obj.Depth = max(obj.TetrahedralGridData.Cell.Centroid(:,3)) - obj.TetrahedralGridData.Cell.Centroid(:,3);
            end
        end
        function out = tridot( obj, vec1, tensor, vec2 )
            nface = size(vec1,1); N = zeros(3,3,nface); N(1,1,:) = reshape(vec2(:,1),[1 1 nface]);
            N(1,2,:) = reshape(vec2(:,2),[1 1 nface]); N(1,3,:) = reshape(vec2(:,3),[1 1 nface]);
            N(2,:,:) = N(1,:,:); N(3,:,:) = N(1,:,:);
            nT = zeros(3,1,nface); nT(1,1,:) = reshape(vec1(:,1),[1 1 nface]);
            nT(2,1,:) = reshape(vec1(:,2),[1 1 nface]); nT(3,1,:) = reshape(vec1(:,3),[1 1 nface]);
            out = reshape(sum(nT.*sum(tensor.*N,2),1),[nface 1 1])./sum(vec1.*vec1,2);
        end
    end           
end