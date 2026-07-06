% Compositional Formulation base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Compositional_formulation < formulation
    properties
        NofComponents

        Utot
    end
    methods
        function obj = Compositional_formulation(n_components)
            obj@formulation();
            obj.NofComponents = n_components;
            %             obj.Tph = cell(n_components, 2);
            %             obj.Gph = cell(n_components, 2);

            obj.Tph = cell(4, 1); % c1ph1, c1ph2, c2ph1, c2ph2
            obj.Gph = cell(4, 1);
            obj.Dph = cell(4, 1); % Diffusion

            if isa(obj,'Overall_Composition_Energy_formulation')
                obj.Thph = cell(2,1);
                obj.Ghph = cell(2,1);
                obj.Tk   = cell(1,1);
            end

        end

        function q = ComputeSourceTerms(obj, N, Wells)
            q = zeros(N, obj.NofComponents);
            %Injectors
            for i=1:Wells.NofInj
                c = Wells.Inj{i}.Cells;
                q(c, :) = q(c, :) + Wells.Inj{i}.QComponents(:,:);
            end
            %Producers
            for i=1:Wells.NofProd
                c = Wells.Prod{i}.Cells;
                q(c, :) = q(c, :) + Wells.Prod{i}.QComponents(:,:);
            end
        end

        function q = ComputeSourceTerms_heat(obj, N, Wells)
            q = zeros(N, obj.NofPhases);
            %Injectors
            for i=1:Wells.NofInj
                c = Wells.Inj{i}.Cells;
                q(c, :) = q(c, :) + Wells.Inj{i}.Qh(:,:);
            end
            %Producers
            for i=1:Wells.NofProd
                c = Wells.Prod{i}.Cells;
                q(c, :) = q(c, :) + Wells.Prod{i}.Qh(:,:);
            end
        end

        % Formulation.AverageMassOnCoarseBlocks(ProductionSystem, FineGrid, FluidModel, ADMRest);
        %%
        function AverageMassOnCoarseBlocks_old(obj, Status, FluidModel, R, P)
            % Perform Average for ADM
            for i=1:obj.NofComponents
                z = Status.Properties(['z_',num2str(i)]);
                z_rest = R * z.Value;
                z.Value = P * z_rest;
            end

            % Update other unknwons as well
            obj.UpdatePhaseCompositions(Status, FluidModel);
        end
        %%
        function AverageMassOnCoarseBlocks(obj, DiscretizationModel, ProductionSystem, FineGrid, FluidModel, R)
            z = ProductionSystem.CreateGlobalVariables(FineGrid, obj.NofPhases, 'z_');

            % Perform Average for ADM
            delta = zeros(sum([FineGrid.N]), obj.NofComponents);
            for c = 1:obj.NofComponents - 1
                z_rest = R * z(:,c);
                z_av = R' * (z_rest ./ sum(R, 2));
                delta(:, c) = z_av - z(:,c);
            end
            delta(:, end) = -sum(delta(:,1:end-1), 2);

            for c=1:obj.NofComponents
                Start=1;
                End = FineGrid(1).N;
                z = ProductionSystem.Reservoir.State.Properties(['z_', num2str(c)]);
                z.update(delta(Start:End, c));
                for frac = 1:ProductionSystem.FracturesNetwork.NumOfFrac
                    Start = End + 1;
                    End = Start + FineGrid(frac+1).N - 1;
                    z = ProductionSystem.FracturesNetwork.Fractures(frac).State.Properties(['z_', num2str(c)]);
                    z.update(delta(Start:End, c));
                end
            end
            % Here FluidModel is useless but for compositional formulations it
            % is important
            % Update Remaining properties (x, S, rho, rhoT, Pc)
            obj.UpdatePhaseCompositions(ProductionSystem.Reservoir, DiscretizationModel, FluidModel, 1);
        end

        function UpdatePhaseCompositions(obj, Medium, DiscretizationModel, FluidModel, label)
            %% 2. Perform composition update ((x, ni) = f(T, p, z))
            % obj.SinglePhase = FluidModel.Flash(Reservoir.State);

            SinglePhase = FluidModel.Flash(Medium.State);
            status_SinglePhase = Medium.State.Properties('SinglePhase');
            status_SinglePhase.Value = SinglePhase;

            %% 3. Compute Densities (rho = rho(p, x, T))
            % FluidModel.ComputePhaseDensities(Reservoir.State, obj.SinglePhase);

            FluidModel.ComputePhaseDensities(Medium.State, Medium.State.Properties('SinglePhase').Value);

            %% 4. Compute Saturations (S = S(z, x))
            % FluidModel.ComputePhaseSaturation(Reservoir.State, obj.SinglePhase);

            FluidModel.ComputePhaseSaturation(Medium.State, Medium.State.Properties('SinglePhase').Value);

            %% 5. Compute Total Density (rhoT = rhoT(S, rho))
            FluidModel.ComputeTotalDensity(Medium.State);

            %% 6. Compute Pc (Pc = Pc(S))
            FluidModel.ComputePc(Medium, DiscretizationModel);

            if strcmpi(FluidModel.name,'Compositional_Energy')
            
                %% 7. Compute Hr (Hr = Cpr*T)
                FluidModel.ComputeRockEnthalpy(Medium);

                %% 8. Compute Deff
                FluidModel.ComputeThermalConductivity(Medium);

                %% 9. Compute Enthalpy (H = H(p, T))
                FluidModel.ComputePhaseEnthalpies(Medium.State);

            end

        end

        % corner_point_grid case is added
        function ComputeTotalFluxes(obj, ProductionSystem, DiscretizationModel)

            % Nx = DiscretizationModel.ReservoirGrid.Nx;
            % Ny = DiscretizationModel.ReservoirGrid.Ny;
            % Nz = DiscretizationModel.ReservoirGrid.Nz;
            N  = DiscretizationModel.ReservoirGrid.N;

            switch class(DiscretizationModel.ReservoirGrid)
                case('unstructured_grid')
                    N_Faces = length(obj.U{1,1});
                    obj.Utot = zeros(N_Faces,1);
                    for ph = 1:obj.NofPhases
                        %                         MobUpwind = obj.UpWind{ph,1}  * (obj.Mob(1:N, ph) .* obj.rho(:,ph) .* obj.x(:,ph));
                        MobUpwind = obj.UpWind{ph,1}  * (obj.Mob(1:N, ph));
                        obj.Utot = obj.Utot + obj.U{ph,1} .* MobUpwind;
                    end

                case('corner_point_grid')
                    N_Faces = length(obj.U{1,1});
                    obj.Utot = zeros(N_Faces,1);
                    for ph = 1:obj.NofPhases
                        %                         MobUpwind = obj.UpWind{ph,1}  * (obj.Mob(1:N, ph) .* obj.rho(:,ph) .* obj.x(:,ph));
                        MobUpwind = obj.UpWind{ph,1}  * (obj.Mob(1:N, ph));
                        obj.Utot = obj.Utot + obj.U{ph,1} .* MobUpwind;
                    end
                    %                     if ProductionSystem.Reservoir.State.Properties('V_tot').Plot
                    %                        % Not yet implemented
                    %                     end
                case('cartesian_grid')
                    Nx = DiscretizationModel.ReservoirGrid.Nx;
                    Ny = DiscretizationModel.ReservoirGrid.Ny;
                    Nz = DiscretizationModel.ReservoirGrid.Nz;
                    obj.Utot.x = zeros(Nx+1, Ny, Nz);
                    obj.Utot.y = zeros(Nx, Ny+1, Nz);
                    obj.Utot.z = zeros(Nx, Ny, Nz+1);

                    for ph = 1:obj.NofPhases
                        obj.Utot.x(2:Nx+1,:,:) = obj.Utot.x(2:Nx+1,:,:) + obj.U{ph, 1}.x(2:Nx+1,:,:) .* reshape(obj.UpWind{ph,1}.x *  (obj.Mob(1:N, ph)), Nx, Ny, Nz);  %- Ucap.x(2:Nx,:);
                        obj.Utot.y(:,2:Ny+1,:) = obj.Utot.y(:,2:Ny+1,:) + obj.U{ph, 1}.y(:,2:Ny+1,:) .* reshape(obj.UpWind{ph,1}.y *  (obj.Mob(1:N, ph)), Nx, Ny, Nz);  %- Ucap.y(:,2:Ny);
                        obj.Utot.z(:,:,2:Nz+1) = obj.Utot.z(:,:,2:Nz+1) + obj.U{ph, 1}.z(:,:,2:Nz+1) .* reshape(obj.UpWind{ph,1}.z *  (obj.Mob(1:N, ph)), Nx, Ny, Nz);  %- Ucap.y(:,2:Ny);
                    end

                otherwise
                    error("Discretization not available")
            end


            %             obj.Utot.x(2:Nx+1,:,:) = obj.Utot.x(2:Nx+1,:,:) + obj.U{2, 1}.x(2:Nx+1,:,:) .* reshape(obj.UpWind{2,1}.x *  (obj.Mob(1:N, 2)), Nx, Ny, Nz);  %- Ucap.x(2:Nx,:);
            %             % obj.Utot.y(:,2:Ny+1,:) = obj.Utot.y(:,2:Ny+1,:) + obj.U{ph, 1}.y(:,2:Ny+1,:) .* reshape(obj.UpWind{ph,1}.y *  (obj.Mob(1:N, ph)), Nx, Ny, Nz);  %- Ucap.y(:,2:Ny);
            %             obj.Utot.z(:,:,2:Nz+1) = obj.Utot.z(:,:,2:Nz+1) + obj.U{2, 1}.z(:,:,2:Nz+1) .* reshape(obj.UpWind{2,1}.z *  (obj.Mob(1:N, 2)), Nx, Ny, Nz);  %- Ucap.y(:,2:Ny);

            %             ux_temp = obj.Utot.x(2:Nx+1,:,:);
            %             status_ux = ProductionSystem.Reservoir.State.Properties('ux');
            %             status_ux.Value = ux_temp(:);
            %
            %             uz_temp = obj.Utot.z(:,:,2:Nz+1);
            %             status_uz = ProductionSystem.Reservoir.State.Properties('uz');
            %             status_uz.Value = uz_temp(:);
        end
        %%
        function conservative = CheckMassConservation(obj, DiscretizationModel)
            %Checks mass balance in all cells
            Nx = DiscretizationModel.ReservoirGrid.Nx;
            Ny = DiscretizationModel.ReservoirGrid.Ny;
            Nz = DiscretizationModel.ReservoirGrid.Nz;
            N  = DiscretizationModel.ReservoirGrid.N;

            conservative = 'Pass';
            maxUx = max(max(max(obj.Utot.x)));
            maxUy = max(max(max(obj.Utot.y)));
            maxUz = max(max(max(obj.Utot.z)));
            maxU = max([maxUx, maxUy, maxUz]);
            ux = reshape(obj.Utot.x(1:Nx,:,:) - obj.Utot.x(2:Nx+1,:,:), N, 1);
            uy = reshape(obj.Utot.y(:,1:Ny,:) - obj.Utot.y(:,2:Ny+1,:), N, 1);
            uz = reshape(obj.Utot.z(:,:,1:Nz) - obj.Utot.z(:,:,2:Nz+1), N, 1);

            Balance = ux + uy + uz; % wells are not included
            if norm(Balance, inf) > 1e-5
                conservative = 'Failed: Mass is NOT conserved!';
            end
            disp(['Mass conservation check:' num2str(conservative)]);
        end
        function CFL = ComputeCFLNumber(obj, ProductionSystem, DiscretizationModel, dt)
            %             N = DiscretizationModel.ReservoirGrid.N;
            %             pv = ProductionSystem.Reservoir.Por*DiscretizationModel.ReservoirGrid.Volume;
            %             P = zeros(N, obj.NofPhases);
            %             z = zeros(N, obj.NofComponents);
            %             rho = zeros(N, obj.NofPhases);
            %             x = zeros(N, obj.NofComponents*obj.NofPhases);
            %             % Copy values in local variables
            %             for i=1:obj.NofPhases
            %                 P(:, i) = ProductionSystem.Reservoir.State.Properties(['P_', num2str(i)]).Value;
            %                 rho(:, i) = ProductionSystem.Reservoir.State.Properties(['rho_', num2str(i)]).Value;
            %             end
            %             for i=1:obj.NofComponents
            %                 z(:, i) = ProductionSystem.Reservoir.State.Properties(['z_', num2str(i)]).Value;
            %                 for j=1:obj.NofComponents
            %                     x(:,(i-1)*obj.NofPhases + j) = ProductionSystem.Reservoir.State.Properties(['x_', num2str(i),'ph',num2str(j)]).Value;
            %                 end
            %             end
            %             rhoT =  ProductionSystem.Reservoir.State.Properties('rhoT').Value;
            %
            %              % Depths
            %             depth = DiscretizationModel.ReservoirGrid.Depth;
            %
            %             % Source terms
            %             q = obj.ComputeSourceTerms(N, ProductionSystem.Wells);
            %
            %             ThroughPut = zeros(N, obj.NofComponents);
            %             Mass = zeros(N, obj.NofComponents);
            %             for i=1:obj.NofComponents
            %                 obj.TransmissibilityMatrix(DiscretizationModel.ReservoirGrid, rho, obj.GravityModel.RhoInt, x(:,(i-1)*2+1:(i-1)*2+2), i);
            %                 % d1 = tril(obj.Tph{i, 1}, -1);
            %                 % d2 = tril(obj.Tph{i, 2}, -1);
            %                 % u1 = triu(obj.Tph{i, 1},  1);
            %                 % u2 = triu(obj.Tph{i, 2},  1);
            %                 d1 = tril(obj.Tph{2*i-1, 1}, -1);
            %                 d2 = tril(obj.Tph{2*i, 1}, -1);
            %                 u1 = triu(obj.Tph{2*i-1, 1},  1);
            %                 u2 = triu(obj.Tph{2*i, 1},  1);
            %
            %                 % Diag1 = diag(obj.Tph{i, 1});
            %                 % Diag2 = diag(obj.Tph{i, 2});
            %                 Diag1 = diag(obj.Tph{2*i-1, 1});
            %                 Diag2 = diag(obj.Tph{2*i, 1});
            %
            %                 % Assemble
            %                 D1 = diag(Diag1 + sum(u1, 2)) + d1;
            %                 D2 = diag(Diag2 + sum(u2, 2)) + d2;
            %                 U1 = diag(Diag1 + sum(d1, 2)) + u1;
            %                 U2 = diag(Diag2 + sum(d2, 2)) + u2;
            %                 ThroughPut(:,i) = ...
            %                            - min(D1 * P(:,1), 0)...    % Convective term
            %                            - min(D2 * P(:,2), 0)...
            %                            - min(U1 * P(:,1), 0)...    % Convective term
            %                            - min(U2 * P(:,2), 0)...
            %                            + max(q(:,i), 0);                       % Wells
            %                 Mass(:,i) = rhoT .* z(:,i);
            %             end
            %             Mass = max(Mass, 1e-10);
            %             ThroughPut(ThroughPut < 1e-10) = 0;
            %             Ratio = ThroughPut ./ Mass;
            %             CFL = dt * max(max(Ratio));
            CFL = 1;
        end

        %% corner_point_grid case is added to TransmissibilityMatrix function
        function TransmissibilityMatrix(obj, Grid, Reservoir, Rho, x, Index, i, f) % i represents component here; x are mole frictons of component i in different phases
            % obj.TransmissibilityMatrix(DiscretizationModel.ReservoirGrid, obj.rho(Index.Start:Index.End,:), obj.x(Index.Start:Index.End,(c-1)*2+1:(c-1)*2+2), Index, c, 0);
            switch class(Grid)
                case 'unstructured_grid'
                    %% Mass Flux Transmissibility and Gravity Matrix
                    N_Faces = length(Grid.T);
                    %% PHASE 1
                    % Apply Upwind Operator
                    Mob_rho_Upwind = obj.UpWind{1,1+f} * (obj.Mob(Index.Start:Index.End,1) .* Rho(:,1) .* x(:,1));
                    % Transmisibility Matrix
                    obj.Tph{2*i-1,1+f} = Grid.ConnectivityMatrix * sparse(1:N_Faces, 1:N_Faces, Grid.T .* Mob_rho_Upwind) * Grid.ConnectivityMatrix';
                    % Gravity Matrix
                    obj.Gph{2*i-1,1+f} = Grid.ConnectivityMatrix * sparse(1:N_Faces, 1:N_Faces, Grid.T .* Mob_rho_Upwind .* obj.GravityModel.RhoInt{1, 1+f}) * Grid.ConnectivityMatrix';
                    
                    % Diffusion Matrix
                    Diff1 = obj.S(Index.Start:Index.End, 1) .* Reservoir.Diff(:,1) .* Grid.Por;
                    Diff2 = obj.S(Index.Start:Index.End, 2) .* Reservoir.Diff(:,2) .* Grid.Por;

                    CellNeighbor1Index = Grid.mesh.face.inner.upstream;
                    CellNeighbor2Index = Grid.mesh.face.inner.downstream;

                    Grid.ComputeDiffusionTransmissibilities([Diff1, Diff2]);
                    obj.Dph{2*i-1,1+f} = Grid.ConnectivityMatrix * sparse(1:N_Faces, 1:N_Faces, (Rho(CellNeighbor1Index,1) + Rho(CellNeighbor2Index,1)) ./2 .* Grid.Diff(:,1)) * Grid.ConnectivityMatrix';

                    %% PHASE 2
                    % Apply Upwind Operator
                    Mob_rho_Upwind = obj.UpWind{2,1+f} * (obj.Mob(Index.Start:Index.End,2) .* Rho(:,2) .* x(:,2));
                    % Transmisibility Matrix
                    obj.Tph{2*i,1+f} = Grid.ConnectivityMatrix * sparse(1:N_Faces, 1:N_Faces, Grid.T.* Mob_rho_Upwind) * Grid.ConnectivityMatrix';
                    % Gravity Matrix
                    obj.Gph{2*i,1+f} = Grid.ConnectivityMatrix * sparse(1:N_Faces, 1:N_Faces, Grid.T.* Mob_rho_Upwind .* obj.GravityModel.RhoInt{2, 1+f}) * Grid.ConnectivityMatrix';
                    
                    % Diffusion Matrix
                    obj.Dph{2*i,1+f} = Grid.ConnectivityMatrix * sparse(1:N_Faces, 1:N_Faces, (Rho(CellNeighbor1Index,2) + Rho(CellNeighbor2Index,2)) ./2 .* Grid.Diff(:,2)) * Grid.ConnectivityMatrix';

                case('corner_point_grid')
                    N_Faces = length(Grid.Trans);

                    %% Mass Flux Transmissibility and Gravity Matrix
                    %% PHASE 1
                    % Apply Upwind Operator
                    Mob_rho_Upwind = obj.UpWind{1,1+f} * (obj.Mob(Index.Start:Index.End,1) .* Rho(:,1) .* x(:,1));

                    % Transmisibility Matrix
                    obj.Tph{2*i-1,1+f} = Grid.ConnectivityMatrix * sparse(1:N_Faces, 1:N_Faces, Grid.Trans .* Mob_rho_Upwind) * Grid.ConnectivityMatrix';

                    % Gravity Matrix
                    obj.Gph{2*i-1,1+f} = Grid.ConnectivityMatrix * sparse(1:N_Faces, 1:N_Faces, Grid.Trans .* Mob_rho_Upwind .* obj.GravityModel.RhoInt{1, 1+f}) * Grid.ConnectivityMatrix';

                    % Diffusion Matrix
                    Diff1 = obj.S(Index.Start:Index.End, 1) .* Reservoir.Diff(:,1) .* Grid.Por;
                    Diff2 = obj.S(Index.Start:Index.End, 2) .* Reservoir.Diff(:,2) .* Grid.Por;

                    CellNeighbor1Index = Grid.CornerPointGridData.Internal_Faces.CellNeighbor1Index;
                    CellNeighbor2Index = Grid.CornerPointGridData.Internal_Faces.CellNeighbor2Index;

                    Grid.ComputeDiffusionTransmissibilities([Diff1, Diff2]);
                    obj.Dph{2*i-1,1+f} = Grid.ConnectivityMatrix * sparse(1:N_Faces, 1:N_Faces, (Rho(CellNeighbor1Index,1) + Rho(CellNeighbor2Index,1)) ./2 .* Grid.Diff(:,1)) * Grid.ConnectivityMatrix';

                    %% PHASE 2
                    % Apply Upwind Operator
                    Mob_rho_Upwind = obj.UpWind{2,1+f} * (obj.Mob(Index.Start:Index.End,2) .* Rho(:,2) .* x(:,2));

                    % Transmisibility Matrix
                    obj.Tph{2*i,1+f} = Grid.ConnectivityMatrix * sparse(1:N_Faces, 1:N_Faces, Grid.Trans .* Mob_rho_Upwind) * Grid.ConnectivityMatrix';
                    % obj.Tph{2*i,1+f} = Grid.ConnectivityMatrix * spdiags(Grid.Trans .* Mob_rho_Upwind, 0, N_Faces, N_Faces) * Grid.ConnectivityMatrix';

                    % Gravity Matrix
                    obj.Gph{2*i,1+f} = Grid.ConnectivityMatrix * sparse(1:N_Faces, 1:N_Faces, Grid.Trans .* Mob_rho_Upwind .* obj.GravityModel.RhoInt{2, 1+f}) * Grid.ConnectivityMatrix';
                    
                    % Diffusion Matrix
                    obj.Dph{2*i,1+f} = Grid.ConnectivityMatrix * sparse(1:N_Faces, 1:N_Faces, (Rho(CellNeighbor1Index,2) + Rho(CellNeighbor2Index,2)) ./2 .* Grid.Diff(:,2)) * Grid.ConnectivityMatrix';

                case('cartesian_grid')
                    %%%Transmissibility matrix construction
                    Tx = zeros(Grid.Nx+1, Grid.Ny, Grid.Nz);
                    Ty = zeros(Grid.Nx, Grid.Ny+1, Grid.Nz);
                    Tz = zeros(Grid.Nx, Grid.Ny, Grid.Nz+1);

                    % add diffusion flux 
                    DTx = zeros(Grid.Nx+1, Grid.Ny, Grid.Nz);
                    DTy = zeros(Grid.Nx, Grid.Ny+1, Grid.Nz);
                    DTz = zeros(Grid.Nx, Grid.Ny, Grid.Nz+1);

                    %% PHASE 1
                    %Apply upwind operator
                    Mupx = obj.UpWind{1,1+f}.x*(obj.Mob(Index.Start:Index.End,1) .* Rho(:,1) .* x(:,1)); % upstream weighted phase-based transimissibility
                    Mupy = obj.UpWind{1,1+f}.y*(obj.Mob(Index.Start:Index.End,1) .* Rho(:,1) .* x(:,1)); % (upwind is a N x N matrix, filled with  a pair of 1 and 0 at each row, phase.m file)
                    Mupz = obj.UpWind{1,1+f}.z*(obj.Mob(Index.Start:Index.End,1) .* Rho(:,1) .* x(:,1)); % N x 1: ρ_α^up * (k_(r,α)^up)/(μ_α^up) * x_cα
                    Mupx = reshape(Mupx, Grid.Nx, Grid.Ny, Grid.Nz);
                    Mupy = reshape(Mupy, Grid.Nx, Grid.Ny, Grid.Nz);
                    Mupz = reshape(Mupz, Grid.Nx, Grid.Ny, Grid.Nz);

                    Tx(2:Grid.Nx,:,:)= Grid.Tx(2:Grid.Nx,:,:).*Mupx(1:Grid.Nx-1,:,:); % multiplying geometric part of transmissibility over inner interfaces. Nx in Mupx is for right boundary.
                    Ty(:,2:Grid.Ny,:)= Grid.Ty(:,2:Grid.Ny,:).*Mupy(:,1:Grid.Ny-1,:);
                    Tz(:,:,2:Grid.Nz)= Grid.Tz(:,:,2:Grid.Nz).*Mupz(:,:,1:Grid.Nz-1);

                    %Construct matrix
                    x1 = reshape(Tx(1:Grid.Nx,:,:), Grid.N, 1);   % West boundary,  transmissbiility between I and IW   (x: W->E, y: S->N, z: Bot->Top)
                    x2 = reshape(Tx(2:Grid.Nx+1,:,:), Grid.N, 1); % East boundary,  transmissbiility between I and IE
                    y1 = reshape(Ty(:,1:Grid.Ny,:), Grid.N, 1);   % South boundary, transmissbiility between I and IS
                    y2 = reshape(Ty(:,2:Grid.Ny+1,:), Grid.N, 1); % North boundary, transmissbiility between I and IN
                    z1 = reshape(Tz(:,:,1:Grid.Nz), Grid.N, 1);   % Bottom boundary,transmissbiility between I and IBot
                    z2 = reshape(Tz(:,:,2:Grid.Nz+1), Grid.N, 1); % Top boundary,   transmissbiility between I and ITop

                    DiagVecs = [-z2,-y2,-x2,z2+y2+x2+y1+x1+z1,-x1,-y1,-z1];
                    DiagIndx = [-Grid.Nx*Grid.Ny, -Grid.Nx, -1, 0, 1, Grid.Nx, Grid.Nx*Grid.Ny];
                    % obj.Tph{i,1} = spdiags(DiagVecs, DiagIndx, Grid.N, Grid.N);
                    obj.Tph{2*i-1,1+f} = spdiags(DiagVecs, DiagIndx, Grid.N, Grid.N); % i is component

                    % Gravity Matrix
                    Tx(2:Grid.Nx,:,:)= Tx(2:Grid.Nx,:,:) .* obj.GravityModel.RhoInt{1,1+f}.x(2:Grid.Nx,:,:); % inner interfaces: obj.GravityModel.RhoInt: rho_average * g
                    Ty(:,2:Grid.Ny,:)= Ty(:,2:Grid.Ny,:) .* obj.GravityModel.RhoInt{1,1+f}.y(:,2:Grid.Ny,:);
                    Tz(:,:,2:Grid.Nz)= Tz(:,:,2:Grid.Nz) .* obj.GravityModel.RhoInt{1,1+f}.z(:,:,2:Grid.Nz);

                    %Construct matrix
                    x1 = reshape(Tx(1:Grid.Nx,:,:), Grid.N, 1);
                    x2 = reshape(Tx(2:Grid.Nx+1,:,:), Grid.N, 1);
                    y1 = reshape(Ty(:,1:Grid.Ny,:), Grid.N, 1);
                    y2 = reshape(Ty(:,2:Grid.Ny+1,:), Grid.N, 1);
                    z1 = reshape(Tz(:,:,1:Grid.Nz), Grid.N, 1);
                    z2 = reshape(Tz(:,:,2:Grid.Nz+1), Grid.N, 1);

                    DiagVecs = [-z2,-y2,-x2,z2+y2+x2+y1+x1+z1,-x1,-y1,-z1];
                    DiagIndx = [-Grid.Nx*Grid.Ny, -Grid.Nx, -1, 0, 1, Grid.Nx, Grid.Nx*Grid.Ny];
                    % obj.Gph{i,1} = spdiags(DiagVecs, DiagIndx, Grid.N, Grid.N);
                    obj.Gph{2*i-1,1+f} = spdiags(DiagVecs, DiagIndx, Grid.N, Grid.N);
                    
                    % Diffusion matrix
                    rho_1 = reshape(Rho(:,1), Grid.Nx, Grid.Ny, Grid.Nz);
                    Diff1 = obj.S(Index.Start:Index.End, 1) .* Reservoir.Diff(:,1) .* Grid.Por;
                    Diff2 = obj.S(Index.Start:Index.End, 2) .* Reservoir.Diff(:,2) .* Grid.Por;
                    Grid.ComputeDiffusionTransmissibilities([Diff1, Diff2]);
                    Grid.ComputeDispersionTransmissibilities(Reservoir.Disp, obj.S(Index.Start:Index.End, :), obj.UpWind, obj.U, obj.Mob(Index.Start:Index.End, :), f);
                   
                    DTx(2:Grid.Nx,:,:) = (rho_1(1:Grid.Nx-1,:,:) + rho_1(2:Grid.Nx,:,:)) / 2 .* (Grid.Diffx(2:Grid.Nx,:,:,1) + Grid.DispT{1,1+f}.x(2:Grid.Nx,:,:));
                    DTy(:,2:Grid.Ny,:) = (rho_1(:,1:Grid.Ny-1,:) + rho_1(:,2:Grid.Ny,:)) / 2 .* (Grid.Diffy(:,2:Grid.Ny,:,1) + Grid.DispT{1,1+f}.y(:,2:Grid.Ny,:));
                    DTz(:,:,2:Grid.Nz) = (rho_1(:,:,1:Grid.Nz-1) + rho_1(:,:,2:Grid.Nz)) / 2 .* (Grid.Diffz(:,:,2:Grid.Nz,1) + Grid.DispT{1,1+f}.z(:,:,2:Grid.Nz));

                    x1 = reshape(DTx(1:Grid.Nx,:,:), Grid.N, 1);  
                    x2 = reshape(DTx(2:Grid.Nx+1,:,:), Grid.N, 1);
                    y1 = reshape(DTy(:,1:Grid.Ny,:), Grid.N, 1); 
                    y2 = reshape(DTy(:,2:Grid.Ny+1,:), Grid.N, 1);
                    z1 = reshape(DTz(:,:,1:Grid.Nz), Grid.N, 1); 
                    z2 = reshape(DTz(:,:,2:Grid.Nz+1), Grid.N, 1);

                    DiagVecs = [-z2,-y2,-x2,z2+y2+x2+y1+x1+z1,-x1,-y1,-z1];
                    DiagIndx = [-Grid.Nx*Grid.Ny, -Grid.Nx, -1, 0, 1, Grid.Nx, Grid.Nx*Grid.Ny];
                    obj.Dph{2*i-1,1+f} = spdiags(DiagVecs, DiagIndx, Grid.N, Grid.N); % i is component

                    %% PHASE 2
                    %Apply upwind operator
                    Mupx = obj.UpWind{2,1+f}.x*(obj.Mob(Index.Start:Index.End,2) .* Rho(:,2) .* x(:,2));
                    Mupy = obj.UpWind{2,1+f}.y*(obj.Mob(Index.Start:Index.End,2) .* Rho(:,2) .* x(:,2));
                    Mupz = obj.UpWind{2,1+f}.z*(obj.Mob(Index.Start:Index.End,2) .* Rho(:,2) .* x(:,2));
                    Mupx = reshape(Mupx, Grid.Nx, Grid.Ny, Grid.Nz);
                    Mupy = reshape(Mupy, Grid.Nx, Grid.Ny, Grid.Nz);
                    Mupz = reshape(Mupz, Grid.Nx, Grid.Ny, Grid.Nz);

                    Tx(2:Grid.Nx,:,:) = Grid.Tx(2:Grid.Nx,:,:).*Mupx(1:Grid.Nx-1,:,:);
                    Ty(:,2:Grid.Ny,:) = Grid.Ty(:,2:Grid.Ny,:).*Mupy(:,1:Grid.Ny-1,:);
                    Tz(:,:,2:Grid.Nz) = Grid.Tz(:,:,2:Grid.Nz).*Mupz(:,:,1:Grid.Nz-1);

                    %Construct matrix
                    x1 = reshape(Tx(1:Grid.Nx,:,:), Grid.N, 1);
                    x2 = reshape(Tx(2:Grid.Nx+1,:,:), Grid.N, 1);
                    y1 = reshape(Ty(:,1:Grid.Ny,:), Grid.N, 1);
                    y2 = reshape(Ty(:,2:Grid.Ny+1,:), Grid.N, 1);
                    z1 = reshape(Tz(:,:,1:Grid.Nz), Grid.N, 1);
                    z2 = reshape(Tz(:,:,2:Grid.Nz+1), Grid.N, 1);

                    DiagVecs = [-z2,-y2,-x2,z2+y2+x2+y1+x1+z1,-x1,-y1,-z1];
                    DiagIndx = [-Grid.Nx*Grid.Ny, -Grid.Nx, -1, 0, 1, Grid.Nx, Grid.Nx*Grid.Ny];
                    % obj.Tph{i,2} = spdiags(DiagVecs, DiagIndx, Grid.N, Grid.N);
                    obj.Tph{2*i,1+f} = spdiags(DiagVecs, DiagIndx, Grid.N, Grid.N);

                    % Gravity Matrix
                    Tx(2:Grid.Nx,:,:) = Tx(2:Grid.Nx,:,:) .* obj.GravityModel.RhoInt{2,1+f}.x(2:Grid.Nx,:,:);
                    Ty(:,2:Grid.Ny,:) = Ty(:,2:Grid.Ny,:) .* obj.GravityModel.RhoInt{2,1+f}.y(:,2:Grid.Ny,:);
                    Tz(:,:,2:Grid.Nz) = Tz(:,:,2:Grid.Nz) .* obj.GravityModel.RhoInt{2,1+f}.z(:,:,2:Grid.Nz);

                    %Construct matrix
                    x1 = reshape(Tx(1:Grid.Nx,:,:), Grid.N, 1);
                    x2 = reshape(Tx(2:Grid.Nx+1,:,:), Grid.N, 1);
                    y1 = reshape(Ty(:,1:Grid.Ny,:), Grid.N, 1);
                    y2 = reshape(Ty(:,2:Grid.Ny+1,:), Grid.N, 1);
                    z1 = reshape(Tz(:,:,1:Grid.Nz), Grid.N, 1);
                    z2 = reshape(Tz(:,:,2:Grid.Nz+1), Grid.N, 1);

                    DiagVecs = [-z2,-y2,-x2,z2+y2+x2+y1+x1+z1,-x1,-y1,-z1];
                    DiagIndx = [-Grid.Nx*Grid.Ny, -Grid.Nx, -1, 0, 1, Grid.Nx, Grid.Nx*Grid.Ny];
                    % obj.Gph{i,2} = spdiags(DiagVecs, DiagIndx, Grid.N, Grid.N);
                    obj.Gph{2*i,1+f} = spdiags(DiagVecs, DiagIndx, Grid.N, Grid.N);

                    % Diffusion matrix
                    rho_2 = reshape(Rho(:,2), Grid.Nx, Grid.Ny, Grid.Nz);
                   
                    DTx(2:Grid.Nx,:,:) = (rho_2(1:Grid.Nx-1,:,:) + rho_2(2:Grid.Nx,:,:)) / 2 .* (Grid.Diffx(2:Grid.Nx,:,:,2) + Grid.DispT{2,1+f}.x(2:Grid.Nx,:,:));
                    DTy(:,2:Grid.Ny,:) = (rho_2(:,1:Grid.Ny-1,:) + rho_2(:,2:Grid.Ny,:)) / 2 .* (Grid.Diffy(:,2:Grid.Ny,:,2) + Grid.DispT{2,1+f}.y(:,2:Grid.Ny,:));
                    DTz(:,:,2:Grid.Nz) = (rho_2(:,:,1:Grid.Nz-1) + rho_2(:,:,2:Grid.Nz)) / 2 .* (Grid.Diffz(:,:,2:Grid.Nz,2) + Grid.DispT{2,1+f}.z(:,:,2:Grid.Nz));

                    x1 = reshape(DTx(1:Grid.Nx,:,:), Grid.N, 1);  
                    x2 = reshape(DTx(2:Grid.Nx+1,:,:), Grid.N, 1);
                    y1 = reshape(DTy(:,1:Grid.Ny,:), Grid.N, 1); 
                    y2 = reshape(DTy(:,2:Grid.Ny+1,:), Grid.N, 1);
                    z1 = reshape(DTz(:,:,1:Grid.Nz), Grid.N, 1); 
                    z2 = reshape(DTz(:,:,2:Grid.Nz+1), Grid.N, 1);

                    DiagVecs = [-z2,-y2,-x2,z2+y2+x2+y1+x1+z1,-x1,-y1,-z1];
                    DiagIndx = [-Grid.Nx*Grid.Ny, -Grid.Nx, -1, 0, 1, Grid.Nx, Grid.Nx*Grid.Ny];
                    obj.Dph{2*i,1+f} = spdiags(DiagVecs, DiagIndx, Grid.N, Grid.N); % i is component

                case('unstructured_fracture_grid')
                    %% Mass Flux Transmissibility and Gravity Matrix
                    N_Faces = length(Grid.T);
                    CellNeighbor1Index = Grid.fracmesh.faces.inner.upstream;
                    CellNeighbor2Index = Grid.fracmesh.faces.inner.downstream;
                    %% PHASE 1
                    % Apply Upwind Operator
                    Mob_rho_Upwind = obj.UpWind{1,1+f} * (obj.Mob(Index.Start:Index.End,1) .* Rho(:,1) .* x(:,1));
                    % Transmisibility Matrix
                    obj.Tph{2*i-1,1+f} = Grid.ConnectivityMatrix * sparse(1:N_Faces, 1:N_Faces, Grid.T .* Mob_rho_Upwind) * Grid.ConnectivityMatrix';
                    % Gravity Matrix
                    obj.Gph{2*i-1,1+f} = Grid.ConnectivityMatrix * sparse(1:N_Faces, 1:N_Faces, Grid.T .* Mob_rho_Upwind .* obj.GravityModel.RhoInt{1, 1+f}) * Grid.ConnectivityMatrix';
                    %% PHASE 2
                    % Apply Upwind Operator
                    Mob_rho_Upwind = obj.UpWind{2,1+f} * (obj.Mob(Index.Start:Index.End,2) .* Rho(:,2) .* x(:,2));
                    % Transmisibility Matrix
                    obj.Tph{2*i,1+f} = Grid.ConnectivityMatrix * sparse(1:N_Faces, 1:N_Faces, Grid.T.* Mob_rho_Upwind) * Grid.ConnectivityMatrix';
                    % Gravity Matrix
                    obj.Gph{2*i,1+f} = Grid.ConnectivityMatrix * sparse(1:N_Faces, 1:N_Faces, Grid.T.* Mob_rho_Upwind .* obj.GravityModel.RhoInt{2, 1+f}) * Grid.ConnectivityMatrix';
                    
                    Diff1 = obj.S(Index.Start:Index.End, 1) .* Reservoir.Diff(:,1) .* Grid.Por;
                    Diff2 = obj.S(Index.Start:Index.End, 2) .* Reservoir.Diff(:,2) .* Grid.Por;

                    Grid.ComputeDiffusionTransmissibilities([Diff1, Diff2]);

                    obj.Dph{2*i-1,1+f} = Grid.ConnectivityMatrix * sparse(1:N_Faces, 1:N_Faces, (Rho(CellNeighbor1Index,1) + Rho(CellNeighbor2Index,1)) ./2 .* Grid.Diff(:,1)) * Grid.ConnectivityMatrix';
                    obj.Dph{2*i,1+f} = Grid.ConnectivityMatrix * sparse(1:N_Faces, 1:N_Faces, (Rho(CellNeighbor1Index,2) + Rho(CellNeighbor2Index,2)) ./2 .* Grid.Diff(:,2)) * Grid.ConnectivityMatrix';

                otherwise
                    error("Discretization not defined.")
            end
        end

        function ConvectiveHeatTransmissibilityMatrix(obj, Grid, Index, rho, h, ph, f)
            % Th : Heat convection transmissibility (Tph*rho)
            % Gh : Heat convection gravity
            switch class(Grid)
                case('corner_point_grid')
                    N_Faces = length(Grid.Trans);

                    %% Heat Convection Flux Transmissibility Matrix
                    % Apply Upwind Operator
                    Mob_rho_h_Upwind = obj.UpWind{ph,1} * (obj.Mob(Index.Start:Index.End, ph) .* rho(:, ph) .* h(Index.Start:Index.End, ph));
                    % Transmisibility Matrix
                    obj.Thph{ph,1+f} = Grid.ConnectivityMatrix * spdiags(Grid.Trans .* Mob_rho_h_Upwind, 0, N_Faces, N_Faces) * Grid.ConnectivityMatrix';

                    % Gravity Matrix
                    obj.Ghph{ph,1+f} = Grid.ConnectivityMatrix * spdiags(Grid.Trans .* Mob_rho_h_Upwind .* obj.GravityModel.RhoInt{ph, 1+f}, 0, N_Faces, N_Faces) * Grid.ConnectivityMatrix';

                case('cartesian_grid')
                    Nx = Grid.Nx;
                    Ny = Grid.Ny;
                    Nz = Grid.Nz;
                    N  = Grid.N;

                    %% Th Matrix
                    Tx = zeros(Nx+1, Ny, Nz);
                    Ty = zeros(Nx, Ny+1, Nz);
                    Tz = zeros(Nx, Ny, Nz+1);

                    % Apply upwind operator; Phase enthalpy is included in Upwind
                    Mupx = obj.UpWind{ph,1+f}.x * (obj.Mob(Index.Start:Index.End,ph) .* rho(:,ph) .* h(:,ph));
                    Mupy = obj.UpWind{ph,1+f}.y * (obj.Mob(Index.Start:Index.End,ph) .* rho(:,ph) .* h(:,ph));
                    Mupz = obj.UpWind{ph,1+f}.z * (obj.Mob(Index.Start:Index.End,ph) .* rho(:,ph) .* h(:,ph));
                    Mupx = reshape(Mupx, Nx, Ny, Nz);
                    Mupy = reshape(Mupy, Nx, Ny, Nz);
                    Mupz = reshape(Mupz, Nx, Ny, Nz);

                    % Transmissibility matrix
                    Tx(2:Nx,:,:)= Grid.Tx(2:Nx,:,:).*Mupx(1:Nx-1,:,:);
                    Ty(:,2:Ny,:)= Grid.Ty(:,2:Ny,:).*Mupy(:,1:Ny-1,:);
                    Tz(:,:,2:Nz)= Grid.Tz(:,:,2:Nz).*Mupz(:,:,1:Nz-1);

                    %Construct matrix
                    x1 = reshape(Tx(1:Grid.Nx,:,:), Grid.N, 1);   % West boundary,  transmissbiility between I and IW   (x: W->E, y: S->N, z: Bot->Top)
                    x2 = reshape(Tx(2:Grid.Nx+1,:,:), Grid.N, 1); % East boundary,  transmissbiility between I and IE
                    y1 = reshape(Ty(:,1:Grid.Ny,:), Grid.N, 1);   % South boundary, transmissbiility between I and IS
                    y2 = reshape(Ty(:,2:Grid.Ny+1,:), Grid.N, 1); % North boundary, transmissbiility between I and IN
                    z1 = reshape(Tz(:,:,1:Grid.Nz), Grid.N, 1);   % Bottom boundary,transmissbiility between I and IBot
                    z2 = reshape(Tz(:,:,2:Grid.Nz+1), Grid.N, 1); % Top boundary,   transmissbiility between I and ITop

                    DiagVecs = [-z2,-y2,-x2,z2+y2+x2+y1+x1+z1,-x1,-y1,-z1];
                    DiagIndx = [-Nx*Ny,-Nx,-1,0,1,Nx,Nx*Ny];
                    obj.Thph{ph,1+f} = spdiags(DiagVecs,DiagIndx,N,N);

                    % Gravity Matrix
                    Tx(2:Grid.Nx,:,:)= Tx(2:Grid.Nx,:,:) .* obj.GravityModel.RhoInt{ph,1+f}.x(2:Grid.Nx,:,:);
                    Ty(:,2:Grid.Ny,:)= Ty(:,2:Grid.Ny,:) .* obj.GravityModel.RhoInt{ph,1+f}.y(:,2:Grid.Ny,:);
                    Tz(:,:,2:Grid.Nz)= Tz(:,:,2:Grid.Nz) .* obj.GravityModel.RhoInt{ph,1+f}.z(:,:,2:Grid.Nz);

                    %Construct matrix
                    x1 = reshape(Tx(1:Grid.Nx,:,:), Grid.N, 1);
                    x2 = reshape(Tx(2:Grid.Nx+1,:,:), Grid.N, 1);
                    y1 = reshape(Ty(:,1:Grid.Ny,:), Grid.N, 1);
                    y2 = reshape(Ty(:,2:Grid.Ny+1,:), Grid.N, 1);
                    z1 = reshape(Tz(:,:,1:Grid.Nz), Grid.N, 1);
                    z2 = reshape(Tz(:,:,2:Grid.Nz+1), Grid.N, 1);

                    DiagVecs = [-z2,-y2,-x2,z2+y2+x2+y1+x1+z1,-x1,-y1,-z1];
                    DiagIndx = [-Grid.Nx*Grid.Ny, -Grid.Nx, -1, 0, 1, Grid.Nx, Grid.Nx*Grid.Ny];
                    obj.Ghph{ph,1+f} = spdiags(DiagVecs, DiagIndx, Grid.N, Grid.N);
            end
        end

        function ConductiveHeatTransmissibilityMatrix(obj, Grid, State, f)
            switch class(Grid)
                case('corner_point_grid')
                    N_Faces = length(Grid.Trans);

                    CondEff = State.Properties('CondEff').Value;
                    Grid.ComputeRockHeatConductivities(CondEff);
                    %% Heat Conduction Flux Transmissibility matrix
                    obj.Tk{1,1+f} = Grid.ConnectivityMatrix * spdiags(Grid.HeatTrans,0,N_Faces,N_Faces) * Grid.ConnectivityMatrix';

                case('cartesian_grid')
                    % Tk : Heat conduction transmissibility
                    Nx = Grid.Nx;
                    Ny = Grid.Ny;
                    Nz = Grid.Nz;

                    CondEff = State.Properties('CondEff').Value;
                    Grid.ComputeHeatConductivitiyTransmissibilities(CondEff);

                    %% Tk Matrix
                    condTrans_x = zeros(Nx+1, Ny, Nz);
                    condTrans_y = zeros(Nx, Ny+1, Nz);
                    condTrans_z = zeros(Nx, Ny, Nz+1);

                    condTrans_x(2:Nx,:,:)= Grid.condTrans_x(2:Nx,:,:);
                    condTrans_y(:,2:Ny,:)= Grid.condTrans_y(:,2:Ny,:);
                    condTrans_z(:,:,2:Nz)= Grid.condTrans_z(:,:,2:Nz);

                    x1 = reshape(condTrans_x(1:Grid.Nx,:,:), Grid.N, 1);
                    x2 = reshape(condTrans_x(2:Grid.Nx+1,:,:), Grid.N, 1);
                    y1 = reshape(condTrans_y(:,1:Grid.Ny,:), Grid.N, 1);
                    y2 = reshape(condTrans_y(:,2:Grid.Ny+1,:), Grid.N, 1);
                    z1 = reshape(condTrans_z(:,:,1:Grid.Nz), Grid.N, 1);
                    z2 = reshape(condTrans_z(:,:,2:Grid.Nz+1), Grid.N, 1);

                    DiagVecs = [-z2,-y2,-x2,z2+y2+x2+y1+x1+z1,-x1,-y1,-z1];
                    DiagIndx = [-Grid.Nx*Grid.Ny, -Grid.Nx, -1, 0, 1, Grid.Nx, Grid.Nx*Grid.Ny];
                    obj.Tk{1,1+f} = spdiags(DiagVecs, DiagIndx, Grid.N, Grid.N);
            end
        end

    end
end