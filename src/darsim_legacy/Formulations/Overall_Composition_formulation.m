% Overall Composition Formulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini, Yuhang Wang
%TU Delft
%Created: 13 December 2017
%Last modified: 26 July 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Overall_Composition_formulation < Compositional_formulation
    properties
        dxdz
        dxdp
        dSdp
        dSdz
        drhoTdp
        drhodz
        drhoTdz
        dMobdp
        dMobdz
        dPcdp
        dPcdz

        rho_old
        rho
        P
        z_old
        z
        x
    end
    methods
        function obj = Overall_Composition_formulation(SimulationInput,n_phases,n_components)
            obj@Compositional_formulation(n_components);

            % allocate memory
            if (SimulationInput.FracturesProperties.isFractured)
                N = SimulationInput.ReservoirProperties.Grid.N_ActiveCells + SimulationInput.FracturesProperties.NumOfAllFracGrids; % total number of cells
            else
                N = SimulationInput.ReservoirProperties.Grid.N_ActiveCells;
            end
            obj.rho_old = zeros(N, n_phases);
            obj.rho = zeros(N, n_phases);
            obj.P = zeros(N, n_phases);
            obj.z_old = zeros(N, n_components);
            obj.z = zeros(N, n_components);
            obj.x = zeros(N, n_components*n_phases);
        end
        function x = GetPrimaryUnknowns(obj, ProductionSystem, DiscretizationModel)
            Nt = DiscretizationModel.N;
            Nm = DiscretizationModel.ReservoirGrid.N;
            x = zeros(obj.NofComponents * Nt, 1);
            %% 1. reservoir
            Start = 1;
            End = Nm;
            x(Start:End) = ProductionSystem.Reservoir.State.Properties('P_2').Value;
            for i=1:obj.NofComponents-1
                Start = End + 1;
                End = Start + Nm - 1;
                x(Start:End) = ProductionSystem.Reservoir.State.Properties(['z_', num2str(i)]).Value;
            end
            %% 2. fracture
            if ProductionSystem.FracturesNetwork.Active
                Nf = DiscretizationModel.FracturesGrid.N;
                for f = 1:ProductionSystem.FracturesNetwork.NumOfFrac
                    Start = End + 1;
                    End = Start + Nf(f) - 1;
                    x(Start:End) = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('P_2').Value;
                    for i=1:obj.NofComponents-1
                        Start = End + 1;
                        End = Start + Nf(f) - 1;
                        x(Start:End) = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['z_', num2str(i)]).Value;
                    end
                end
            end
        end
        %%
        function ComputePropertiesAndDerivatives(obj, ProductionSystem, FluidModel, ~)
            %% 1. reservoir
            % This is the bitchy part! (Indeed)
            [obj.dxdp, obj.dxdz] = FluidModel.DxDpDz(ProductionSystem.Reservoir.State, ProductionSystem.Reservoir.State.Properties('SinglePhase').Value);
            % %%%%%%%%%%%%%%%%%%%%%%%%
            obj.drhodp = FluidModel.ComputeDrhoDp(ProductionSystem.Reservoir.State, ProductionSystem.Reservoir.State.Properties('SinglePhase').Value, obj.dxdp);
            obj.drhodz = FluidModel.ComputeDrhoDz(ProductionSystem.Reservoir.State, ProductionSystem.Reservoir.State.Properties('SinglePhase').Value, obj.dxdz);
            obj.dSdp = FluidModel.ComputeDSDp(ProductionSystem.Reservoir.State, obj.drhodp, obj.dxdp(:, end));
            obj.dSdz = FluidModel.ComputeDSDz(ProductionSystem.Reservoir.State, obj.drhodz, obj.dxdz(:, end, :));
            obj.drhoTdp = FluidModel.ComputeDrhotDp(ProductionSystem.Reservoir.State,obj.drhodp, obj.dSdp);
            obj.drhoTdz = FluidModel.ComputeDrhotDz(ProductionSystem.Reservoir.State, obj.drhodz, obj.dSdz);
            obj.dMobdp = FluidModel.ComputeDMobDp(ProductionSystem.Reservoir, [obj.dSdp, obj.dSdp]);
            obj.dMob = FluidModel.ComputeDMobDz(ProductionSystem.Reservoir, obj.dSdz, obj.dxdz);

            obj.dPcdz = FluidModel.ComputeDPcDz(ProductionSystem.Reservoir.State, obj.dSdz, 1);
            obj.dPcdp = FluidModel.ComputeDPcDp(ProductionSystem.Reservoir.State, obj.dSdp, 1);
            %% 2. fracture
            for f = 1:ProductionSystem.FracturesNetwork.NumOfFrac

                [dxdp_temp, dxdz_temp] = FluidModel.DxDpDz(ProductionSystem.FracturesNetwork.Fractures(f).State, ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('SinglePhase').Value);

                % %%%%%%%%%%%%%%%%%%%%%%%%
                % Mengjie:
                % obj.dxdp > dxdp_temp   & obj.dxdz > dxdz_temp
                drhodp_temp = FluidModel.ComputeDrhoDp(ProductionSystem.FracturesNetwork.Fractures(f).State, ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('SinglePhase').Value, dxdp_temp);
                drhodz_temp = FluidModel.ComputeDrhoDz(ProductionSystem.FracturesNetwork.Fractures(f).State, ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('SinglePhase').Value, dxdz_temp);
                dSdp_temp = FluidModel.ComputeDSDp(ProductionSystem.FracturesNetwork.Fractures(f).State, drhodp_temp, dxdp_temp(:, end));
                dSdz_temp = FluidModel.ComputeDSDz(ProductionSystem.FracturesNetwork.Fractures(f).State, drhodz_temp, dxdz_temp(:, end, :));

                obj.drhoTdz = [obj.drhoTdz; FluidModel.ComputeDrhotDz(ProductionSystem.FracturesNetwork.Fractures(f).State, drhodz_temp, dSdz_temp)];
                obj.drhoTdp = [obj.drhoTdp; FluidModel.ComputeDrhotDp(ProductionSystem.FracturesNetwork.Fractures(f).State, drhodp_temp, dSdp_temp)];
                obj.dMobdp = [obj.dMobdp; FluidModel.ComputeDMobDp(ProductionSystem.FracturesNetwork.Fractures(f), [dSdp_temp, dSdp_temp])];
                obj.dMob = [obj.dMob; FluidModel.ComputeDMobDz(ProductionSystem.FracturesNetwork.Fractures(f), dSdz_temp, dxdz_temp)];
                obj.dPcdz = [obj.dPcdz; FluidModel.ComputeDPcDz(ProductionSystem.FracturesNetwork.Fractures(f).State, dSdz_temp, 0)];
                obj.dPcdp = [obj.dPcdp; FluidModel.ComputeDPcDp(ProductionSystem.FracturesNetwork.Fractures(f).State, dSdp_temp, 0)];

                obj.dxdp = [obj.dxdp; dxdp_temp];
                obj.dxdz = [obj.dxdz; dxdz_temp];
                obj.dSdp = [obj.dSdp; dSdp_temp];
                obj.dSdz = [obj.dSdz; dSdz_temp];
                obj.drhodp = [obj.drhodp; drhodp_temp];
                obj.drhodz = [obj.drhodz; drhodz_temp];
            end

        end
        %%
        function [Residual, RHS] = BuildMediumResidual(obj, Medium, Grid, dt, State0, Index, qw, qf, f, c)
            N = Grid.N;
            pv = Medium.Por .* Grid.Volume;

            rhoT_old = State0.Properties('rhoT').Value(Index.Start:Index.End);
            rhoT =  Medium.State.Properties('rhoT').Value;

            % Accumulation term
            n = rhoT .* obj.z(Index.Start:Index.End,c);
            n_old = rhoT_old .* obj.z_old(Index.Start:Index.End,c);

            % A = speye(N).*pv/dt;
            % A = spdiags(pv/dt, 0, N, N);
            A = sparse(1:N, 1:N, pv/dt);

            % Residual
            RHS = qw(Index.Start:Index.End, c);
            %             Residual = A * (n - n_old) ... % accumulation
            %                      + obj.Tph{c, 1} *  obj.P(:,1) ... % Convective term
            %                      + obj.Tph{c, 2} *  obj.P(:,2) ...
            %                      - obj.Gph{c,1} * Grid.Depth ... % Gravity
            %                      - obj.Gph{c,2} * Grid.Depth ...
            %                      - qw(Index.Start:Index.End,c); % Source/Sink

            Residual = A * (n - n_old) ... % accumulation
                + obj.Tph{2*c-1, 1+f} *  obj.P(Index.Start:Index.End,1) ... % Convective term
                + obj.Tph{2*c, 1+f} *  obj.P(Index.Start:Index.End,2) ...
                - obj.Gph{2*c-1,1+f} * Grid.Depth ... % Gravity
                - obj.Gph{2*c,1+f} * Grid.Depth ...
                + obj.Dph{2*c-1, 1+f} * obj.x(Index.Start:Index.End, 2*c-1) ...
                + obj.Dph{2*c, 1+f} * obj.x(Index.Start:Index.End, 2*c) ...
                - qw(Index.Start:Index.End, c) ... % Source/Sink
                - qf(Index.Start:Index.End, c); % frac-matrix

        end
        %%
        function CopytoLocalVariables(obj, ProductionSystem, DiscretizationModel, State0)
            Index.Start = 1;
            Index.End = DiscretizationModel.ReservoirGrid.N;
            % Reservoir
            % Copy values to local variables (phase-dependent)
            for i = 1:obj.NofPhases
                obj.rho_old(Index.Start:Index.End, i) = State0.Properties(['rho_', num2str(i)]).Value(Index.Start:Index.End);
                obj.P(Index.Start:Index.End, i) = ProductionSystem.Reservoir.State.Properties(['P_', num2str(i)]).Value;
                obj.rho(Index.Start:Index.End, i) = ProductionSystem.Reservoir.State.Properties(['rho_', num2str(i)]).Value;
                obj.S(Index.Start:Index.End, i) = ProductionSystem.Reservoir.State.Properties(['S_', num2str(i)]).Value;
            end
            % Copy values to local variables (component-dependent)
            for i=1:obj.NofComponents
                obj.z_old(Index.Start:Index.End,i) = State0.Properties(['z_', num2str(i)]).Value(Index.Start:Index.End);
                obj.z(Index.Start:Index.End, i) = ProductionSystem.Reservoir.State.Properties(['z_', num2str(i)]).Value;
                for j=1:obj.NofPhases
                    obj.x(Index.Start:Index.End,(i-1)*obj.NofPhases + j) = ProductionSystem.Reservoir.State.Properties(['x_', num2str(i),'ph',num2str(j)]).Value;
                end
            end
            % Fractures
            for f = 1:ProductionSystem.FracturesNetwork.NumOfFrac
                Index.Start = Index.End+1;
                Index.End = Index.Start + DiscretizationModel.FracturesGrid.N(f) - 1;
                % Copy values to local variables (phase-dependent)
                for i = 1:obj.NofPhases
                    obj.rho_old(Index.Start:Index.End, i) = State0.Properties(['rho_', num2str(i)]).Value(Index.Start:Index.End);
                    obj.P(Index.Start:Index.End, i) = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['P_', num2str(i)]).Value;
                    obj.rho(Index.Start:Index.End, i) = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['rho_', num2str(i)]).Value;
                    obj.S(Index.Start:Index.End, i) = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['S_', num2str(i)]).Value;
                end
                % Copy values to local variables (component-dependent)
                for i=1:obj.NofComponents
                    obj.z_old(Index.Start:Index.End,i) = State0.Properties(['z_', num2str(i)]).Value(Index.Start:Index.End);
                    obj.z(Index.Start:Index.End, i) = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['z_', num2str(i)]).Value;
                    for j=1:obj.NofPhases
                        obj.x(Index.Start:Index.End,(i-1)*obj.NofPhases + j) = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['x_', num2str(i),'ph',num2str(j)]).Value;
                    end
                end
            end
        end
        %%
        function [Residual, RHS] = BuildFullResidual(obj, ProductionSystem, DiscretizationModel, Time, dt, State0)

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Initialise residual vector
            Nt = DiscretizationModel.N;
            Residual = zeros(Nt * obj.NofComponents, 1);
            RHS = zeros(Nt * obj.NofComponents, 1);
            Nm = DiscretizationModel.ReservoirGrid.N;

            % Copy values to local variables
            CopytoLocalVariables(obj, ProductionSystem, DiscretizationModel, State0);

            % Source terms
            % well
            qw = obj.ComputeSourceTerms(Nt, ProductionSystem.Wells);
            % This part is updated by including the wells schedule
            % if (Time <= obj.Injection_time) % well shutdown control
            %     qw = obj.ComputeSourceTerms(Nt, ProductionSystem.Wells);
            % else
            %     qw = zeros(Nt, obj.NofComponents);
            % end
            % fracture_matrix
            qf = zeros(Nt, obj.NofComponents);
            if ProductionSystem.FracturesNetwork.Active
                qf = ComputeSourceTerms_frac_mat(obj, ProductionSystem, DiscretizationModel);
            end

            % unique(qf)
            % Build residual
            for c = 1:obj.NofComponents
                % Reservoir
                % Index
                Index.Start = 1;
                Index.End = Nm;
                % Phase Transmissibilities
                obj.TransmissibilityMatrix(DiscretizationModel.ReservoirGrid, ProductionSystem.Reservoir, obj.rho(Index.Start:Index.End,:), obj.x(Index.Start:Index.End,(c-1)*2+1:(c-1)*2+2), Index, c, 0);
                % Residual
                [Residualm, RHSm] = BuildMediumResidual(obj, ProductionSystem.Reservoir, DiscretizationModel.ReservoirGrid, dt, State0, Index, qw, qf, 0, c);
                Residual((c-1)*Nt + Index.Start: (c-1)*Nt + Index.End) = Residualm;
                RHS(     (c-1)*Nt + Index.Start: (c-1)*Nt + Index.End) = RHSm;

                % Fractures
                for f = 1:ProductionSystem.FracturesNetwork.NumOfFrac
                    % Index
                    Index.Start = Index.End+1;
                    Index.End = Index.Start + DiscretizationModel.FracturesGrid.N(f) - 1;
                    % Phase Transmissibilities
                    obj.TransmissibilityMatrix(DiscretizationModel.FracturesGrid.Grids(f), ProductionSystem.FracturesNetwork.Fractures(f), obj.rho(Index.Start:Index.End,:), obj.x(Index.Start:Index.End,(c-1)*2+1:(c-1)*2+2), Index, c, f);
                    % Residual
                    Residual_frac_f = BuildMediumResidual(obj, ProductionSystem.FracturesNetwork.Fractures(f), DiscretizationModel.FracturesGrid.Grids(f), dt, State0, Index, qw, qf, f, c);
                    Residual((c-1)*Nt + Index.Start: (c-1)*Nt + Index.End) = Residual_frac_f;
                end
            end
        end
        %%
        % function ConstrainedPressureResidual(obj, Formulation, FluidModel, ProductionSystem, DiscretizationModel, dt, time, State0)
        %     % Initialise residual vector
        %     Nt = DiscretizationModel.N;
        %     Residual = zeros(Nt * obj.NofComponents, 1);
        %     RHS = zeros(Nt * obj.NofComponents, 1);
        %     Nm = DiscretizationModel.ReservoirGrid.N;
        % 
        %     % Copy values to local variables
        %     CopytoLocalVariables(obj, ProductionSystem, DiscretizationModel, State0);
        % 
        %     %% Source terms
        %     % well
        %     qw = zeros(Nt,2);
        %     if (time <= obj.Injection_time) % well shutdown control
        %         qw = obj.ComputeSourceTerms(Nt, ProductionSystem.Wells);
        %     else
        %         qw = zeros(Nt,2);
        %     end
        % 
        %     % fracture_matrix
        %     qf = zeros(Nt, obj.NofComponents);
        %     if ProductionSystem.FracturesNetwork.Active
        %         qf = ComputeSourceTerms_frac_mat(obj, ProductionSystem, DiscretizationModel);
        %     end
        % 
        % 
        %     for c = 1:obj.NofComponents
        %         %% Build residual
        %         % Reservoir
        %         % Index
        %         Index.Start = 1;
        %         Index.End = Nm;
        %         % Phase Transmissibilities
        %         obj.TransmissibilityMatrix(DiscretizationModel.ReservoirGrid, obj.rho(Index.Start:Index.End,:), obj.x(Index.Start:Index.End,(c-1)*2+1:(c-1)*2+2), Index, c, 0);
        %         % Residual
        %         [Residualm, RHSm] = BuildMediumResidual(obj, ProductionSystem.Reservoir, DiscretizationModel.ReservoirGrid, dt, State0, Index, qw, qf, 0, c);
        %         Residual((c-1)*Nt + Index.Start: (c-1)*Nt + Index.End) = Residualm;
        %         RHS(     (c-1)*Nt + Index.Start: (c-1)*Nt + Index.End) = RHSm;
        % 
        % 
        %         %% Jacobian's assembly
        %         % Reservoir
        %         Index.Start = 1;
        %         Index.End = Nm;
        %         [Jp_res, Jz_res] = BuildMediumJacobian(obj, ProductionSystem.Reservoir, DiscretizationModel.ReservoirGrid, dt, Index, 0, c);
        % 
        % 
        %         % Wells
        %         if (time <= obj.Injection_time) % well shutdown control
        %             [Jp_res, ~] = obj.AddWellsToJacobian(Jp_res, Jz_res, obj.P(:,obj.NofPhases), obj.rho, obj.x, ProductionSystem.Wells, ProductionSystem.Reservoir.K, c);
        %         end
        % 
        %         %% Solve for delta
        %         if c == 1
        %             deltaz = Jz_res\(-Residualm);
        %         else
        %             deltaP2 = Jp_res\(-Residualm);
        %         end
        %     end
        % 
        %     %% Update Reservoir State
        %     Pm = ProductionSystem.Reservoir.State.Properties('P_2');
        %     Pm.update(deltaP2(1:Nm));
        %     Pm = ProductionSystem.Reservoir.State.Properties('P_1');
        %     Pm.update(deltaP2(1:Nm));
        %     %             % 1.b Update z
        %     %             DeltaLast = zeros(Nm, 1);
        %     %             for c = 1:obj.NofComponents-1
        %     %                 Zm = ProductionSystem.Reservoir.State.Properties(['z_', num2str(c)]);
        %     %                 % Deltaz = delta(c*Nt + 1:c*Nt + Nm);
        %     %                 Deltaz = deltaz(1:Nm);
        %     %
        %     %                 Zm.update(Deltaz);
        %     %                 Zm.Value = max(Zm.Value, 1.0e-8);
        %     %                 Zm.Value = min(Zm.Value, 1-1.0e-8);
        %     %                 DeltaLast = DeltaLast + Deltaz;
        %     %             end
        %     %             Zm = ProductionSystem.Reservoir.State.Properties(['z_', num2str(obj.NofComponents)]);
        %     %             Zm.update(-DeltaLast);
        %     %             Zm.Value = max(Zm.Value, 1.0e-8);
        %     %             Zm.Value = min(Zm.Value, 1-1.0e-8);
        % 
        %     %% Update the peoperties and well state
        %     % Update Remaining properties (x, S, rho, rhoT, Pc)
        %     obj.UpdatePhaseCompositions(ProductionSystem.Reservoir, DiscretizationModel, FluidModel, 1);
        %     ProductionSystem.Wells.UpdateState(ProductionSystem, Formulation, FluidModel);
        % end
        %% corner_point_grid case is added to BuildMediumJacobian function
        function [Jp, Jz] = BuildMediumJacobian(obj, Medium, Grid, dt, Index, f, c)
            % Create local variables
            % Nx = Grid.Nx;
            % Ny = Grid.Ny;
            % Nz = Grid.Nz;
            N = Grid.N;
            % N_Face = length(obj.UpWind{c,1+f});
            N_Face = size(obj.UpWind{1,1+f}, 1);
            pv = Grid.Volume .* Medium.Por;

            rhoT = Medium.State.Properties('rhoT').Value;

            % BUILD FIM JACOBIAN BLOCK BY BLOCK
            % Jz = sparse(N, N * (obj.NofComponents-1));
            Jz_comp = cell(1, obj.NofComponents-1);
            % 1.a: divergence
            Jp = obj.Tph{2*c-1,1+f}  + obj.Tph{2*c, 1+f};
            switch class(Grid)
                case 'unstructured_grid'
                    nc = Grid.N;
                    nf = length(Grid.T);
                    C = [Grid.mesh.face.inner.upstream Grid.mesh.face.inner.downstream];
                    D1Ph1 = [ -double(obj.U{1,1+f}>=0)+double(obj.U{1,1+f}<0)  , double(obj.U{1,1+f}>=0)-double(obj.U{1,1+f}<0)]; D1Ph1(D1Ph1==1)=0; %D1 filled with pairs of -1 0 (U>=0) or 0 -1 (U<0)
                    D2Ph1 = [ -double(obj.U{1,1+f}>=0)+double(obj.U{1,1+f}<0)  , double(obj.U{1,1+f}>=0)-double(obj.U{1,1+f}<0)]; % D2 filled with pairs of -1 1 (U>=0) or 1 -1 (U<0)
                    D1Ph2 = [ -double(obj.U{2,1+f}>=0)+double(obj.U{2,1+f}<0)  , double(obj.U{2,1+f}>=0)-double(obj.U{2,1+f}<0)]; D1Ph2(D1Ph2==1)=0; %D1 filled with pairs of -1 0 (U>=0) or 0 -1 (U<0)
                    D2Ph2 = [ -double(obj.U{2,1+f}>=0)+double(obj.U{2,1+f}<0)  , double(obj.U{2,1+f}>=0)-double(obj.U{2,1+f}<0)]; % D2 filled with pairs of -1 1 (U>=0) or 1 -1 (U<0)
                    UpwindPermutation1Ph1 = sparse([(1:nf)'; (1:nf)'], C, D1Ph1, nf, nc)'; % nc x nf
                    UpwindPermutation2Ph1 = sparse([(1:nf)'; (1:nf)'], C, D2Ph1, nf, nc)';
                    UpwindPermutation1Ph2 = sparse([(1:nf)'; (1:nf)'], C, D1Ph2, nf, nc)';
                    UpwindPermutation2Ph2 = sparse([(1:nf)'; (1:nf)'], C, D2Ph2, nf, nc)';
                    % 1.b: compressibility part
                    dMupPh1 = obj.UpWind{1,1+f} * (obj.Mob(Index.Start:Index.End, 1) .* obj.x(Index.Start:Index.End,(c-1)*2+1) .* obj.drhodp(Index.Start:Index.End,1) + obj.Mob(Index.Start:Index.End, 1) .* obj.dxdp(Index.Start:Index.End,(c-1)*2+1) .* obj.rho(Index.Start:Index.End,1) + obj.dMobdp(Index.Start:Index.End, 1) .* obj.x(Index.Start:Index.End,(c-1)*2+1) .* obj.rho(Index.Start:Index.End,1));
                    dMupPh2 = obj.UpWind{2,1+f} * (obj.Mob(Index.Start:Index.End, 2) .* obj.x(Index.Start:Index.End,(c-1)*2+2) .* obj.drhodp(Index.Start:Index.End,2) + obj.Mob(Index.Start:Index.End, 2) .* obj.dxdp(Index.Start:Index.End,(c-1)*2+2) .* obj.rho(Index.Start:Index.End,2) + obj.dMobdp(Index.Start:Index.End, 2) .* obj.x(Index.Start:Index.End,(c-1)*2+2) .* obj.rho(Index.Start:Index.End,2));
                    FluxDerivativePh1 = dMupPh1 .* abs(obj.U{1,1+f}); % nf x 1
                    FluxDerivativePh2 = dMupPh2 .* abs(obj.U{2,1+f});
                    acc = pv/dt .* (obj.drhoTdp(Index.Start:Index.End) .* obj.z(Index.Start:Index.End,c));
                    Jp = Jp + sparse(1:N, 1:N, acc)...
                        + (UpwindPermutation1Ph1 * sparse(1:N_Face, 1:N_Face, FluxDerivativePh1) * UpwindPermutation2Ph1')'... % ncxnf * nfxnf * nfxnc = ncxnc
                        + (UpwindPermutation1Ph2 * sparse(1:N_Face, 1:N_Face, FluxDerivativePh2) * UpwindPermutation2Ph2')';

                    % Diffusion
                    CellNeighbor1Index = Grid.Neighbours(:,1);
                    CellNeighbor2Index = Grid.Neighbours(:,2);
                    CellNeighbor1Vec   = Grid.mesh.face.inner.centroid - Grid.mesh.element.centroid(CellNeighbor1Index,:);
                    CellNeighbor2Vec   = Grid.mesh.element.centroid(CellNeighbor2Index,:) - Grid.mesh.face.inner.centroid;
                    Nvec               = Grid.mesh.face.inner.normal;
                    
                    S_D_ph1 = obj.S(Index.Start:Index.End,1) .* Medium.Diff(:,1) .* Grid.Por;
                    S_D_ph2 = obj.S(Index.Start:Index.End,2) .* Medium.Diff(:,2) .* Grid.Por;
                    dS_rho_D_ph1 =  obj.dSdp(Index.Start:Index.End,1) .* Medium.Diff(:,1) .* Grid.Por;
                    dS_rho_D_ph2 = -obj.dSdp(Index.Start:Index.End,1) .* Medium.Diff(:,2) .* Grid.Por;
                    NeighborPermutation1 = sparse(repmat((1:nf)', 1, 2), C, [ones(nf,1), zeros(nf,1)], nf, nc)'; % nc x nf
                    NeighborPermutation2 = sparse(repmat((1:nf)', 1, 2), C, [zeros(nf,1), ones(nf,1)], nf, nc)'; % nc x nf

                    Half1_ph1 = S_D_ph1(CellNeighbor1Index) .* abs(sum(CellNeighbor1Vec .* Nvec, 2) ./ (sum(CellNeighbor1Vec .* CellNeighbor1Vec, 2) + 1e-30));
                    Half2_ph1 = S_D_ph1(CellNeighbor2Index) .* abs(sum(CellNeighbor2Vec .* Nvec, 2) ./ (sum(CellNeighbor2Vec .* CellNeighbor2Vec, 2) + 1e-30));
                    Half1_ph2 = S_D_ph2(CellNeighbor1Index) .* abs(sum(CellNeighbor1Vec .* Nvec, 2) ./ (sum(CellNeighbor1Vec .* CellNeighbor1Vec, 2) + 1e-30));
                    Half2_ph2 = S_D_ph2(CellNeighbor2Index) .* abs(sum(CellNeighbor2Vec .* Nvec, 2) ./ (sum(CellNeighbor2Vec .* CellNeighbor2Vec, 2) + 1e-30));

                    rho_1 = obj.rho(Index.Start:Index.End,1);
                    rho_2 = obj.rho(Index.Start:Index.End,2);
                    drhodp1 = obj.drhodp(Index.Start:Index.End,1);
                    drhodp2 = obj.drhodp(Index.Start:Index.End,2);
                    dDT_ph1_L = 0.5 * (rho_1(CellNeighbor1Index) + rho_1(CellNeighbor2Index)) .* abs(sum(CellNeighbor1Vec .* Nvec, 2) ./ (sum(CellNeighbor1Vec .* CellNeighbor1Vec, 2) + 1e-30)) .* Half2_ph1.^2 ./ (Half1_ph1 + Half2_ph1 + 1e-30).^2;
                    dDT_ph1_R = 0.5 * (rho_1(CellNeighbor1Index) + rho_1(CellNeighbor2Index)) .* abs(sum(CellNeighbor2Vec .* Nvec, 2) ./ (sum(CellNeighbor2Vec .* CellNeighbor2Vec, 2) + 1e-30)) .* Half1_ph1.^2 ./ (Half1_ph1 + Half2_ph1 + 1e-30).^2;
                    dDT_ph2_L = 0.5 * (rho_2(CellNeighbor1Index) + rho_2(CellNeighbor2Index)) .* abs(sum(CellNeighbor1Vec .* Nvec, 2) ./ (sum(CellNeighbor1Vec .* CellNeighbor1Vec, 2) + 1e-30)) .* Half2_ph2.^2 ./ (Half1_ph2 + Half2_ph2 + 1e-30).^2;
                    dDT_ph2_R = 0.5 * (rho_2(CellNeighbor1Index) + rho_2(CellNeighbor2Index)) .* abs(sum(CellNeighbor2Vec .* Nvec, 2) ./ (sum(CellNeighbor2Vec .* CellNeighbor2Vec, 2) + 1e-30)) .* Half1_ph2.^2 ./ (Half1_ph2 + Half2_ph2 + 1e-30).^2;
                    
                    vecL = obj.U_D{2*c-1,1+f} .* dDT_ph1_L .* dS_rho_D_ph1(CellNeighbor1Index) + obj.U_D{2*c,  1+f} .* dDT_ph2_L .* dS_rho_D_ph2(CellNeighbor1Index);
                    vecR = obj.U_D{2*c-1,1+f} .* dDT_ph1_R .* dS_rho_D_ph1(CellNeighbor2Index) + obj.U_D{2*c,  1+f} .* dDT_ph2_R .* dS_rho_D_ph2(CellNeighbor2Index);

                    Jp = Jp + (NeighborPermutation1 * spdiags(vecL, 0, nf, nf) * Grid.ConnectivityMatrix')' ...
                        + (NeighborPermutation2 * spdiags(vecR, 0, nf, nf) * Grid.ConnectivityMatrix')' ...
                        + obj.Dph{2*c-1,1+f} * spdiags(obj.dxdp(Index.Start:Index.End,2*c-1), 0, nc, nc) ...
                        + obj.Dph{2*c,  1+f} * spdiags(obj.dxdp(Index.Start:Index.End,2*c  ), 0, nc, nc);

                    vecL = obj.U_D{2*c-1,1+f} .* (drhodp1(CellNeighbor1Index) / 2) .* Grid.Diff(:,1) + obj.U_D{2*c,  1+f} .* (drhodp2(CellNeighbor1Index) / 2) .* Grid.Diff(:,2);
                    vecR = obj.U_D{2*c-1,1+f} .* (drhodp1(CellNeighbor2Index) / 2) .* Grid.Diff(:,1) + obj.U_D{2*c,  1+f} .* (drhodp2(CellNeighbor2Index) / 2) .* Grid.Diff(:,2);

                    Jp = Jp + (NeighborPermutation1 * spdiags(vecL, 0, nf, nf) * Grid.ConnectivityMatrix')' ...
                        + (NeighborPermutation2 * spdiags(vecR, 0, nf, nf) * Grid.ConnectivityMatrix')';  

                    % 2. Component c composition block
                    for j=1:obj.NofComponents-1
                        Jz_comp{j} = sparse(N, N);
                        dMupPh1 = obj.UpWind{1,1+f} * ( obj.dMob(Index.Start:Index.End, 1, j) .* obj.x(Index.Start:Index.End,(c-1)*2+1) .* obj.rho(Index.Start:Index.End,1) + obj.Mob(Index.Start:Index.End, 1) .* obj.dxdz(Index.Start:Index.End,(c-1)*2+1, j) .* obj.rho(Index.Start:Index.End,1) + obj.Mob(Index.Start:Index.End, 1) .* obj.x(Index.Start:Index.End, (c-1)*2+1) .* obj.drhodz(Index.Start:Index.End, 1, j));
                        dMupPh2 = obj.UpWind{2,1+f} * ( obj.dMob(Index.Start:Index.End, 2, j) .* obj.x(Index.Start:Index.End,(c-1)*2+2) .* obj.rho(Index.Start:Index.End,2) + obj.Mob(Index.Start:Index.End, 2) .* obj.dxdz(Index.Start:Index.End,(c-1)*2+2, j) .* obj.rho(Index.Start:Index.End,2) + obj.Mob(Index.Start:Index.End, 2) .* obj.x(Index.Start:Index.End, (c-1)*2+2) .* obj.drhodz(Index.Start:Index.End, 2, j));

                        FluxDerivativePh1 = dMupPh1 .* abs(obj.U{1,1+f});
                        FluxDerivativePh2 = dMupPh2 .* abs(obj.U{2,1+f});

                        if c == obj.NofComponents
                            acc = - pv/dt .* rhoT + pv/dt .* obj.z(Index.Start:Index.End,c) .* obj.drhoTdz(Index.Start:Index.End, j);
                        elseif c == j
                            acc = pv/dt .* rhoT + pv/dt .* obj.z(Index.Start:Index.End,c) .* obj.drhoTdz(Index.Start:Index.End, j);
                        else
                            acc = 0 + pv/dt .* obj.z(Index.Start:Index.End,c) .* obj.drhoTdz(Index.Start:Index.End, j);
                        end

                        Jz_comp{j} = sparse(1:N, 1:N, acc)...
                            + (UpwindPermutation1Ph1 * sparse(1:N_Face, 1:N_Face, FluxDerivativePh1) * UpwindPermutation2Ph1')'...
                            + (UpwindPermutation1Ph2 * sparse(1:N_Face, 1:N_Face, FluxDerivativePh2) * UpwindPermutation2Ph2')';

                                               % Diffusion
                        dS_rho_D_ph1 =  obj.dSdz(Index.Start:Index.End,j) .* Medium.Diff(:,1) .* Grid.Por;
                        dS_rho_D_ph2 = -obj.dSdz(Index.Start:Index.End,j) .* Medium.Diff(:,2) .* Grid.Por;
                        vecL = obj.U_D{2*c-1,1+f} .* dDT_ph1_L .* dS_rho_D_ph1(CellNeighbor1Index) + obj.U_D{2*c,  1+f} .* dDT_ph2_L .* dS_rho_D_ph2(CellNeighbor1Index);
                        vecR = obj.U_D{2*c-1,1+f} .* dDT_ph1_R .* dS_rho_D_ph1(CellNeighbor2Index) + obj.U_D{2*c,  1+f} .* dDT_ph2_R .* dS_rho_D_ph2(CellNeighbor2Index);

                        Jz_comp{j} = Jz_comp{j} + (NeighborPermutation1 * sparse(1:nf, 1:nf, vecL, nf, nf) * Grid.ConnectivityMatrix')' ...
                            + (NeighborPermutation2 * sparse(1:nf, 1:nf, vecR, nf, nf) * Grid.ConnectivityMatrix')' ...
                            + obj.Dph{2*c-1,1+f} * spdiags(obj.dxdz(Index.Start:Index.End,2*c-1,j), 0, N, N) ...
                            + obj.Dph{2*c,  1+f} * spdiags(obj.dxdz(Index.Start:Index.End,2*c,  j), 0, N, N);

                        drhodz1 = obj.drhodz(Index.Start:Index.End,1,j);
                        drhodz2 = obj.drhodz(Index.Start:Index.End,2,j);
                        vecL = obj.U_D{2*c-1,1+f} .* (drhodz1(CellNeighbor1Index) / 2) .* Grid.Diff(:,1) + obj.U_D{2*c,  1+f} .* (drhodz2(CellNeighbor1Index) / 2) .* Grid.Diff(:,2);
                        vecR = obj.U_D{2*c-1,1+f} .* (drhodz1(CellNeighbor2Index) / 2) .* Grid.Diff(:,1) + obj.U_D{2*c,  1+f} .* (drhodz2(CellNeighbor2Index) / 2) .* Grid.Diff(:,2);

                        Jz_comp{j} = Jz_comp{j} + (NeighborPermutation1 * sparse(1:nf, 1:nf, vecL, nf, nf) * Grid.ConnectivityMatrix')' ...
                            + (NeighborPermutation2 * sparse(1:nf, 1:nf, vecR, nf, nf) * Grid.ConnectivityMatrix')';
                        
                        % Capillarity
                        Jz_comp{j} = Jz_comp{j} - obj.Tph{2*c-1,1+f} * sparse(1:N, 1:N, obj.dPcdz(Index.Start:Index.End, j));
                        % Jz = Jz - obj.Tph{2*c-1,1+f} * spdiags(obj.dPcdz(Index.Start:Index.End), 0, N, N);
                    end
                case('corner_point_grid')
                    nc = Grid.N;
                    nf = length(Grid.Trans);
                    C = [ Grid.CornerPointGridData.Internal_Faces.CellNeighbor1Index , Grid.CornerPointGridData.Internal_Faces.CellNeighbor2Index ];
                    D1Ph1 = [ -double(obj.U{1,1+f}>=0)+double(obj.U{1,1+f}<0)  , double(obj.U{1,1+f}>=0)-double(obj.U{1,1+f}<0)]; D1Ph1(D1Ph1==1)=0; %D1 filled with pairs of -1 0 (U>=0) or 0 -1 (U<0)
                    D2Ph1 = [ -double(obj.U{1,1+f}>=0)+double(obj.U{1,1+f}<0)  , double(obj.U{1,1+f}>=0)-double(obj.U{1,1+f}<0)]; % D2 filled with pairs of -1 1 (U>=0) or 1 -1 (U<0)
                    D1Ph2 = [ -double(obj.U{2,1+f}>=0)+double(obj.U{2,1+f}<0)  , double(obj.U{2,1+f}>=0)-double(obj.U{2,1+f}<0)]; D1Ph2(D1Ph2==1)=0; %D1 filled with pairs of -1 0 (U>=0) or 0 -1 (U<0)
                    D2Ph2 = [ -double(obj.U{2,1+f}>=0)+double(obj.U{2,1+f}<0)  , double(obj.U{2,1+f}>=0)-double(obj.U{2,1+f}<0)]; % D2 filled with pairs of -1 1 (U>=0) or 1 -1 (U<0)

                    UpwindPermutation1Ph1 = sparse([(1:nf)'; (1:nf)'], C, D1Ph1, nf, nc)'; % nc x nf
                    UpwindPermutation2Ph1 = sparse([(1:nf)'; (1:nf)'], C, D2Ph1, nf, nc)';
                    UpwindPermutation1Ph2 = sparse([(1:nf)'; (1:nf)'], C, D1Ph2, nf, nc)';
                    UpwindPermutation2Ph2 = sparse([(1:nf)'; (1:nf)'], C, D2Ph2, nf, nc)';

                    % 1.b: compressibility part

                    dMupPh1 = obj.UpWind{1,1+f} * (obj.Mob(Index.Start:Index.End, 1) .* obj.x(Index.Start:Index.End,(c-1)*2+1) .* obj.drhodp(Index.Start:Index.End,1) + obj.Mob(Index.Start:Index.End, 1) .* obj.dxdp(Index.Start:Index.End,(c-1)*2+1) .* obj.rho(Index.Start:Index.End,1) + obj.dMobdp(Index.Start:Index.End, 1) .* obj.x(Index.Start:Index.End,(c-1)*2+1) .* obj.rho(Index.Start:Index.End,1));
                    dMupPh2 = obj.UpWind{2,1+f} * (obj.Mob(Index.Start:Index.End, 2) .* obj.x(Index.Start:Index.End,(c-1)*2+2) .* obj.drhodp(Index.Start:Index.End,2) + obj.Mob(Index.Start:Index.End, 2) .* obj.dxdp(Index.Start:Index.End,(c-1)*2+2) .* obj.rho(Index.Start:Index.End,2) + obj.dMobdp(Index.Start:Index.End, 2) .* obj.x(Index.Start:Index.End,(c-1)*2+2) .* obj.rho(Index.Start:Index.End,2));

                    FluxDerivativePh1 = dMupPh1 .* abs(obj.U{1,1+f}); % nf x 1
                    FluxDerivativePh2 = dMupPh2 .* abs(obj.U{2,1+f});

                    acc = pv/dt .* (obj.drhoTdp(Index.Start:Index.End) .* obj.z(Index.Start:Index.End,c));
                    Jp = Jp + sparse(1:N, 1:N, acc)...
                        + (UpwindPermutation1Ph1 * sparse(1:N_Face, 1:N_Face, FluxDerivativePh1) * UpwindPermutation2Ph1')'... % ncxnf * nfxnf * nfxnc = ncxnc
                        + (UpwindPermutation1Ph2 * sparse(1:N_Face, 1:N_Face, FluxDerivativePh2) * UpwindPermutation2Ph2')';

                    % Diffusion
                    CellNeighbor1Index = Grid.CornerPointGridData.Internal_Faces.CellNeighbor1Index;
                    CellNeighbor2Index = Grid.CornerPointGridData.Internal_Faces.CellNeighbor2Index;
                    CellNeighbor1Vec = Grid.CornerPointGridData.Internal_Faces.CellNeighbor1Vec;
                    CellNeighbor2Vec = Grid.CornerPointGridData.Internal_Faces.CellNeighbor2Vec;
                    Nvec = Grid.CornerPointGridData.Internal_Faces.Nvec;
                    
                    S_D_ph1 = obj.S(Index.Start:Index.End,1) .* Medium.Diff(:,1) .* Grid.Por;
                    S_D_ph2 = obj.S(Index.Start:Index.End,2) .* Medium.Diff(:,2) .* Grid.Por;
                    dS_rho_D_ph1 =  obj.dSdp(Index.Start:Index.End,1) .* Medium.Diff(:,1) .* Grid.Por;
                    dS_rho_D_ph2 = -obj.dSdp(Index.Start:Index.End,1) .* Medium.Diff(:,2) .* Grid.Por;
                    NeighborPermutation1 = sparse(repmat((1:nf)', 1, 2), C, [ones(nf,1), zeros(nf,1)], nf, nc)'; % nc x nf
                    NeighborPermutation2 = sparse(repmat((1:nf)', 1, 2), C, [zeros(nf,1), ones(nf,1)], nf, nc)'; % nc x nf

                    Half1_ph1 = S_D_ph1(CellNeighbor1Index) .* abs(sum(CellNeighbor1Vec .* Nvec, 2) ./ (sum(CellNeighbor1Vec .* CellNeighbor1Vec, 2) + 1e-30));
                    Half2_ph1 = S_D_ph1(CellNeighbor2Index) .* abs(sum(CellNeighbor2Vec .* Nvec, 2) ./ (sum(CellNeighbor2Vec .* CellNeighbor2Vec, 2) + 1e-30));
                    Half1_ph2 = S_D_ph2(CellNeighbor1Index) .* abs(sum(CellNeighbor1Vec .* Nvec, 2) ./ (sum(CellNeighbor1Vec .* CellNeighbor1Vec, 2) + 1e-30));
                    Half2_ph2 = S_D_ph2(CellNeighbor2Index) .* abs(sum(CellNeighbor2Vec .* Nvec, 2) ./ (sum(CellNeighbor2Vec .* CellNeighbor2Vec, 2) + 1e-30));

                    rho_1 = obj.rho(Index.Start:Index.End,1);
                    rho_2 = obj.rho(Index.Start:Index.End,2);
                    drhodp1 = obj.drhodp(Index.Start:Index.End,1);
                    drhodp2 = obj.drhodp(Index.Start:Index.End,2);
                    dDT_ph1_L = 0.5 * (rho_1(CellNeighbor1Index) + rho_1(CellNeighbor2Index)) .* abs(sum(CellNeighbor1Vec .* Nvec, 2) ./ (sum(CellNeighbor1Vec .* CellNeighbor1Vec, 2) + 1e-30)) .* Half2_ph1.^2 ./ (Half1_ph1 + Half2_ph1 + 1e-30).^2;
                    dDT_ph1_R = 0.5 * (rho_1(CellNeighbor1Index) + rho_1(CellNeighbor2Index)) .* abs(sum(CellNeighbor2Vec .* Nvec, 2) ./ (sum(CellNeighbor2Vec .* CellNeighbor2Vec, 2) + 1e-30)) .* Half1_ph1.^2 ./ (Half1_ph1 + Half2_ph1 + 1e-30).^2;
                    dDT_ph2_L = 0.5 * (rho_2(CellNeighbor1Index) + rho_2(CellNeighbor2Index)) .* abs(sum(CellNeighbor1Vec .* Nvec, 2) ./ (sum(CellNeighbor1Vec .* CellNeighbor1Vec, 2) + 1e-30)) .* Half2_ph2.^2 ./ (Half1_ph2 + Half2_ph2 + 1e-30).^2;
                    dDT_ph2_R = 0.5 * (rho_2(CellNeighbor1Index) + rho_2(CellNeighbor2Index)) .* abs(sum(CellNeighbor2Vec .* Nvec, 2) ./ (sum(CellNeighbor2Vec .* CellNeighbor2Vec, 2) + 1e-30)) .* Half1_ph2.^2 ./ (Half1_ph2 + Half2_ph2 + 1e-30).^2;
                    
                    vecL = obj.U_D{2*c-1,1+f} .* dDT_ph1_L .* dS_rho_D_ph1(CellNeighbor1Index) + obj.U_D{2*c,  1+f} .* dDT_ph2_L .* dS_rho_D_ph2(CellNeighbor1Index);
                    vecR = obj.U_D{2*c-1,1+f} .* dDT_ph1_R .* dS_rho_D_ph1(CellNeighbor2Index) + obj.U_D{2*c,  1+f} .* dDT_ph2_R .* dS_rho_D_ph2(CellNeighbor2Index);

                    Jp = Jp + (NeighborPermutation1 * spdiags(vecL, 0, nf, nf) * Grid.ConnectivityMatrix')' ...
                        + (NeighborPermutation2 * spdiags(vecR, 0, nf, nf) * Grid.ConnectivityMatrix')' ...
                        + obj.Dph{2*c-1,1+f} * spdiags(obj.dxdp(Index.Start:Index.End,2*c-1), 0, nc, nc) ...
                        + obj.Dph{2*c,  1+f} * spdiags(obj.dxdp(Index.Start:Index.End,2*c  ), 0, nc, nc);

                    vecL = obj.U_D{2*c-1,1+f} .* (drhodp1(CellNeighbor1Index) / 2) .* Grid.Diff(:,1) + obj.U_D{2*c,  1+f} .* (drhodp2(CellNeighbor1Index) / 2) .* Grid.Diff(:,2);
                    vecR = obj.U_D{2*c-1,1+f} .* (drhodp1(CellNeighbor2Index) / 2) .* Grid.Diff(:,1) + obj.U_D{2*c,  1+f} .* (drhodp2(CellNeighbor2Index) / 2) .* Grid.Diff(:,2);

                    Jp = Jp + (NeighborPermutation1 * spdiags(vecL, 0, nf, nf) * Grid.ConnectivityMatrix')' ...
                        + (NeighborPermutation2 * spdiags(vecR, 0, nf, nf) * Grid.ConnectivityMatrix')';    

                    % 2. Component c composition block
                    for j=1:obj.NofComponents-1
                        Jz_comp{j} = sparse(N, N);
                        dMupPh1 = obj.UpWind{1,1+f} * ( obj.dMob(Index.Start:Index.End, 1, j) .* obj.x(Index.Start:Index.End,(c-1)*2+1) .* obj.rho(Index.Start:Index.End,1) + obj.Mob(Index.Start:Index.End, 1) .* obj.dxdz(Index.Start:Index.End,(c-1)*2+1, j) .* obj.rho(Index.Start:Index.End,1) + obj.Mob(Index.Start:Index.End, 1) .* obj.x(Index.Start:Index.End, (c-1)*2+1) .* obj.drhodz(Index.Start:Index.End, 1, j));
                        dMupPh2 = obj.UpWind{2,1+f} * ( obj.dMob(Index.Start:Index.End, 2, j) .* obj.x(Index.Start:Index.End,(c-1)*2+2) .* obj.rho(Index.Start:Index.End,2) + obj.Mob(Index.Start:Index.End, 2) .* obj.dxdz(Index.Start:Index.End,(c-1)*2+2, j) .* obj.rho(Index.Start:Index.End,2) + obj.Mob(Index.Start:Index.End, 2) .* obj.x(Index.Start:Index.End, (c-1)*2+2) .* obj.drhodz(Index.Start:Index.End, 2, j));

                        FluxDerivativePh1 = dMupPh1 .* abs(obj.U{1,1+f});
                        FluxDerivativePh2 = dMupPh2 .* abs(obj.U{2,1+f});

                        if c == obj.NofComponents
                            acc = - pv/dt .* rhoT + pv/dt .* obj.z(Index.Start:Index.End,c) .* obj.drhoTdz(Index.Start:Index.End, j);
                        elseif c == j
                            acc = pv/dt .* rhoT + pv/dt .* obj.z(Index.Start:Index.End,c) .* obj.drhoTdz(Index.Start:Index.End, j);
                        else
                            acc = 0 + pv/dt .* obj.z(Index.Start:Index.End,c) .* obj.drhoTdz(Index.Start:Index.End, j);
                        end

                        Jz_comp{j} = sparse(1:N, 1:N, acc) + (UpwindPermutation1Ph1 * sparse(1:N_Face, 1:N_Face, FluxDerivativePh1) * UpwindPermutation2Ph1')'...
                            + (UpwindPermutation1Ph2 * sparse(1:N_Face, 1:N_Face, FluxDerivativePh2) * UpwindPermutation2Ph2')';

                        % Diffusion
                        dS_rho_D_ph1 =  obj.dSdz(Index.Start:Index.End,j) .* Medium.Diff(:,1) .* Grid.Por;
                        dS_rho_D_ph2 = -obj.dSdz(Index.Start:Index.End,j) .* Medium.Diff(:,2) .* Grid.Por;
                        vecL = obj.U_D{2*c-1,1+f} .* dDT_ph1_L .* dS_rho_D_ph1(CellNeighbor1Index) + obj.U_D{2*c,  1+f} .* dDT_ph2_L .* dS_rho_D_ph2(CellNeighbor1Index);
                        vecR = obj.U_D{2*c-1,1+f} .* dDT_ph1_R .* dS_rho_D_ph1(CellNeighbor2Index) + obj.U_D{2*c,  1+f} .* dDT_ph2_R .* dS_rho_D_ph2(CellNeighbor2Index);

                        Jz_comp{j} = Jz_comp{j} + (NeighborPermutation1 * sparse(1:nf, 1:nf, vecL, nf, nf) * Grid.ConnectivityMatrix')' ...
                            + (NeighborPermutation2 * sparse(1:nf, 1:nf, vecR, nf, nf) * Grid.ConnectivityMatrix')' ...
                            + obj.Dph{2*c-1,1+f} * spdiags(obj.dxdz(Index.Start:Index.End,2*c-1,j), 0, N, N) ...
                            + obj.Dph{2*c,  1+f} * spdiags(obj.dxdz(Index.Start:Index.End,2*c,  j), 0, N, N);

                        drhodz1 = obj.drhodz(Index.Start:Index.End,1,j);
                        drhodz2 = obj.drhodz(Index.Start:Index.End,2,j);
                        vecL = obj.U_D{2*c-1,1+f} .* (drhodz1(CellNeighbor1Index) / 2) .* Grid.Diff(:,1) + obj.U_D{2*c,  1+f} .* (drhodz2(CellNeighbor1Index) / 2) .* Grid.Diff(:,2);
                        vecR = obj.U_D{2*c-1,1+f} .* (drhodz1(CellNeighbor2Index) / 2) .* Grid.Diff(:,1) + obj.U_D{2*c,  1+f} .* (drhodz2(CellNeighbor2Index) / 2) .* Grid.Diff(:,2);

                        Jz_comp{j} = Jz_comp{j} + (NeighborPermutation1 * sparse(1:nf, 1:nf, vecL, nf, nf) * Grid.ConnectivityMatrix')' ...
                            + (NeighborPermutation2 * sparse(1:nf, 1:nf, vecR, nf, nf) * Grid.ConnectivityMatrix')';

                        % Capillarity
                        Jz_comp{j} = Jz_comp{j} - obj.Tph{2*c-1,1+f} * sparse(1:N, 1:N, obj.dPcdz(Index.Start:Index.End, j));
                        % Jz = Jz - obj.Tph{2*c-1,1+f} * spdiags(obj.dPcdz(Index.Start:Index.End), 0, N, N);
                    end

                case('cartesian_grid')
                    Nx = Grid.Nx;
                    Ny = Grid.Ny;
                    Nz = Grid.Nz;

                    dDTx_ph1_L = zeros(Grid.Nx, Grid.Ny, Grid.Nz);
                    dDTy_ph1_L = zeros(Grid.Nx, Grid.Ny, Grid.Nz);
                    dDTz_ph1_L = zeros(Grid.Nx, Grid.Ny, Grid.Nz);
                    dDTx_ph2_L = zeros(Grid.Nx, Grid.Ny, Grid.Nz);
                    dDTy_ph2_L = zeros(Grid.Nx, Grid.Ny, Grid.Nz);
                    dDTz_ph2_L = zeros(Grid.Nx, Grid.Ny, Grid.Nz);
                    dDTx_ph1_R = zeros(Grid.Nx, Grid.Ny, Grid.Nz);
                    dDTy_ph1_R = zeros(Grid.Nx, Grid.Ny, Grid.Nz);
                    dDTz_ph1_R = zeros(Grid.Nx, Grid.Ny, Grid.Nz);
                    dDTx_ph2_R = zeros(Grid.Nx, Grid.Ny, Grid.Nz);
                    dDTy_ph2_R = zeros(Grid.Nx, Grid.Ny, Grid.Nz);
                    dDTz_ph2_R = zeros(Grid.Nx, Grid.Ny, Grid.Nz);

                    Tx_ph1 = zeros(Nx+1, Ny, Nz);
                    Ty_ph1 = zeros(Nx, Ny+1, Nz);
                    Tz_ph1 = zeros(Nx, Ny, Nz+1);
                    Tx_ph2 = zeros(Nx+1, Ny, Nz);
                    Ty_ph2 = zeros(Nx, Ny+1, Nz);
                    Tz_ph2 = zeros(Nx, Ny, Nz+1);

                    % 1.b: compressibility part
                    dMupxPh1 = reshape(obj.UpWind{1,1+f}.x * (obj.Mob(Index.Start:Index.End, 1) .* obj.x(Index.Start:Index.End,(c-1)*2+1) .* obj.drhodp(Index.Start:Index.End,1) + obj.Mob(Index.Start:Index.End, 1) .* obj.dxdp(Index.Start:Index.End,(c-1)*2+1) .* obj.rho(Index.Start:Index.End,1) + obj.dMobdp(Index.Start:Index.End, 1) .* obj.x(Index.Start:Index.End,(c-1)*2+1) .* obj.rho(Index.Start:Index.End,1)), Nx, Ny, Nz);
                    dMupyPh1 = reshape(obj.UpWind{1,1+f}.y * (obj.Mob(Index.Start:Index.End, 1) .* obj.x(Index.Start:Index.End,(c-1)*2+1) .* obj.drhodp(Index.Start:Index.End,1) + obj.Mob(Index.Start:Index.End, 1) .* obj.dxdp(Index.Start:Index.End,(c-1)*2+1) .* obj.rho(Index.Start:Index.End,1) + obj.dMobdp(Index.Start:Index.End, 1) .* obj.x(Index.Start:Index.End,(c-1)*2+1) .* obj.rho(Index.Start:Index.End,1)), Nx, Ny, Nz);
                    dMupzPh1 = reshape(obj.UpWind{1,1+f}.z * (obj.Mob(Index.Start:Index.End, 1) .* obj.x(Index.Start:Index.End,(c-1)*2+1) .* obj.drhodp(Index.Start:Index.End,1) + obj.Mob(Index.Start:Index.End, 1) .* obj.dxdp(Index.Start:Index.End,(c-1)*2+1) .* obj.rho(Index.Start:Index.End,1) + obj.dMobdp(Index.Start:Index.End, 1) .* obj.x(Index.Start:Index.End,(c-1)*2+1) .* obj.rho(Index.Start:Index.End,1)), Nx, Ny, Nz);
                    dMupxPh2 = reshape(obj.UpWind{2,1+f}.x * (obj.Mob(Index.Start:Index.End, 2) .* obj.x(Index.Start:Index.End,(c-1)*2+2) .* obj.drhodp(Index.Start:Index.End,2) + obj.Mob(Index.Start:Index.End, 2) .* obj.dxdp(Index.Start:Index.End,(c-1)*2+2) .* obj.rho(Index.Start:Index.End,2) + obj.dMobdp(Index.Start:Index.End, 2) .* obj.x(Index.Start:Index.End,(c-1)*2+2) .* obj.rho(Index.Start:Index.End,2)), Nx, Ny, Nz);
                    dMupyPh2 = reshape(obj.UpWind{2,1+f}.y * (obj.Mob(Index.Start:Index.End, 2) .* obj.x(Index.Start:Index.End,(c-1)*2+2) .* obj.drhodp(Index.Start:Index.End,2) + obj.Mob(Index.Start:Index.End, 2) .* obj.dxdp(Index.Start:Index.End,(c-1)*2+2) .* obj.rho(Index.Start:Index.End,2) + obj.dMobdp(Index.Start:Index.End, 2) .* obj.x(Index.Start:Index.End,(c-1)*2+2) .* obj.rho(Index.Start:Index.End,2)), Nx, Ny, Nz);
                    dMupzPh2 = reshape(obj.UpWind{2,1+f}.z * (obj.Mob(Index.Start:Index.End, 2) .* obj.x(Index.Start:Index.End,(c-1)*2+2) .* obj.drhodp(Index.Start:Index.End,2) + obj.Mob(Index.Start:Index.End, 2) .* obj.dxdp(Index.Start:Index.End,(c-1)*2+2) .* obj.rho(Index.Start:Index.End,2) + obj.dMobdp(Index.Start:Index.End, 2) .* obj.x(Index.Start:Index.End,(c-1)*2+2) .* obj.rho(Index.Start:Index.End,2)), Nx, Ny, Nz);

                    Tx_ph1(2:Grid.Nx,:,:)= dMupxPh1(1:Grid.Nx-1,:,:); % multiplying geometric part of transmissibility over inner interfaces. Nx in Mupx is for right boundary.
                    Ty_ph1(:,2:Grid.Ny,:)= dMupyPh1(:,1:Grid.Ny-1,:);
                    Tz_ph1(:,:,2:Grid.Nz)= dMupzPh1(:,:,1:Grid.Nz-1);
                    Tx_ph2(2:Grid.Nx,:,:)= dMupxPh2(1:Grid.Nx-1,:,:); % multiplying geometric part of transmissibility over inner interfaces. Nx in Mupx is for right boundary.
                    Ty_ph2(:,2:Grid.Ny,:)= dMupyPh2(:,1:Grid.Ny-1,:);
                    Tz_ph2(:,:,2:Grid.Nz)= dMupzPh2(:,:,1:Grid.Nz-1);

                    vecX1 = min(reshape(obj.U{1,1+f}.x(1:Nx,:,:),N,1), 0).* reshape(Tx_ph1(1:Nx,:,:),N,1) + min(reshape(obj.U{2,1+f}.x(1:Nx,:,:),N,1), 0).* reshape(Tx_ph2(1:Nx,:,:),N,1);
                    vecX2 = max(reshape(obj.U{1,1+f}.x(2:Nx+1,:,:),N,1), 0).* reshape(Tx_ph1(2:Nx+1,:,:),N,1) + max(reshape(obj.U{2,1+f}.x(2:Nx+1,:,:),N,1), 0).*reshape(Tx_ph2(2:Nx+1,:,:),N,1);
                    vecY1 = min(reshape(obj.U{1,1+f}.y(:,1:Ny,:),N,1), 0).* reshape(Ty_ph1(:,1:Ny,:),N,1) + min(reshape(obj.U{2,1+f}.y(:,1:Ny,:),N,1), 0).*reshape(Ty_ph2(:,1:Ny,:),N,1);
                    vecY2 = max(reshape(obj.U{1,1+f}.y(:,2:Ny+1,:),N,1), 0).* reshape(Ty_ph1(:,2:Ny+1,:),N,1) + max(reshape(obj.U{2,1+f}.y(:,2:Ny+1,:),N,1), 0).*reshape(Ty_ph2(:,2:Ny+1,:),N,1);
                    vecZ1 = min(reshape(obj.U{1,1+f}.z(:,:,1:Nz),N,1), 0).* reshape(Tz_ph1(:,:,1:Nz),N,1) + min(reshape(obj.U{2,1+f}.z(:,:,1:Nz),N,1), 0).*reshape(Tz_ph2(:,:,1:Nz),N,1);
                    vecZ2 = max(reshape(obj.U{1,1+f}.z(:,:,2:Nz+1),N,1), 0).* reshape(Tz_ph1(:,:,2:Nz+1),N,1) + max(reshape(obj.U{2,1+f}.z(:,:,2:Nz+1),N,1), 0).*reshape(Tz_ph2(:,:,2:Nz+1),N,1);
                    acc = pv/dt .* (obj.drhoTdp(Index.Start:Index.End) .* obj.z(Index.Start:Index.End,c));
                    DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2 + vecY2+vecX2-vecY1-vecX1-vecZ1+acc, vecX1, vecY1, vecZ1];
                    DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
                    Jp = Jp + spdiags(DiagVecs, DiagIndx, N, N);

                    S_D_ph1 = reshape(obj.S(Index.Start:Index.End,1) .* Medium.Diff(:,1) .* Grid.Por, Grid.Nx, Grid.Ny, Grid.Nz);
                    S_D_ph2 = reshape(obj.S(Index.Start:Index.End,2) .* Medium.Diff(:,2) .* Grid.Por, Grid.Nx, Grid.Ny, Grid.Nz);

                    dS_rho_D_ph1_x = obj.dSdp(Index.Start:Index.End,1)  .* Medium.Diff(:,1) .* Grid.Por .* Grid.Ax / Grid.dx;
                    dS_rho_D_ph1_y = obj.dSdp(Index.Start:Index.End,1)  .* Medium.Diff(:,1) .* Grid.Por .* Grid.Ay / Grid.dy;
                    dS_rho_D_ph1_z = obj.dSdp(Index.Start:Index.End,1)  .* Medium.Diff(:,1) .* Grid.Por .* Grid.Az / Grid.dz;
                    dS_rho_D_ph2_x = -obj.dSdp(Index.Start:Index.End,1) .* Medium.Diff(:,2) .* Grid.Por .* Grid.Ax / Grid.dx;
                    dS_rho_D_ph2_y = -obj.dSdp(Index.Start:Index.End,1) .* Medium.Diff(:,2) .* Grid.Por .* Grid.Ay / Grid.dy;
                    dS_rho_D_ph2_z = -obj.dSdp(Index.Start:Index.End,1) .* Medium.Diff(:,2) .* Grid.Por .* Grid.Az / Grid.dz;

                    rho_1 = reshape(obj.rho(Index.Start:Index.End,1), Grid.Nx, Grid.Ny, Grid.Nz);
                    rho_2 = reshape(obj.rho(Index.Start:Index.End,2), Grid.Nx, Grid.Ny, Grid.Nz);

                    dDTx_ph1_L(2:Nx,:,:)   = (rho_1(1:Nx-1,:,:) + rho_1(2:Nx,:,:)) / 2 .* 2 .* S_D_ph1(1:Nx-1,:,:) .^ 2 ./ (S_D_ph1(1:Nx-1,:,:) + S_D_ph1(2:Nx,:,:) + 1e-10).^2;
                    dDTx_ph1_R(1:Nx-1,:,:) = (rho_1(1:Nx-1,:,:) + rho_1(2:Nx,:,:)) / 2 .* 2 .* S_D_ph1(2:Nx,:,:)   .^ 2 ./ (S_D_ph1(1:Nx-1,:,:) + S_D_ph1(2:Nx,:,:) + 1e-10).^2;
                    dDTy_ph1_L(:,2:Ny,:)   = (rho_1(:,1:Ny-1,:) + rho_1(:,2:Ny,:)) / 2 .* 2 .* S_D_ph1(:,1:Ny-1,:) .^ 2 ./ (S_D_ph1(:,1:Ny-1,:) + S_D_ph1(:,2:Ny,:) + 1e-10).^2;
                    dDTy_ph1_R(:,1:Ny-1,:) = (rho_1(:,1:Ny-1,:) + rho_1(:,2:Ny,:)) / 2 .* 2 .* S_D_ph1(:,2:Ny,:)   .^ 2 ./ (S_D_ph1(:,1:Ny-1,:) + S_D_ph1(:,2:Ny,:) + 1e-10).^2;
                    dDTz_ph1_L(:,:,2:Nz)   = (rho_1(:,:,1:Nz-1) + rho_1(:,:,2:Nz)) / 2 .* 2 .* S_D_ph1(:,:,1:Nz-1) .^ 2 ./ (S_D_ph1(:,:,1:Nz-1) + S_D_ph1(:,:,2:Nz) + 1e-10).^2;
                    dDTz_ph1_R(:,:,1:Nz-1) = (rho_1(:,:,1:Nz-1) + rho_1(:,:,2:Nz)) / 2 .* 2 .* S_D_ph1(:,:,2:Nz)   .^ 2 ./ (S_D_ph1(:,:,1:Nz-1) + S_D_ph1(:,:,2:Nz) + 1e-10).^2;

                    dDTx_ph2_L(2:Nx,:,:)   = (rho_2(1:Nx-1,:,:) + rho_2(2:Nx,:,:)) / 2 .* 2 .* S_D_ph2(1:Nx-1,:,:) .^ 2 ./ (S_D_ph2(1:Nx-1,:,:) + S_D_ph2(2:Nx,:,:) + 1e-10).^2;
                    dDTx_ph2_R(1:Nx-1,:,:) = (rho_2(1:Nx-1,:,:) + rho_2(2:Nx,:,:)) / 2 .* 2 .* S_D_ph2(2:Nx,:,:)   .^ 2 ./ (S_D_ph2(1:Nx-1,:,:) + S_D_ph2(2:Nx,:,:) + 1e-10).^2;
                    dDTy_ph2_L(:,2:Ny,:)   = (rho_2(:,1:Ny-1,:) + rho_2(:,2:Ny,:)) / 2 .* 2 .* S_D_ph2(:,1:Ny-1,:) .^ 2 ./ (S_D_ph2(:,1:Ny-1,:) + S_D_ph2(:,2:Ny,:) + 1e-10).^2;
                    dDTy_ph2_R(:,1:Ny-1,:) = (rho_2(:,1:Ny-1,:) + rho_2(:,2:Ny,:)) / 2 .* 2 .* S_D_ph2(:,2:Ny,:)   .^ 2 ./ (S_D_ph2(:,1:Ny-1,:) + S_D_ph2(:,2:Ny,:) + 1e-10).^2;
                    dDTz_ph2_L(:,:,2:Nz)   = (rho_2(:,:,1:Nz-1) + rho_2(:,:,2:Nz)) / 2 .* 2 .* S_D_ph2(:,:,1:Nz-1) .^ 2 ./ (S_D_ph2(:,:,1:Nz-1) + S_D_ph2(:,:,2:Nz) + 1e-10).^2;
                    dDTz_ph2_R(:,:,1:Nz-1) = (rho_2(:,:,1:Nz-1) + rho_2(:,:,2:Nz)) / 2 .* 2 .* S_D_ph2(:,:,2:Nz)   .^ 2 ./ (S_D_ph2(:,:,1:Nz-1) + S_D_ph2(:,:,2:Nz) + 1e-10).^2;
                    
                    vecX1 = reshape(obj.U_D{2*c-1,1+f}.x(1:Nx,:,:),N,1)   .* dS_rho_D_ph1_x .* reshape(dDTx_ph1_L,N,1) + reshape(obj.U_D{2*c,1+f}.x(1:Nx,:,:),N,1)     .* dS_rho_D_ph2_x .* reshape(dDTx_ph2_L,N,1);
                    vecX2 = reshape(obj.U_D{2*c-1,1+f}.x(2:Nx+1,:,:),N,1) .* dS_rho_D_ph1_x .* reshape(dDTx_ph1_R,N,1) + reshape(obj.U_D{2*c,1+f}.x(2:Nx+1,:,:),N,1)   .* dS_rho_D_ph2_x .* reshape(dDTx_ph2_R,N,1);
                    vecY1 = reshape(obj.U_D{2*c-1,1+f}.y(:,1:Ny,:),N,1)   .* dS_rho_D_ph1_y .* reshape(dDTy_ph1_L,N,1) + reshape(obj.U_D{2*c,1+f}.y(:,1:Ny,:),N,1)     .* dS_rho_D_ph2_y .* reshape(dDTy_ph2_L,N,1);
                    vecY2 = reshape(obj.U_D{2*c-1,1+f}.y(:,2:Ny+1,:),N,1) .* dS_rho_D_ph1_y .* reshape(dDTy_ph1_R,N,1) + reshape(obj.U_D{2*c,1+f}.y(:,2:Ny+1,:),N,1)   .* dS_rho_D_ph2_y .* reshape(dDTy_ph2_R,N,1);
                    vecZ1 = reshape(obj.U_D{2*c-1,1+f}.z(:,:,1:Nz),N,1)   .* dS_rho_D_ph1_z .* reshape(dDTz_ph1_L,N,1) + reshape(obj.U_D{2*c,1+f}.z(:,:,1:Nz),N,1)     .* dS_rho_D_ph2_z .* reshape(dDTz_ph2_L,N,1);
                    vecZ2 = reshape(obj.U_D{2*c-1,1+f}.z(:,:,2:Nz+1),N,1) .* dS_rho_D_ph1_z .* reshape(dDTz_ph1_R,N,1) + reshape(obj.U_D{2*c,1+f}.z(:,:,2:Nz+1),N,1)   .* dS_rho_D_ph2_z .* reshape(dDTz_ph2_R,N,1);
                    
                    DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2 + vecY2 + vecX2 - vecY1 - vecX1 - vecZ1, vecX1, vecY1, vecZ1];
                    DiagIndx = [-Grid.Nx*Grid.Ny, -Grid.Nx, -1, 0, 1, Grid.Nx, Grid.Nx*Grid.Ny];
                    Jp = Jp + spdiags(DiagVecs, DiagIndx, Grid.N, Grid.N) + obj.Dph{2*c-1,1+f} * spdiags(obj.dxdp(Index.Start:Index.End,2*c-1), 0, N, N) + obj.Dph{2*c,1+f} * spdiags(obj.dxdp(Index.Start:Index.End,2*c), 0, N, N);
                    
                    TotTx_ph1 = Grid.Diffx(:,:,:,1) + Grid.DispT{1,1+f}.x;
                    TotTy_ph1 = Grid.Diffy(:,:,:,1) + Grid.DispT{1,1+f}.y;
                    TotTz_ph1 = Grid.Diffz(:,:,:,1) + Grid.DispT{1,1+f}.z;

                    TotTx_ph2 = Grid.Diffx(:,:,:,2) + Grid.DispT{2,1+f}.x;
                    TotTy_ph2 = Grid.Diffy(:,:,:,2) + Grid.DispT{2,1+f}.y;
                    TotTz_ph2 = Grid.Diffz(:,:,:,2) + Grid.DispT{2,1+f}.z;

                    vecX1 = reshape(obj.U_D{2*c-1,1+f}.x(1:Nx,:,:),N,1).* obj.drhodp(Index.Start:Index.End,1) / 2 .* reshape(TotTx_ph1(1:Nx,:,:),N,1) + reshape(obj.U_D{2*c,1+f}.x(1:Nx,:,:),N,1) .* obj.drhodp(Index.Start:Index.End,2) / 2  .* reshape(TotTx_ph2(1:Nx,:,:),N,1);
                    vecX2 = reshape(obj.U_D{2*c-1,1+f}.x(2:Nx+1,:,:),N,1).* obj.drhodp(Index.Start:Index.End,1) / 2 .* reshape(TotTx_ph1(2:Nx+1,:,:),N,1) + reshape(obj.U_D{2*c,1+f}.x(2:Nx+1,:,:),N,1) .* obj.drhodp(Index.Start:Index.End,2) / 2  .* reshape(TotTx_ph2(2:Nx+1,:,:),N,1);
                    vecY1 = reshape(obj.U_D{2*c-1,1+f}.y(:,1:Ny,:),N,1).* obj.drhodp(Index.Start:Index.End,1) / 2.* reshape(TotTy_ph1(:,1:Ny,:),N,1) + reshape(obj.U_D{2*c,1+f}.y(:,1:Ny,:),N,1) .* obj.drhodp(Index.Start:Index.End,2) / 2  .* reshape(TotTy_ph2(:,1:Ny,:),N,1);
                    vecY2 = reshape(obj.U_D{2*c-1,1+f}.y(:,2:Ny+1,:),N,1).* obj.drhodp(Index.Start:Index.End,1) / 2 .* reshape(TotTy_ph1(:,2:Ny+1,:),N,1) + reshape(obj.U_D{2*c,1+f}.y(:,2:Ny+1,:),N,1) .* obj.drhodp(Index.Start:Index.End,2) / 2  .* reshape(TotTy_ph2(:,2:Ny+1,:),N,1);
                    vecZ1 = reshape(obj.U_D{2*c-1,1+f}.z(:,:,1:Nz),N,1) .* obj.drhodp(Index.Start:Index.End,1) / 2 .* reshape(TotTz_ph1(:,:,1:Nz),N,1) + reshape(obj.U_D{2*c,1+f}.z(:,:,1:Nz),N,1) .* obj.drhodp(Index.Start:Index.End,2) / 2  .* reshape(TotTz_ph2(:,:,1:Nz),N,1);
                    vecZ2 = reshape(obj.U_D{2*c-1,1+f}.z(:,:,2:Nz+1),N,1) .* obj.drhodp(Index.Start:Index.End,1) / 2 .* reshape(TotTz_ph1(:,:,2:Nz+1),N,1) + reshape(obj.U_D{2*c,1+f}.z(:,:,2:Nz+1),N,1) .* obj.drhodp(Index.Start:Index.End,2) / 2  .* reshape(TotTz_ph2(:,:,2:Nz+1),N,1);
                    
                    DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2 + vecY2 + vecX2 - vecY1 - vecX1 - vecZ1, vecX1, vecY1, vecZ1];
                    DiagIndx = [-Grid.Nx*Grid.Ny, -Grid.Nx, -1, 0, 1, Grid.Nx, Grid.Nx*Grid.Ny];
                    Jp = Jp + spdiags(DiagVecs, DiagIndx, Grid.N, Grid.N);

                    % 2. Component c composition block
                    for j=1:obj.NofComponents-1
                        Jz_comp{j} = sparse(N, N);
                        dMupxPh1 = reshape(obj.UpWind{1,1+f}.x*...
                            ( obj.dMob(Index.Start:Index.End, 1, j) .* obj.x(Index.Start:Index.End,(c-1)*2+1) .* obj.rho(Index.Start:Index.End,1)...
                            + obj.Mob(Index.Start:Index.End, 1) .* obj.dxdz(Index.Start:Index.End,(c-1)*2+1, j) .* obj.rho(Index.Start:Index.End,1)...
                            + obj.Mob(Index.Start:Index.End, 1) .* obj.x(Index.Start:Index.End, (c-1)*2+1) .* obj.drhodz(Index.Start:Index.End, 1, j)), Nx, Ny, Nz);
                        dMupyPh1 = reshape(obj.UpWind{1,1+f}.y*...
                            ( obj.dMob(Index.Start:Index.End, 1, j) .* obj.x(Index.Start:Index.End, (c-1)*2+1) .* obj.rho(Index.Start:Index.End,1)...
                            + obj.Mob(Index.Start:Index.End, 1) .* obj.dxdz(Index.Start:Index.End, (c-1)*2+1, j) .* obj.rho(Index.Start:Index.End,1)...
                            + obj.Mob(Index.Start:Index.End, 1) .* obj.x(Index.Start:Index.End,(c-1)*2+1) .* obj.drhodz(Index.Start:Index.End, 1, j)), Nx, Ny, Nz);
                        dMupzPh1 = reshape(obj.UpWind{1,1+f}.z*...
                            ( obj.dMob(Index.Start:Index.End, 1, j) .* obj.x(Index.Start:Index.End,(c-1)*2+1) .* obj.rho(Index.Start:Index.End,1)...
                            + obj.Mob(Index.Start:Index.End,1) .* obj.dxdz(Index.Start:Index.End,(c-1)*2+1, j) .* obj.rho(Index.Start:Index.End,1)...
                            + obj.Mob(Index.Start:Index.End,1) .* obj.x(Index.Start:Index.End,(c-1)*2+1) .* obj.drhodz(Index.Start:Index.End, 1, j)), Nx, Ny, Nz);
                        dMupxPh2 = reshape(obj.UpWind{2,1+f}.x*...
                            ( obj.dMob(Index.Start:Index.End, 2, j) .* obj.x(Index.Start:Index.End,(c-1)*2+2) .* obj.rho(Index.Start:Index.End,2)...
                            + obj.Mob(Index.Start:Index.End,2) .* obj.dxdz(Index.Start:Index.End,(c-1)*2+2, j) .* obj.rho(Index.Start:Index.End,2)...
                            + obj.Mob(Index.Start:Index.End,2) .* obj.x(Index.Start:Index.End,(c-1)*2+2) .* obj.drhodz(Index.Start:Index.End, 2, j)), Nx, Ny, Nz);
                        dMupyPh2 = reshape(obj.UpWind{2,1+f}.y*...
                            ( obj.dMob(Index.Start:Index.End, 2, j) .* obj.x(Index.Start:Index.End,(c-1)*2+2) .* obj.rho(Index.Start:Index.End,2)...
                            + obj.Mob(Index.Start:Index.End, 2) .* obj.dxdz(Index.Start:Index.End,(c-1)*2+2, j) .* obj.rho(Index.Start:Index.End,2)...
                            + obj.Mob(Index.Start:Index.End, 2) .* obj.x(Index.Start:Index.End, (c-1)*2+2) .* obj.drhodz(Index.Start:Index.End, 2, j)), Nx, Ny, Nz);
                        dMupzPh2 = reshape(obj.UpWind{2,1+f}.z*...
                            ( obj.dMob(Index.Start:Index.End, 2, j) .* obj.x(Index.Start:Index.End, (c-1)*2+2) .* obj.rho(Index.Start:Index.End,2)...
                            + obj.Mob(Index.Start:Index.End,2) .* obj.dxdz(Index.Start:Index.End, (c-1)*2+2, j) .* obj.rho(Index.Start:Index.End,2)...
                            + obj.Mob(Index.Start:Index.End,2) .* obj.x(Index.Start:Index.End, (c-1)*2+2) .* obj.drhodz(Index.Start:Index.End, 2, j)), Nx, Ny, Nz);

                        Tx_ph1(2:Grid.Nx,:,:)= dMupxPh1(1:Grid.Nx-1,:,:); % multiplying geometric part of transmissibility over inner interfaces. Nx in Mupx is for right boundary.
                        Ty_ph1(:,2:Grid.Ny,:)= dMupyPh1(:,1:Grid.Ny-1,:);
                        Tz_ph1(:,:,2:Grid.Nz)= dMupzPh1(:,:,1:Grid.Nz-1);
                        Tx_ph2(2:Grid.Nx,:,:)= dMupxPh2(1:Grid.Nx-1,:,:); % multiplying geometric part of transmissibility over inner interfaces. Nx in Mupx is for right boundary.
                        Ty_ph2(:,2:Grid.Ny,:)= dMupyPh2(:,1:Grid.Ny-1,:);
                        Tz_ph2(:,:,2:Grid.Nz)= dMupzPh2(:,:,1:Grid.Nz-1);

                        vecX1 = min(reshape(obj.U{1,1+f}.x(1:Nx,:,:),N,1), 0).* reshape(Tx_ph1(1:Nx,:,:),N,1) + min(reshape(obj.U{2,1+f}.x(1:Nx,:,:),N,1), 0).* reshape(Tx_ph2(1:Nx,:,:),N,1);
                        vecX2 = max(reshape(obj.U{1,1+f}.x(2:Nx+1,:,:),N,1), 0).* reshape(Tx_ph1(2:Nx+1,:,:),N,1) + max(reshape(obj.U{2,1+f}.x(2:Nx+1,:,:),N,1), 0).*reshape(Tx_ph2(2:Nx+1,:,:),N,1);
                        vecY1 = min(reshape(obj.U{1,1+f}.y(:,1:Ny,:),N,1), 0).* reshape(Ty_ph1(:,1:Ny,:),N,1) + min(reshape(obj.U{2,1+f}.y(:,1:Ny,:),N,1), 0).*reshape(Ty_ph2(:,1:Ny,:),N,1);
                        vecY2 = max(reshape(obj.U{1,1+f}.y(:,2:Ny+1,:),N,1), 0).* reshape(Ty_ph1(:,2:Ny+1,:),N,1) + max(reshape(obj.U{2,1+f}.y(:,2:Ny+1,:),N,1), 0).*reshape(Ty_ph2(:,2:Ny+1,:),N,1);
                        vecZ1 = min(reshape(obj.U{1,1+f}.z(:,:,1:Nz),N,1), 0).* reshape(Tz_ph1(:,:,1:Nz),N,1) + min(reshape(obj.U{2,1+f}.z(:,:,1:Nz),N,1), 0).*reshape(Tz_ph2(:,:,1:Nz),N,1);
                        vecZ2 = max(reshape(obj.U{1,1+f}.z(:,:,2:Nz+1),N,1), 0).* reshape(Tz_ph1(:,:,2:Nz+1),N,1) + max(reshape(obj.U{2,1+f}.z(:,:,2:Nz+1),N,1), 0).*reshape(Tz_ph2(:,:,2:Nz+1),N,1);
                        if c == obj.NofComponents
                            acc = - pv/dt .* rhoT + pv/dt .* obj.z(Index.Start:Index.End,c) .* obj.drhoTdz(Index.Start:Index.End, j);
                        elseif c == j
                            acc = pv/dt .* rhoT + pv/dt .* obj.z(Index.Start:Index.End,c) .* obj.drhoTdz(Index.Start:Index.End, j);
                        else
                            acc = 0 + pv/dt .* obj.z(Index.Start:Index.End,c) .* obj.drhoTdz(Index.Start:Index.End, j);
                        end

                        DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2+vecY2+vecX2-vecZ1-vecY1-vecX1+acc, vecX1, vecY1, vecZ1];
                        DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
                        Jz_comp{j} = spdiags(DiagVecs,DiagIndx, N, N);

                        dS_rho_D_ph1_x = (obj.dSdz(Index.Start:Index.End, j) .* Medium.Diff(:,1) .* Grid.Por) .* Grid.Ax / Grid.dx;
                        dS_rho_D_ph1_y = (obj.dSdz(Index.Start:Index.End, j) .* Medium.Diff(:,1) .* Grid.Por) .* Grid.Ay / Grid.dy;
                        dS_rho_D_ph1_z = (obj.dSdz(Index.Start:Index.End, j) .* Medium.Diff(:,1) .* Grid.Por) .* Grid.Az / Grid.dz;
                        dS_rho_D_ph2_x = (-obj.dSdz(Index.Start:Index.End, j) .* Medium.Diff(:,2) .* Grid.Por) .* Grid.Ax / Grid.dx;
                        dS_rho_D_ph2_y = (-obj.dSdz(Index.Start:Index.End, j) .* Medium.Diff(:,2) .* Grid.Por) .* Grid.Ay / Grid.dy;
                        dS_rho_D_ph2_z = (-obj.dSdz(Index.Start:Index.End, j) .* Medium.Diff(:,2) .* Grid.Por) .* Grid.Az / Grid.dz;

                        vecX1 = reshape(obj.U_D{2*c-1,1+f}.x(1:Nx,:,:),N,1).* dS_rho_D_ph1_x .* reshape(dDTx_ph1_L,N,1) + reshape(obj.U_D{2*c,1+f}.x(1:Nx,:,:),N,1) .* dS_rho_D_ph2_x .* reshape(dDTx_ph2_L,N,1);
                        vecX2 = reshape(obj.U_D{2*c-1,1+f}.x(2:Nx+1,:,:),N,1).* dS_rho_D_ph1_x .* reshape(dDTx_ph1_R,N,1) + reshape(obj.U_D{2*c,1+f}.x(2:Nx+1,:,:),N,1) .* dS_rho_D_ph2_x .* reshape(dDTx_ph2_R,N,1);
                        vecY1 = reshape(obj.U_D{2*c-1,1+f}.y(:,1:Ny,:),N,1).* dS_rho_D_ph1_y .* reshape(dDTy_ph1_L,N,1) + reshape(obj.U_D{2*c,1+f}.y(:,1:Ny,:),N,1) .* dS_rho_D_ph2_y .* reshape(dDTy_ph2_L,N,1);
                        vecY2 = reshape(obj.U_D{2*c-1,1+f}.y(:,2:Ny+1,:),N,1).* dS_rho_D_ph1_y .* reshape(dDTy_ph1_R,N,1) + reshape(obj.U_D{2*c,1+f}.y(:,2:Ny+1,:),N,1) .* dS_rho_D_ph2_y .* reshape(dDTy_ph2_R,N,1);
                        vecZ1 = reshape(obj.U_D{2*c-1,1+f}.z(:,:,1:Nz),N,1) .* dS_rho_D_ph1_z .* reshape(dDTz_ph1_L,N,1) + reshape(obj.U_D{2*c,1+f}.z(:,:,1:Nz),N,1) .* dS_rho_D_ph2_z .* reshape(dDTz_ph2_L,N,1);
                        vecZ2 = reshape(obj.U_D{2*c-1,1+f}.z(:,:,2:Nz+1),N,1) .* dS_rho_D_ph1_z .* reshape(dDTz_ph1_R,N,1) + reshape(obj.U_D{2*c,1+f}.z(:,:,2:Nz+1),N,1) .* dS_rho_D_ph2_z .* reshape(dDTz_ph2_R,N,1);
                        
                        DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2 + vecY2 + vecX2 - vecY1 - vecX1 - vecZ1, vecX1, vecY1, vecZ1];
                        DiagIndx = [-Grid.Nx*Grid.Ny, -Grid.Nx, -1, 0, 1, Grid.Nx, Grid.Nx*Grid.Ny];
                        Jz_comp{j} = Jz_comp{j} + spdiags(DiagVecs, DiagIndx, Grid.N, Grid.N) + obj.Dph{2*c-1,1+f} * spdiags(obj.dxdz(Index.Start:Index.End,2*c-1,j), 0, N, N) + obj.Dph{2*c,1+f} * spdiags(obj.dxdz(Index.Start:Index.End,2*c,j), 0, N, N);

                        vecX1 = reshape(obj.U_D{2*c-1,1+f}.x(1:Nx,:,:),N,1).* obj.drhodz(Index.Start:Index.End,1,j) / 2 .* reshape(TotTx_ph1(1:Nx,:,:),N,1) + reshape(obj.U_D{2*c,1+f}.x(1:Nx,:,:),N,1) .* obj.drhodz(Index.Start:Index.End,2,j) / 2  .* reshape(TotTx_ph2(1:Nx,:,:),N,1);
                        vecX2 = reshape(obj.U_D{2*c-1,1+f}.x(2:Nx+1,:,:),N,1).* obj.drhodz(Index.Start:Index.End,1,j) / 2 .* reshape(TotTx_ph1(2:Nx+1,:,:),N,1) + reshape(obj.U_D{2*c,1+f}.x(2:Nx+1,:,:),N,1) .* obj.drhodz(Index.Start:Index.End,2,j) / 2  .* reshape(TotTx_ph2(2:Nx+1,:,:),N,1);
                        vecY1 = reshape(obj.U_D{2*c-1,1+f}.y(:,1:Ny,:),N,1).* obj.drhodz(Index.Start:Index.End,1,j) / 2.* reshape(TotTy_ph1(:,1:Ny,:),N,1) + reshape(obj.U_D{2*c,1+f}.y(:,1:Ny,:),N,1) .* obj.drhodz(Index.Start:Index.End,2,j) / 2  .* reshape(TotTy_ph2(:,1:Ny,:),N,1);
                        vecY2 = reshape(obj.U_D{2*c-1,1+f}.y(:,2:Ny+1,:),N,1).* obj.drhodz(Index.Start:Index.End,1,j) / 2 .* reshape(TotTy_ph1(:,2:Ny+1,:),N,1) + reshape(obj.U_D{2*c,1+f}.y(:,2:Ny+1,:),N,1) .* obj.drhodz(Index.Start:Index.End,2,j) / 2  .* reshape(TotTy_ph2(:,2:Ny+1,:),N,1);
                        vecZ1 = reshape(obj.U_D{2*c-1,1+f}.z(:,:,1:Nz),N,1) .* obj.drhodz(Index.Start:Index.End,1,j) / 2 .* reshape(TotTz_ph1(:,:,1:Nz),N,1) + reshape(obj.U_D{2*c,1+f}.z(:,:,1:Nz),N,1) .* obj.drhodz(Index.Start:Index.End,2,j) / 2  .* reshape(TotTz_ph2(:,:,1:Nz),N,1);
                        vecZ2 = reshape(obj.U_D{2*c-1,1+f}.z(:,:,2:Nz+1),N,1) .* obj.drhodz(Index.Start:Index.End,1,j) / 2 .* reshape(TotTz_ph1(:,:,2:Nz+1),N,1) + reshape(obj.U_D{2*c,1+f}.z(:,:,2:Nz+1),N,1) .* obj.drhodz(Index.Start:Index.End,2,j) / 2  .* reshape(TotTz_ph2(:,:,2:Nz+1),N,1);

                        DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2 + vecY2 + vecX2 - vecY1 - vecX1 - vecZ1, vecX1, vecY1, vecZ1];
                        DiagIndx = [-Grid.Nx*Grid.Ny, -Grid.Nx, -1, 0, 1, Grid.Nx, Grid.Nx*Grid.Ny];
                        Jz_comp{j} = Jz_comp{j} + spdiags(DiagVecs, DiagIndx, Grid.N, Grid.N);

                        % Capillarity
                        % Jz{i,j} = Jz{i,j} - obj.Tph{i,1} * spdiags(obj.dPcdz, 0, N, N);
                        Jz_comp{j} = Jz_comp{j} - obj.Tph{2*c-1,1+f} * spdiags(obj.dPcdz(Index.Start:Index.End, j), 0, N, N);
                    end
                case 'unstructured_fracture_grid'
                    nc = Grid.N;
                    nf = length(Grid.T);
                    C = [Grid.fracmesh.faces.inner.upstream Grid.fracmesh.faces.inner.downstream];
                    D1Ph1 = [ -double(obj.U{1,1+f}>=0)+double(obj.U{1,1+f}<0)  , double(obj.U{1,1+f}>=0)-double(obj.U{1,1+f}<0)]; D1Ph1(D1Ph1==1)=0; %D1 filled with pairs of -1 0 (U>=0) or 0 -1 (U<0)
                    D2Ph1 = [ -double(obj.U{1,1+f}>=0)+double(obj.U{1,1+f}<0)  , double(obj.U{1,1+f}>=0)-double(obj.U{1,1+f}<0)]; % D2 filled with pairs of -1 1 (U>=0) or 1 -1 (U<0)
                    D1Ph2 = [ -double(obj.U{2,1+f}>=0)+double(obj.U{2,1+f}<0)  , double(obj.U{2,1+f}>=0)-double(obj.U{2,1+f}<0)]; D1Ph2(D1Ph2==1)=0; %D1 filled with pairs of -1 0 (U>=0) or 0 -1 (U<0)
                    D2Ph2 = [ -double(obj.U{2,1+f}>=0)+double(obj.U{2,1+f}<0)  , double(obj.U{2,1+f}>=0)-double(obj.U{2,1+f}<0)]; % D2 filled with pairs of -1 1 (U>=0) or 1 -1 (U<0)
                    UpwindPermutation1Ph1 = sparse([(1:nf)'; (1:nf)'], C, D1Ph1, nf, nc)'; % nc x nf
                    UpwindPermutation2Ph1 = sparse([(1:nf)'; (1:nf)'], C, D2Ph1, nf, nc)';
                    UpwindPermutation1Ph2 = sparse([(1:nf)'; (1:nf)'], C, D1Ph2, nf, nc)';
                    UpwindPermutation2Ph2 = sparse([(1:nf)'; (1:nf)'], C, D2Ph2, nf, nc)';
                    % 1.b: compressibility part
                    dMupPh1 = obj.UpWind{1,1+f} * (obj.Mob(Index.Start:Index.End, 1) .* obj.x(Index.Start:Index.End,(c-1)*2+1) .* obj.drhodp(Index.Start:Index.End,1) + obj.Mob(Index.Start:Index.End, 1) .* obj.dxdp(Index.Start:Index.End,(c-1)*2+1) .* obj.rho(Index.Start:Index.End,1) + obj.dMobdp(Index.Start:Index.End, 1) .* obj.x(Index.Start:Index.End,(c-1)*2+1) .* obj.rho(Index.Start:Index.End,1));
                    dMupPh2 = obj.UpWind{2,1+f} * (obj.Mob(Index.Start:Index.End, 2) .* obj.x(Index.Start:Index.End,(c-1)*2+2) .* obj.drhodp(Index.Start:Index.End,2) + obj.Mob(Index.Start:Index.End, 2) .* obj.dxdp(Index.Start:Index.End,(c-1)*2+2) .* obj.rho(Index.Start:Index.End,2) + obj.dMobdp(Index.Start:Index.End, 2) .* obj.x(Index.Start:Index.End,(c-1)*2+2) .* obj.rho(Index.Start:Index.End,2));
                    FluxDerivativePh1 = dMupPh1 .* abs(obj.U{1,1+f}); % nf x 1
                    FluxDerivativePh2 = dMupPh2 .* abs(obj.U{2,1+f});
                    acc = pv/dt .* (obj.drhoTdp(Index.Start:Index.End) .* obj.z(Index.Start:Index.End,c));
                    Jp = Jp + sparse(1:N, 1:N, acc)...
                        + (UpwindPermutation1Ph1 * sparse(1:N_Face, 1:N_Face, FluxDerivativePh1) * UpwindPermutation2Ph1')'... % ncxnf * nfxnf * nfxnc = ncxnc
                        + (UpwindPermutation1Ph2 * sparse(1:N_Face, 1:N_Face, FluxDerivativePh2) * UpwindPermutation2Ph2')';
                    % 2. Component c composition block
                    for j=1:obj.NofComponents-1
                        Jz_comp{j} = sparse(N, N);
                        dMupPh1 = obj.UpWind{1,1+f} * ( obj.dMob(Index.Start:Index.End, 1, j) .* obj.x(Index.Start:Index.End,(c-1)*2+1) .* obj.rho(Index.Start:Index.End,1) + obj.Mob(Index.Start:Index.End, 1) .* obj.dxdz(Index.Start:Index.End,(c-1)*2+1, j) .* obj.rho(Index.Start:Index.End,1) + obj.Mob(Index.Start:Index.End, 1) .* obj.x(Index.Start:Index.End, (c-1)*2+1) .* obj.drhodz(Index.Start:Index.End, 1, j));
                        dMupPh2 = obj.UpWind{2,1+f} * ( obj.dMob(Index.Start:Index.End, 2, j) .* obj.x(Index.Start:Index.End,(c-1)*2+2) .* obj.rho(Index.Start:Index.End,2) + obj.Mob(Index.Start:Index.End, 2) .* obj.dxdz(Index.Start:Index.End,(c-1)*2+2, j) .* obj.rho(Index.Start:Index.End,2) + obj.Mob(Index.Start:Index.End, 2) .* obj.x(Index.Start:Index.End, (c-1)*2+2) .* obj.drhodz(Index.Start:Index.End, 2, j));

                        FluxDerivativePh1 = dMupPh1 .* abs(obj.U{1,1+f});
                        FluxDerivativePh2 = dMupPh2 .* abs(obj.U{2,1+f});

                        if c == obj.NofComponents
                            acc = - pv/dt .* rhoT + pv/dt .* obj.z(Index.Start:Index.End,c) .* obj.drhoTdz(Index.Start:Index.End, j);
                        elseif c == j
                            acc = pv/dt .* rhoT + pv/dt .* obj.z(Index.Start:Index.End,c) .* obj.drhoTdz(Index.Start:Index.End, j);
                        else
                            acc = 0 + pv/dt .* obj.z(Index.Start:Index.End,c) .* obj.drhoTdz(Index.Start:Index.End, j);
                        end

                        Jz_comp{j} = sparse(1:N, 1:N, acc)...
                            + (UpwindPermutation1Ph1 * sparse(1:N_Face, 1:N_Face, FluxDerivativePh1) * UpwindPermutation2Ph1')'...
                            + (UpwindPermutation1Ph2 * sparse(1:N_Face, 1:N_Face, FluxDerivativePh2) * UpwindPermutation2Ph2')';
                        % Capillarity
                        Jz_comp{j} = Jz_comp{j} - obj.Tph{2*c-1,1+f} * sparse(1:N, 1:N, obj.dPcdz(Index.Start:Index.End));
                    end
                otherwise
                    error("Discretization not found.")
            end
            % Assembling Jz
            Jz = Jz_comp{1};
            for j=2:obj.NofComponents-1
                Jz = horzcat(Jz, Jz_comp{j});
            end
        end
        %%
        % function Jacobian = BuildFullJacobian(obj, ProductionSystem, DiscretizationModel, Time, dt)
            % BUILD FIM JACOBIAN BLOCK BY BLOCK
            % Initialise local variables
            %             Nx = DiscretizationModel.ReservoirGrid.Nx;
            %             Ny = DiscretizationModel.ReservoirGrid.Ny;
            %             Nz = DiscretizationModel.ReservoirGrid.Nz;
            %             N = DiscretizationModel.ReservoirGrid.N;
            %             %
            %             pv = DiscretizationModel.ReservoirGrid.Volume*ProductionSystem.Reservoir.Por;
            %             z = zeros(N, obj.NofComponents);
            %             rho = zeros(N, obj.NofPhases);
            %             P = zeros(N, obj.NofPhases);
            %             x = zeros(N, obj.NofComponents*obj.NofPhases);
            % Copy values in local variables
            %             for i=1:obj.NofPhases
            %                 obj.P(:, i) = ProductionSystem.Reservoir.State.Properties(['P_', num2str(i)]).Value;
            %                 obj.rho(:, i) = ProductionSystem.Reservoir.State.Properties(['rho_', num2str(i)]).Value;
            %             end
            %             for i=1:obj.NofComponents
            %                 obj.z(:, i) = ProductionSystem.Reservoir.State.Properties(['z_', num2str(i)]).Value;
            %                 for j=1:obj.NofPhases
            %                     obj.x(:,(i-1)*obj.NofPhases + j) = ProductionSystem.Reservoir.State.Properties(['x_', num2str(i),'ph',num2str(j)]).Value;
            %                 end
            %             end
            %             rhoT =  ProductionSystem.Reservoir.State.Properties('rhoT').Value;
            % Fill in block by block
            %             Jp = cell(obj.NofComponents, 1);
            %             Jz = cell(obj.NofComponents, obj.NofComponents-1);
            %             for i=1:obj.NofComponents
            %                 %% 1. Component i pressure block
            %                 % 1.a: divergence
            %                 % Jp{i} = obj.Tph{i,1}  + obj.Tph{i, 2};
            %                 Jp{i} = obj.Tph{2*i-1,1}  + obj.Tph{2*i, 1};
            %
            %
            %                 % 1.b: compressibility part
            %                 dMupxPh1 = obj.UpWind{1}.x * (obj.Mob(:, 1) .* obj.x(:,(i-1)*2+1) .* obj.drhodp(:,1) + obj.Mob(:, 1) .* obj.dxdp(:,(i-1)*2+1) .* obj.rho(:,1) + obj.dMobdp(:, 1) .* obj.x(:,(i-1)*2+1) .* obj.rho(:,1));
            %                 dMupyPh1 = obj.UpWind{1}.y * (obj.Mob(:, 1) .* obj.x(:,(i-1)*2+1) .* obj.drhodp(:,1) + obj.Mob(:, 1) .* obj.dxdp(:,(i-1)*2+1) .* obj.rho(:,1) + obj.dMobdp(:, 1) .* obj.x(:,(i-1)*2+1) .* obj.rho(:,1));
            %                 dMupzPh1 = obj.UpWind{1}.z * (obj.Mob(:, 1) .* obj.x(:,(i-1)*2+1) .* obj.drhodp(:,1) + obj.Mob(:, 1) .* obj.dxdp(:,(i-1)*2+1) .* obj.rho(:,1) + obj.dMobdp(:, 1) .* obj.x(:,(i-1)*2+1) .* obj.rho(:,1));
            %                 dMupxPh2 = obj.UpWind{2}.x * (obj.Mob(:, 2) .* obj.x(:,(i-1)*2+2) .* obj.drhodp(:,2) + obj.Mob(:, 2) .* obj.dxdp(:,(i-1)*2+2) .* obj.rho(:,2) + obj.dMobdp(:, 2) .* obj.x(:,(i-1)*2+2) .* obj.rho(:,2));
            %                 dMupyPh2 = obj.UpWind{2}.y * (obj.Mob(:, 2) .* obj.x(:,(i-1)*2+2) .* obj.drhodp(:,2) + obj.Mob(:, 2) .* obj.dxdp(:,(i-1)*2+2) .* obj.rho(:,2) + obj.dMobdp(:, 2) .* obj.x(:,(i-1)*2+2) .* obj.rho(:,2));
            %                 dMupzPh2 = obj.UpWind{2}.z * (obj.Mob(:, 2) .* obj.x(:,(i-1)*2+2) .* obj.drhodp(:,2) + obj.Mob(:, 2) .* obj.dxdp(:,(i-1)*2+2) .* obj.rho(:,2) + obj.dMobdp(:, 2) .* obj.x(:,(i-1)*2+2) .* obj.rho(:,2));
            %
            %                 vecX1 = min(reshape(obj.U{1}.x(1:Nx,:,:),N,1), 0).*dMupxPh1 + min(reshape(obj.U{2}.x(1:Nx,:,:),N,1), 0).*dMupxPh2;
            %                 vecX2 = max(reshape(obj.U{1}.x(2:Nx+1,:,:),N,1), 0).*dMupxPh1 + max(reshape(obj.U{2}.x(2:Nx+1,:,:),N,1), 0).*dMupxPh2;
            %                 vecY1 = min(reshape(obj.U{1}.y(:,1:Ny,:),N,1), 0).*dMupyPh1 + min(reshape(obj.U{2}.y(:,1:Ny,:),N,1), 0).*dMupyPh2;
            %                 vecY2 = max(reshape(obj.U{1}.y(:,2:Ny+1,:),N,1), 0).*dMupyPh1 + max(reshape(obj.U{2}.y(:,2:Ny+1,:),N,1), 0).*dMupyPh2;
            %                 vecZ1 = min(reshape(obj.U{1}.z(:,:,1:Nz),N,1), 0).*dMupzPh1 + min(reshape(obj.U{2}.z(:,:,1:Nz),N,1), 0).*dMupzPh2;
            %                 vecZ2 = max(reshape(obj.U{1}.z(:,:,2:Nz+1),N,1), 0).*dMupzPh1 + max(reshape(obj.U{2}.z(:,:,2:Nz+1),N,1), 0).*dMupzPh2;
            %                 acc = pv/dt .* (obj.drhoTdp .* obj.z(:,i));
            %                 DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2 + vecY2+vecX2-vecY1-vecX1-vecZ1+acc, vecX1, vecY1, vecZ1];
            %                 DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
            %                 Jp{i} = Jp{i} + spdiags(DiagVecs, DiagIndx, N, N);
            %
            %
            %
            %                 %% 2. Component i composition block
            %                 for j=1:obj.NofComponents-1
            %
            %                     dMupxPh1 = obj.UpWind{1}.x*...
            %                         ( obj.dMob(:, 1, j) .* obj.x(:,(i-1)*2+1) .* obj.rho(:,1)...
            %                         + obj.Mob(:, 1) .* obj.dxdz(:,(i-1)*2+1, j) .* obj.rho(:,1)...
            %                         + obj.Mob(:, 1) .* obj.x(:, (i-1)*2+1) .* obj.drhodz(:, 1));
            %                     dMupyPh1 = obj.UpWind{1}.y*...
            %                         ( obj.dMob(:, 1, j) .* obj.x(:, (i-1)*2+1) .* obj.rho(:,1)...
            %                         + obj.Mob(:, 1) .* obj.dxdz(:, (i-1)*2+1, j) .* obj.rho(:,1)...
            %                         + obj.Mob(:, 1) .* obj.x(:,(i-1)*2+1) .* obj.drhodz(:, 1));
            %                     dMupzPh1 = obj.UpWind{1}.z*...
            %                         ( obj.dMob(:, 1, j) .* obj.x(:,(i-1)*2+1) .* obj.rho(:,1)...
            %                         + obj.Mob(:,1) .* obj.dxdz(:,(i-1)*2+1) .* obj.rho(:,1)...
            %                         + obj.Mob(:,1) .* obj.x(:,(i-1)*2+1) .* obj.drhodz(:, 1));
            %                     dMupxPh2 = obj.UpWind{2}.x*...
            %                         ( obj.dMob(:, 2, j) .* obj.x(:,(i-1)*2+2) .* obj.rho(:,2)...
            %                         + obj.Mob(:,2) .* obj.dxdz(:,(i-1)*2+2) .* obj.rho(:,2)...
            %                         + obj.Mob(:,2) .* obj.x(:,(i-1)*2+2) .* obj.drhodz(:, 2));
            %                     dMupyPh2 = obj.UpWind{2}.y*...
            %                         ( obj.dMob(:, 2, j) .* obj.x(:,(i-1)*2+2) .* obj.rho(:,2)...
            %                         + obj.Mob(:, 2) .* obj.dxdz(:,(i-1)*2+2, j) .* obj.rho(:,2)...
            %                         + obj.Mob(:, 2) .* obj.x(:, (i-1)*2+2) .* obj.drhodz(:, 2));
            %                     dMupzPh2 = obj.UpWind{2}.z*...
            %                         ( obj.dMob(:, 2, j) .* obj.x(:, (i-1)*2+2) .* obj.rho(:,2)...
            %                         + obj.Mob(:,2) .* obj.dxdz(:, (i-1)*2+2, j) .* obj.rho(:,2)...
            %                         + obj.Mob(:,2) .* obj.x(:, (i-1)*2+2) .* obj.drhodz(:, 2));
            %
            %                     vecX1 = min(reshape(obj.U{1}.x(1:Nx,:,:),N,1), 0).*dMupxPh1 + min(reshape(obj.U{2}.x(1:Nx,:,:),N,1), 0).*dMupxPh2;
            %                     vecX2 = max(reshape(obj.U{1}.x(2:Nx+1,:,:),N,1), 0).*dMupxPh1 + max(reshape(obj.U{2}.x(2:Nx+1,:,:),N,1), 0).*dMupxPh2;
            %                     vecY1 = min(reshape(obj.U{1}.y(:,1:Ny,:),N,1), 0).*dMupyPh1 + min(reshape(obj.U{2}.y(:,1:Ny,:),N,1), 0).*dMupyPh2;
            %                     vecY2 = max(reshape(obj.U{1}.y(:,2:Ny+1,:),N,1), 0).*dMupyPh1 + max(reshape(obj.U{2}.y(:,2:Ny+1,:),N,1), 0).*dMupyPh2;
            %                     vecZ1 = min(reshape(obj.U{1}.z(:,:,1:Nz),N,1), 0).*dMupzPh1 + min(reshape(obj.U{2}.z(:,:,1:Nz),N,1), 0).*dMupzPh2;
            %                     vecZ2 = max(reshape(obj.U{1}.z(:,:,2:Nz+1),N,1), 0).*dMupzPh1 + max(reshape(obj.U{2}.z(:,:,2:Nz+1),N,1), 0).*dMupzPh2;
            %
            %                     if i == obj.NofComponents
            %                         acc = - pv/dt .* rhoT + pv/dt .* obj.z(:,i) .* obj.drhoTdz(:, j);
            %                     elseif i == j
            %                         acc = pv/dt .* rhoT + pv/dt .* obj.z(:,i) .* obj.drhoTdz(:, j);
            %                     end
            %
            %                     DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2+vecY2+vecX2-vecZ1-vecY1-vecX1+acc, vecX1, vecY1, vecZ1];
            %                     DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
            %                     Jz{i,j} = spdiags(DiagVecs,DiagIndx, N, N);
            %                     % Capillarity
            %                     % Jz{i,j} = Jz{i,j} - obj.Tph{i,1} * spdiags(obj.dPcdz, 0, N, N);
            %                     Jz{i,j} = Jz{i,j} - obj.Tph{2*i-1,1} * spdiags(obj.dPcdz, 0, N, N);
            %
            %                 end
            %             end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Jacobian = BuildFullJacobian(obj, ProductionSystem, DiscretizationModel, Time, dt)
            Jhorizontal = cell(obj.NofComponents, 1);
            Nt = DiscretizationModel.N;
            Nm = DiscretizationModel.ReservoirGrid.N;

            for c = 1:obj.NofComponents
                Jph{c} = sparse(Nt, obj.NofComponents*Nt);
                % Reservoir
                Index.Start = 1;
                Index.End = Nm;
                [Jp_res, Jz_res] = BuildMediumJacobian(obj, ProductionSystem.Reservoir, DiscretizationModel.ReservoirGrid, dt, Index, 0, c);

                % Wells
                for i_well = 1:ProductionSystem.Wells.NofProd
                    Prod = ProductionSystem.Wells.Prod{i_well};
                    % Check if the well is open/close
                    isOpen = mod(sum(Time >= Prod.Schedule), 2);
                    if isOpen
                        [Jp_res, Jz_res] = obj.AddWellsToJacobian(Jp_res, Jz_res, obj.P(:,obj.NofPhases), obj.rho, obj.x, Prod, ProductionSystem.Reservoir.K, c);
                    end
                end

                % Fractures
                Jp = [];  Jz = [];
                for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                    Nf = DiscretizationModel.FracturesGrid.N;
                    Index.Start = DiscretizationModel.Index_Local_to_Global( Nm, f, 1     );
                    Index.End   = DiscretizationModel.Index_Local_to_Global( Nm, f, Nf(f) );
                    [Jp_frac, Jz_frac] = BuildMediumJacobian(obj, ProductionSystem.FracturesNetwork.Fractures(f), DiscretizationModel.FracturesGrid.Grids(f), dt, Index, f, c);
                    Jp  = blkdiag(Jp, Jp_frac);
                    Jz  = blkdiag(Jz, Jz_frac);
                end

                Jp  = blkdiag(Jp_res,Jp);
                Jz  = blkdiag(Jz_res,Jz);

                Depth = DiscretizationModel.ReservoirGrid.Depth(:);    % matrix part
                ofs = Nm;
                for fid = 1:ProductionSystem.FracturesNetwork.NumOfFrac
                    dloc  = DiscretizationModel.FracturesGrid.Grids(fid,1).Depth(:);
                    nfc   = numel(dloc);
                    Depth(ofs+1:ofs+nfc) = dloc;
                    ofs = ofs + nfc;
                end

                Depth = Depth .* 0;

                % frac-matrix and frac-frac connections
                for i=1:length(DiscretizationModel.CrossConnections) % The length is equal to the number of fracture cells
                    if isempty(DiscretizationModel.CrossConnections(i).Cells),  continue;  end
                    T_Geo = DiscretizationModel.CrossConnections(i).T_Geo;
                    UpWind = DiscretizationModel.CrossConnections(i).UpWind; % boolean vector (0 and 1)
                    idx_m = DiscretizationModel.CrossConnections(i).Cells; % idx_m and idx_f are global indexes; idx_m is a vector
                    idx_f = i + Nm;

                    depthDiff    = (Depth(idx_m) - Depth(idx_f));
                    rho_avg_1 = (obj.rho(idx_m,1) + obj.rho(idx_f,1)) ./ 2;
                    rho_avg_2 = (obj.rho(idx_m,2) + obj.rho(idx_f,2)) ./ 2;

                    DeltaP1 = (obj.P(idx_m,1) - obj.P(idx_f,1)) - rho_avg_1 .* obj.GravityModel.g .* depthDiff;
                    dDeltaP1_dPm = 1 - 0.5 * obj.drhodp(idx_m,1) .* obj.GravityModel.g .* depthDiff;
                    dDeltaP1_dPf = 1 - 0.5 * obj.drhodp(idx_f,1) .* obj.GravityModel.g .* depthDiff;
                    dDeltaP1_dzm = - 0.5 * obj.drhodz(idx_m,1) .* obj.GravityModel.g .* depthDiff;
                    dDeltaP1_dzf = - 0.5 * obj.drhodz(idx_f,1) .* obj.GravityModel.g .* depthDiff;
                
                    % Phase 2
                    DeltaP2 = (obj.P(idx_m,2) - obj.P(idx_f,2)) - rho_avg_2 .* obj.GravityModel.g .* depthDiff;
                    dDeltaP2_dPm = 1 - 0.5 * obj.drhodp(idx_m,2) .* obj.GravityModel.g .* depthDiff;
                    dDeltaP2_dPf = 1 - 0.5 * obj.drhodp(idx_f,2) .* obj.GravityModel.g .* depthDiff;
                    dDeltaP2_dzm = - 0.5 * obj.drhodz(idx_m,2) .* obj.GravityModel.g .* depthDiff;
                    dDeltaP2_dzf = - 0.5 * obj.drhodz(idx_f,2) .* obj.GravityModel.g .* depthDiff;

                    % 1. Pressure block
                    Jp_conn = - T_Geo .* ( ...
                        UpWind(:,1) .* ( ...
                        obj.Mob(idx_m,1)    .* obj.rho(idx_m,1).*obj.x(idx_m,(c-1)*2+1) .* dDeltaP1_dPm + ...
                        obj.Mob(idx_m,1)    .* obj.rho(idx_m,1).*obj.dxdp(idx_m,(c-1)*2+1).*DeltaP1 + ...
                        obj.Mob(idx_m,1)    .* obj.drhodp(idx_m,1).*obj.x(idx_m,(c-1)*2+1).*DeltaP1  + ...
                        obj.dMobdp(idx_m,1) .* obj.rho(idx_m,1).*obj.x(idx_m,(c-1)*2+1).*DeltaP1)  + ...
                    ~UpWind(:,1) .*  ( ...
                        obj.Mob(idx_f,1) .*obj.rho(idx_f,1).*obj.x(idx_f,(c-1)*2+1) .* dDeltaP1_dPf + ...
                        obj.Mob(idx_f,1) .*obj.rho(idx_f,1).*obj.dxdp(idx_f,(c-1)*2+1).*DeltaP1 + ...
                        obj.Mob(idx_f,1) .*obj.drhodp(idx_f,1).*obj.x(idx_f,(c-1)*2+1).*DeltaP1 + ...
                        obj.dMobdp(idx_f,1) .*obj.rho(idx_f,1).*obj.x(idx_f,(c-1)*2+1).*DeltaP1) + ...
                    UpWind(:,2) .* ( ...
                        obj.Mob(idx_m,2) .* obj.rho(idx_m,2).*obj.x(idx_m,(c-1)*2+2) .* dDeltaP2_dPm + ...
                        obj.Mob(idx_m,2) .* obj.rho(idx_m,2).*obj.dxdp(idx_m,(c-1)*2+2).*DeltaP2 + ...
                        obj.Mob(idx_m,2) .* obj.drhodp(idx_m,2).*obj.x(idx_m,(c-1)*2+2).*DeltaP2  + ...
                        obj.dMobdp(idx_m,2) .* obj.rho(idx_m,2).*obj.x(idx_m,(c-1)*2+2).*DeltaP2)  + ...
                    ~UpWind(:,2) .*  ( ...
                        obj.Mob(idx_f,2) .*obj.rho(idx_f,2).*obj.x(idx_f,(c-1)*2+2) .* dDeltaP2_dPf + ...
                        obj.Mob(idx_f,2) .*obj.rho(idx_f,2).*obj.dxdp(idx_f,(c-1)*2+2).*DeltaP2 + ...
                        obj.Mob(idx_f,2) .*obj.drhodp(idx_f,2).*obj.x(idx_f,(c-1)*2+2).*DeltaP2 + ...
                        obj.dMobdp(idx_f,2) .*obj.rho(idx_f,2).*obj.x(idx_f,(c-1)*2+2).*DeltaP2 ) );

                    Jz1_conn = T_Geo .* ( ...
                        UpWind(:,1) .* (obj.drhodz(idx_m,1) .* obj.Mob(idx_m,1) .* obj.x(idx_m,(c-1)*2+1) .* (-DeltaP1) + obj.rho(idx_m,1) .* obj.dMob(idx_m,1) .* obj.x(idx_m,(c-1)*2+1) .* (-DeltaP1) + obj.rho(idx_m,1) .* obj.Mob(idx_m,1) .* obj.dxdz(idx_m,(c-1)*2+1) .* (-DeltaP1) + obj.rho(idx_m,1) .* obj.Mob(idx_m,1) .* obj.x(idx_m,(c-1)*2+1) .* (-dDeltaP1_dzm))  + ...
                        UpWind(:,2) .* (obj.drhodz(idx_m,2) .* obj.Mob(idx_m,2) .* obj.x(idx_m,(c-1)*2+2) .* (-DeltaP2) + obj.rho(idx_m,2) .* obj.dMob(idx_m,2) .* obj.x(idx_m,(c-1)*2+2) .* (-DeltaP2) + obj.rho(idx_m,2) .* obj.Mob(idx_m,2) .* obj.dxdz(idx_m,(c-1)*2+2) .* (-DeltaP2) + obj.rho(idx_m,2) .* obj.Mob(idx_m,2) .* obj.x(idx_m,(c-1)*2+2) .* (-dDeltaP2_dzm)) );

                    Jz2_conn = T_Geo .* ( ...
                        ~UpWind(:,1) .* (obj.drhodz(idx_f,1) .* obj.Mob(idx_f,1).* obj.x(idx_f,(c-1)*2+1) .* (-DeltaP1) + obj.rho(idx_f,1) .* obj.dMob(idx_f,1).* obj.x(idx_f,(c-1)*2+1) .*(-DeltaP1) + obj.rho(idx_f,1) .* obj.Mob(idx_f,1) .* obj.dxdz(idx_f,(c-1)*2+1) .* (-DeltaP1) + obj.rho(idx_f,1) .* obj.Mob(idx_f,1) .* obj.x(idx_f,(c-1)*2+1) .* (-dDeltaP1_dzf)) + ...
                        ~UpWind(:,2) .* (obj.drhodz(idx_f,2) .* obj.Mob(idx_f,2).* obj.x(idx_f,(c-1)*2+2) .* (-DeltaP2) + obj.rho(idx_f,2) .* obj.dMob(idx_f,2).* obj.x(idx_f,(c-1)*2+2) .*(-DeltaP2) + obj.rho(idx_f,2) .* obj.Mob(idx_f,2) .* obj.dxdz(idx_f,(c-1)*2+2) .* (-DeltaP2) + obj.rho(idx_f,2) .* obj.Mob(idx_f,2) .* obj.x(idx_f,(c-1)*2+2) .* (-dDeltaP2_dzf)));

                    %%% original one without gravity
                    % 1. Pressure block
                    % Jp_conn = - T_Geo .* (...
                    %     UpWind(:,1) .* obj.Mob(idx_m, 1) .* (obj.rho(idx_m,1).* obj.x(idx_m,(c-1)*2+1) + obj.rho(idx_m,1) .* obj.dxdp(idx_m,(c-1)*2+1) .* ( obj.P(idx_m,1) - obj.P(idx_f,1) ) + obj.drhodp(idx_m,1).*obj.x(idx_m,(c-1)*2+1).*( obj.P(idx_m,1) - obj.P(idx_f,1) )) + ...
                    %     ~UpWind(:,1) .* obj.Mob(idx_f, 1) .* (obj.rho(idx_f,1).* obj.x(idx_f,(c-1)*2+1) + obj.rho(idx_f,1) .* obj.dxdp(idx_f,(c-1)*2+1) .* ( obj.P(idx_m,1) - obj.P(idx_f,1) ) + obj.drhodp(idx_f,1).*obj.x(idx_f,(c-1)*2+1).*( obj.P(idx_m,1) - obj.P(idx_f,1))) + ...
                    %     UpWind(:,2) .* obj.Mob(idx_m, 2) .* (obj.rho(idx_m,2).* obj.x(idx_m,(c-1)*2+2) + obj.rho(idx_m,2) .* obj.dxdp(idx_m,(c-1)*2+2) .* ( obj.P(idx_m,2) - obj.P(idx_f,2) ) + obj.drhodp(idx_m,2).*obj.x(idx_m,(c-1)*2+2).*( obj.P(idx_m,2) - obj.P(idx_f,2) )) + ...
                    %     ~UpWind(:,2) .* obj.Mob(idx_f, 2) .* (obj.rho(idx_f,2).* obj.x(idx_f,(c-1)*2+2) + obj.rho(idx_f,2) .* obj.dxdp(idx_f,(c-1)*2+2) .* ( obj.P(idx_m,2) - obj.P(idx_f,2) ) + obj.drhodp(idx_f,2).*obj.x(idx_f,(c-1)*2+2).*( obj.P(idx_m,2) - obj.P(idx_f,2))) );
                     % 2. Composition block
                    % % Jz1_conn: dqf(idx_f)/dz(idx_m) for upwind cell being matrix
                    % Jz1_conn = T_Geo .* (...
                    %     UpWind(:,1) .* obj.rho(idx_m,1) .* ( obj.dMob(idx_m,1) .* obj.x(idx_m,(c-1)*2+1) + obj.Mob(idx_m,1) .* obj.dxdz(idx_m,(c-1)*2+1) ) .* (obj.P(idx_f,1) - obj.P(idx_m,1)) + ...
                    %     UpWind(:,2) .* obj.rho(idx_m,2) .* ( obj.dMob(idx_m,2) .* obj.x(idx_m,(c-1)*2+2) + obj.Mob(idx_m,2) .* obj.dxdz(idx_m,(c-1)*2+2) ) .* (obj.P(idx_f,2) - obj.P(idx_m,2)) );
                    % % Jz2_conn: dqf(idx_f)/dz(idx_f) for upwind cell being frac
                    % Jz2_conn = T_Geo .* (...
                    %     ~UpWind(:,1) .* obj.rho(idx_f,1) .* ( obj.dMob(idx_f,1) .* obj.x(idx_f,(c-1)*2+1) + obj.Mob(idx_f,1) .* obj.dxdz(idx_f,(c-1)*2+1) ) .* (obj.P(idx_f,1) - obj.P(idx_m,1)) + ...
                    %     ~UpWind(:,2) .* obj.rho(idx_f,2) .* ( obj.dMob(idx_f,2) .* obj.x(idx_f,(c-1)*2+2) + obj.Mob(idx_f,2) .* obj.dxdz(idx_f,(c-1)*2+2) ) .* (obj.P(idx_f,2) - obj.P(idx_m,2)) );
                    % frac - mat or frac1 - frac2
                    Jp(idx_f, idx_m) = Jp_conn;
                    Jp(idx_f, idx_f) = Jp(idx_f, idx_f) - sum(Jp_conn);
                    % mat-frac or frac2 - frac1
                    Jp(idx_m, idx_f) = Jp_conn';
                    Jp(sub2ind([Nt, Nt], idx_m, idx_m)) = Jp(sub2ind([Nt, Nt], idx_m, idx_m)) - Jp_conn;

                   

                    % frac - mat or frac1 - frac2
                    Jz(idx_f, idx_m) = Jz1_conn;
                    Jz(idx_f, idx_f) = Jz(idx_f, idx_f) + sum(Jz2_conn);
                    % mat-frac or frac2 - frac1
                    Jz(idx_m, idx_f) = -Jz2_conn;
                    % diag of mat or frac2
                    Jz(sub2ind([Nt, Nt], idx_m, idx_m)) = Jz(sub2ind([Nt, Nt], idx_m, idx_m)) - Jz1_conn;
                end

                % Combine the Jacobian blocks
                Jhorizontal{c} = horzcat(Jp, Jz);
            end

           % % 3. Add wells to each block
           %  if (Time <= obj.Injection_time) % well shutdown control
           %      [Jp, Jz] = obj.AddWellsToJacobian(Jp, Jz, obj.P(:,obj.NofPhases), obj.rho, obj.x, ProductionSystem.Wells, ProductionSystem.Reservoir.K);
           %  end
           % 
           % % Full Jacobian
           %  if obj.NofComponents <= 2
           %  Jacobian = [Jp{1}, Jz{1,1}
           %              Jp{2}, Jz{2,1}];
           %  else
           %  % for now let's consider a max of 3 components!!
           %  Jacobian = [Jp{1}, Jz{1,1}, Jz{1,2};...
           %              Jp{2}, Jz{2,1}, Jz{2,2};...
           %              Jp{3}, Jz{3,1}, Jz{3,2}];
           %  end

            Jacobian = Jhorizontal{1};
            %Jacobian_noconn = Jphnoconn{1};
            for c=2:obj.NofComponents
                Jacobian = vertcat(Jacobian, Jhorizontal{c});
            end

        end
        %%
        function [Jp, Jz] = AddWellsToJacobian(obj, Jp, Jz, p, rho, x, Prod, K, c)
            N = length(p);
            b = Prod.Cells;
            for j=1:length(b)
                %Pressure blocks
                Jp(b(j),b(j)) = Jp(b(j),b(j))...
                    + Prod.PI(j) * K(b(j)) *...
                    ( obj.Mob(b(j), 1) * rho(b(j),1) * x(b(j), (c-1)*2+1) ...
                    + obj.Mob(b(j), 2) * rho(b(j),2) * x(b(j), (c-1)*2+2)...
                    )...
                    - Prod.PI(j) * K(b(j)) * obj.Mob(b(j), 1) * x(b(j), (c-1)*2+1) * obj.drhodp(b(j),1) * (Prod.p(j) - p(b(j))) ...
                    - Prod.PI(j) * K(b(j)) * obj.Mob(b(j), 2) * x(b(j), (c-1)*2+2) * obj.drhodp(b(j),2) * (Prod.p(j) - p(b(j))) ...
                    - Prod.PI(j) * K(b(j)) * obj.Mob(b(j), 1) * obj.dxdp(b(j), (c-1)*2+1) * rho(b(j),1) * (Prod.p(j) - p(b(j))) ...
                    - Prod.PI(j) * K(b(j)) * obj.Mob(b(j), 2) * obj.dxdp(b(j), (c-1)*2+2) * rho(b(j),2) * (Prod.p(j) - p(b(j))) ...
                    - Prod.PI(j) * K(b(j)) * obj.dMobdp(b(j), 1) * x(b(j), (c-1)*2+1) * rho(b(j),1) * (Prod.p(j) - p(b(j))) ...
                    - Prod.PI(j) * K(b(j)) * obj.dMobdp(b(j), 2) * x(b(j), (c-1)*2+2) * rho(b(j),2) * (Prod.p(j) - p(b(j)));
                for zI = 1:obj.NofComponents-1
                    Jz(b(j),(zI-1)*N+b(j)) = Jz(b(j),(zI-1)*N+b(j))...
                        - Prod.PI(j) * K(b(j)) *...
                        ( obj.dMob(b(j), 1, zI) * rho(b(j), 1) * x(b(j), (c-1)*2+1) ...
                        + obj.dMob(b(j), 2, zI) * rho(b(j), 2) * x(b(j), (c-1)*2+2)...
                        + obj.Mob(b(j), 1) * rho(b(j), 1) * obj.dxdz(b(j), (c-1)*2+1, zI) ...
                        + obj.Mob(b(j), 2) * rho(b(j), 2) * obj.dxdz(b(j), (c-1)*2+2, zI)...
                        + obj.Mob(b(j), 1) * obj.drhodz(b(j), 1, zI) * x(b(j), (c-1)*2+1) ...
                        + obj.Mob(b(j), 2) * obj.drhodz(b(j), 2, zI) * x(b(j), (c-1)*2+2)...
                        ) * (Prod.p(j) - p(b(j)));
                end
            end
        end
        %%
        function delta = UpdateState(obj, delta, ProductionSystem, FluidModel, DiscretizationModel)
            if sum(isnan(delta))
                % if the solution makes no sense, skip this step
                return
            else
                Nm = DiscretizationModel.ReservoirGrid.N;
                Nt = DiscretizationModel.N;
                deltaP = delta(1:Nt);
                deltaz = delta(Nt+1: end);
                %% 1. Update matrix
                % 1.a Update Pressure
                Pm = ProductionSystem.Reservoir.State.Properties(['P_', num2str(obj.NofPhases)]);
                Pm.update(deltaP(1:Nm));
                % 1.b Update z
                DeltaLast = zeros(Nm, 1);
                for c = 1:obj.NofComponents-1
                    Zm = ProductionSystem.Reservoir.State.Properties(['z_', num2str(c)]);
                    Deltaz = delta(c*Nt + 1:c*Nt + Nm);
                    % Deltaz = deltaz(1:Nm);

                    Zm.update(Deltaz);
                    Zm.Value = max(Zm.Value, 1.0e-8);
                    Zm.Value = min(Zm.Value, 1-1.0e-8);
                    DeltaLast = DeltaLast + Deltaz;
                end
                Zm = ProductionSystem.Reservoir.State.Properties(['z_', num2str(obj.NofComponents)]);
                Zm.update(-DeltaLast);
                Zm.Value = max(Zm.Value, 1.0e-8);
                Zm.Value = min(Zm.Value, 1-1.0e-8);
                % Update Remaining properties (x, S, rho, rhoT, Pc)
                obj.UpdatePhaseCompositions(ProductionSystem.Reservoir, DiscretizationModel, FluidModel, 1);

                %% 2. Update fractures pressure and densities
                if ProductionSystem.FracturesNetwork.Active
                    EP = Nm;
                    Nf = DiscretizationModel.FracturesGrid.N;
                    for f = 1:ProductionSystem.FracturesNetwork.NumOfFrac
                        IP = EP + 1;
                        EP = IP + Nf(f) - 1;
                        % 2.a Update Pressure
                        Pf = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['P_', num2str(obj.NofPhases)]);
                        Pf.update(deltaP(IP:EP));
                        % 2.b Update z
                        DeltaLast = zeros(Nf(f), 1);
                        for c = 1:obj.NofComponents-1
                            IZ = Nt*(c-1) + IP;
                            EZ = Nt*(c-1) + EP;
                            Zf = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['z_', num2str(c)]);
                            % Zf.update(delta(c*Nf(f) + 1:(c+1)*Nf(f)));
                            Deltaz = deltaz(IZ:EZ);

                            Zf.update(Deltaz);
                            % Remove values that are not physical
                            Zf.Value = max(Zf.Value, 1.0e-8);
                            Zf.Value = min(Zf.Value, 1-1.0e-8);
                            DeltaLast = DeltaLast + Deltaz;
                        end
                        Zf = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['z_', num2str(obj.NofComponents)]);
                        Zf.update(-DeltaLast);
                        % Remove values that are not physical
                        Zf.Value = max(Zf.Value, 1.0e-8);
                        Zf.Value = min(Zf.Value, 1-1.0e-8);
                        % Update Remaining properties (x, S, rho, rhoT, Pc)
                        obj.UpdatePhaseCompositions(ProductionSystem.FracturesNetwork.Fractures(f), DiscretizationModel, FluidModel, 0);
                    end
                end
            end
        end
        %%
        function qf = ComputeSourceTerms_frac_mat_wog(obj, ProductionSystem, DiscretizationModel)
            %             qf = zeros(DiscretizationModel.N, obj.NofPhases);
            %             Nm = DiscretizationModel.ReservoirGrid.N;
            %             % Global variables
            %             % P = ProductionSystem.CreateGlobalVariables([DiscretizationModel.ReservoirGrid; DiscretizationModel.FracturesGrid.Grids], obj.NofPhases, 'P_'); % useful for cross connections assembly
            %             % rho = ProductionSystem.CreateGlobalVariables([DiscretizationModel.ReservoirGrid; DiscretizationModel.FracturesGrid.Grids], obj.NofPhases, 'rho_'); % useful for cross connections assembly
            %             for ph=1:obj.NofPhases
            %                 % fill in qf
            %                 for c=1:length(DiscretizationModel.CrossConnections) % The length is equal to the number of fracture cells
            %                     j = DiscretizationModel.CrossConnections(c).Cells; % j and i are global indexes
            %                     i = c + Nm;
            %                     T_Geo = DiscretizationModel.CrossConnections(c).T_Geo;
            %                     UpWind = DiscretizationModel.CrossConnections(c).UpWind; % boolean vector (0 and 1)
            %                     qf(i, ph) = qf(i, ph) + sum(T_Geo .* ...
            %                         ( UpWind(:, ph) .* obj.Mob(j, ph) .* obj.rho(j, ph) .* (obj.P(j, ph) - obj.P(i, ph)) + ...
            %                         ~UpWind(:, ph) .* obj.Mob(i, ph) .* obj.rho(i, ph) .* (obj.P(j, ph) - obj.P(i, ph))) );
            %                     qf(j, ph) = qf(j, ph) + T_Geo .* ...
            %                         ( UpWind(:, ph) .* obj.Mob(j, ph) .* obj.rho(j, ph) .* (obj.P(i, ph) - obj.P(j, ph)) + ...
            %                         ~UpWind(:, ph) .* obj.Mob(i, ph) .* obj.rho(i, ph) .* (obj.P(i, ph) - obj.P(j, ph)) );
            %                 end
            %             end

            qf = zeros(DiscretizationModel.N, obj.NofComponents);
            Nm = DiscretizationModel.ReservoirGrid.N;

            for c=1:obj.NofComponents
                % fill in qf
                for i=1:length(DiscretizationModel.CrossConnections) % The length is equal to the number of fracture cells
                    T_Geo = DiscretizationModel.CrossConnections(i).T_Geo;
                    UpWind = DiscretizationModel.CrossConnections(i).UpWind; % boolean vector (0 and 1)
                    idx_m = DiscretizationModel.CrossConnections(i).Cells; % idx_m and idx_f are global indexes; idx_m is a vector
                    idx_f = i + Nm;
                    % qf = T_Geo * Upwind * Mob * rho * x * (dP)
                    qf(idx_f, c) = qf(idx_f, c) + sum(T_Geo .* ...
                        ( UpWind(:,1) .* obj.Mob(idx_m,1) .* obj.rho(idx_m,1) .*obj.x(idx_m,(c-1)*2+1) .* (obj.P(idx_m,1) - obj.P(idx_f,1)) + ...
                        ~UpWind(:,1) .* obj.Mob(idx_f,1) .* obj.rho(idx_f,1) .*obj.x(idx_f,(c-1)*2+1) .* (obj.P(idx_m,1) - obj.P(idx_f,1)))+ ...
                        T_Geo .* ...
                        ( UpWind(:,2) .* obj.Mob(idx_m,2) .* obj.rho(idx_m,2) .*obj.x(idx_m,(c-1)*2+2) .* (obj.P(idx_m,2) - obj.P(idx_f,2)) + ...
                        ~UpWind(:,2) .* obj.Mob(idx_f,2) .* obj.rho(idx_f,2) .*obj.x(idx_f,(c-1)*2+2) .* (obj.P(idx_m,2) - obj.P(idx_f,2))) );

                    qf(idx_m, c) = qf(idx_m, c) + T_Geo .* ...
                        ( UpWind(:,1) .* obj.Mob(idx_m,1) .* obj.rho(idx_m,1) .*obj.x(idx_m,(c-1)*2+1) .* (obj.P(idx_f,1) - obj.P(idx_m,1)) + ...
                        ~UpWind(:,1) .* obj.Mob(idx_f,1) .* obj.rho(idx_f,1) .*obj.x(idx_f,(c-1)*2+1) .* (obj.P(idx_f,1) - obj.P(idx_m,1)))+ ...
                        T_Geo .*...
                        ( UpWind(:,2) .* obj.Mob(idx_m,2) .* obj.rho(idx_m,2) .*obj.x(idx_m,(c-1)*2+2) .* (obj.P(idx_f,2) - obj.P(idx_m,2)) + ...
                        ~UpWind(:,2) .* obj.Mob(idx_f,2) .* obj.rho(idx_f,2) .*obj.x(idx_f,(c-1)*2+2) .* (obj.P(idx_f,2) - obj.P(idx_m,2)));
                end
            end
        end

        function qf = ComputeSourceTerms_frac_mat(obj, ProductionSystem, DiscretizationModel)
            qf = zeros(DiscretizationModel.N, obj.NofComponents);
            Nm = DiscretizationModel.ReservoirGrid.N;
            Depth = DiscretizationModel.ReservoirGrid.Depth(:);
            offset = Nm;
            for frac_id = 1:DiscretizationModel.FracturesGrid.Nfrac
                frac_depth = DiscretizationModel.FracturesGrid.Grids(frac_id, 1).Depth(:);
                num_frac_cells = length(frac_depth);
                Depth(offset+1 : offset+num_frac_cells) = frac_depth;
                offset = offset + num_frac_cells;
            end

            Depth = Depth .* 0;

            for c = 1:obj.NofComponents
                for i = 1:length(DiscretizationModel.CrossConnections)
                    T_Geo = DiscretizationModel.CrossConnections(i).T_Geo;
                    UpWind = DiscretizationModel.CrossConnections(i).UpWind;
                    idx_m = DiscretizationModel.CrossConnections(i).Cells;
                    idx_f = i + Nm;

                    depthDiff    = (Depth(idx_m) - Depth(idx_f));

                    rho_avg_1 = (obj.rho(idx_m,1) + obj.rho(idx_f,1)) ./ 2;
                    rho_avg_2 = (obj.rho(idx_m,2) + obj.rho(idx_f,2)) ./ 2;

                    qf(idx_f, c) = qf(idx_f, c) + sum( T_Geo .* ...
                        ( UpWind(:,1) .* obj.Mob(idx_m,1) .* obj.rho(idx_m,1) .* obj.x(idx_m,(c-1)*2+1) .* ...
                        ( (obj.P(idx_m,1) - obj.P(idx_f,1)) - rho_avg_1 .* obj.GravityModel.g .* depthDiff ) + ...
                        ~UpWind(:,1) .* obj.Mob(idx_f,1) .* obj.rho(idx_f,1) .* obj.x(idx_f,(c-1)*2+1) .* ...
                        ( (obj.P(idx_m,1) - obj.P(idx_f,1)) - rho_avg_1 .* obj.GravityModel.g .* depthDiff ) ) + ...
                        T_Geo .* ...
                        ( UpWind(:,2) .* obj.Mob(idx_m,2) .* obj.rho(idx_m,2) .* obj.x(idx_m,(c-1)*2+2) .* ...
                        ( (obj.P(idx_m,2) - obj.P(idx_f,2)) - rho_avg_2 .* obj.GravityModel.g .* depthDiff ) + ...
                        ~UpWind(:,2) .* obj.Mob(idx_f,2) .* obj.rho(idx_f,2) .* obj.x(idx_f,(c-1)*2+2) .* ...
                        ( (obj.P(idx_m,2) - obj.P(idx_f,2)) - rho_avg_2 .* obj.GravityModel.g .* depthDiff ) ) );
                    qf(idx_m, c) = qf(idx_m, c) + T_Geo .* ...
                        ( UpWind(:,1) .* obj.Mob(idx_m,1) .* obj.rho(idx_m,1) .* obj.x(idx_m,(c-1)*2+1) .* ...
                        ( (obj.P(idx_f,1) - obj.P(idx_m,1)) - rho_avg_1 .* obj.GravityModel.g .* (-depthDiff) ) + ...
                        ~UpWind(:,1) .* obj.Mob(idx_f,1) .* obj.rho(idx_f,1) .* obj.x(idx_f,(c-1)*2+1) .* ...
                        ( (obj.P(idx_f,1) - obj.P(idx_m,1)) - rho_avg_1 .* obj.GravityModel.g .* (-depthDiff) ) ) + ...
                        T_Geo .* ...
                        ( UpWind(:,2) .* obj.Mob(idx_m,2) .* obj.rho(idx_m,2) .* obj.x(idx_m,(c-1)*2+2) .* ...
                        ( (obj.P(idx_f,2) - obj.P(idx_m,2)) - rho_avg_2 .* obj.GravityModel.g .* (-depthDiff) ) + ...
                        ~UpWind(:,2) .* obj.Mob(idx_f,2) .* obj.rho(idx_f,2) .* obj.x(idx_f,(c-1)*2+2) .* ...
                        ( (obj.P(idx_f,2) - obj.P(idx_m,2)) - rho_avg_2 .* obj.GravityModel.g .* (-depthDiff) ) );
                end
            end
        end
        %%
        %         function [Jp, Jz] = AddWellsToJacobian(obj, Jp, Jz, p, rho, x, Wells, K)
        %             Inj = Wells.Inj;
        %             Prod = Wells.Prod;
        %             % Injectors
        %             % for pressure conrol; (commented for rate control)
        % %             for i=1:Wells.NofInj
        % %                 a = Inj{i}.Cells;
        % %                 for j=1:length(a)
        % %                     for c=1:obj.NofComponents
        % %                     Jp{c}(a(j),a(j)) = Jp{c}(a(j),a(j)) ...
        % %                                        + Inj{i}.PI(j) * K(a(j)) * (Inj{i}.Mob(1,1) * Inj{i}.rho(j,1) * Inj{i}.x(1,(c-1)*2+1) ...
        % %                                        + Inj{i}.Mob(1,2) *  Inj{i}.rho(j,2) * Inj{i}.x(1,(c-1)*2+2));
        % %                     end
        % %                 end
        % %             end
        %
        %             % Producers
        %             for i=1:Wells.NofProd
        %                 b = Prod{i}.Cells;
        %                 for j=1:length(b)
        %                     for c=1:obj.NofComponents
        %                         %Pressure blocks
        %                         Jp{c}(b(j),b(j)) = Jp{c}(b(j),b(j))...
        %                             + Prod{i}.PI(j) * K(b(j)) *...
        %                             ( obj.Mob(b(j), 1) * rho(b(j),1) * x(b(j), (c-1)*2+1) ...
        %                             + obj.Mob(b(j), 2) * rho(b(j),2) * x(b(j), (c-1)*2+2)...
        %                             )...
        %                             - Prod{i}.PI(j) * K(b(j)) * obj.Mob(b(j), 1) * x(b(j), (c-1)*2+1) * obj.drhodp(b(j),1) * (Prod{i}.p(j) - p(b(j))) ...
        %                             - Prod{i}.PI(j) * K(b(j)) * obj.Mob(b(j), 1) * obj.dxdp(b(j), (c-1)*2+1) * rho(b(j),1) * (Prod{i}.p(j) - p(b(j))) ...
        %                             - Prod{i}.PI(j) * K(b(j)) * obj.Mob(b(j), 2) * x(b(j), (c-1)*2+2) * obj.drhodp(b(j),2) * (Prod{i}.p(j) - p(b(j)))...
        %                             - Prod{i}.PI(j) * K(b(j)) * obj.Mob(b(j), 2) * obj.dxdp(b(j), (c-1)*2+2) * rho(b(j),2) * (Prod{i}.p(j) - p(b(j)));
        %                         for zI = 1:obj.NofComponents-1
        %                             Jz{c, zI}(b(j),b(j)) = Jz{c, zI}(b(j),b(j))...
        %                                 - Prod{i}.PI(j) * K(b(j)) *...
        %                                 ( obj.dMob(b(j), 1, zI) * rho(b(j), 1) * x(b(j), (c-1)*2+1) ...
        %                                 + obj.dMob(b(j), 2, zI) * rho(b(j), 2) * x(b(j), (c-1)*2+2)...
        %                                 + obj.Mob(b(j), 1) * rho(b(j), 1) * obj.dxdz(b(j), (c-1)*2+1, zI) ...
        %                                 + obj.Mob(b(j), 2) * rho(b(j), 2) * obj.dxdz(b(j), (c-1)*2+2, zI)...
        %                                 + obj.Mob(b(j), 1, zI) * obj.drhodz(b(j), 1) * x(b(j), (c-1)*2+1) ...
        %                                 + obj.Mob(b(j), 2, zI) * obj.drhodz(b(j), 2) * x(b(j), (c-1)*2+2)...
        %                                 ) * (Prod{i}.p(j) - p(b(j)));
        %                         end
        %                     end
        %                 end
        %             end
        %         end

    end
end