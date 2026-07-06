%Composition Energy Formulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Mengjie Zhao
%TU Delft
%Created: 29 Jan 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Overall_Composition_Energy_formulation < Compositional_formulation
    properties
        dxdz
        dxdp
        dxdT

        drhoTdp
        drhoTdz
        drhoTdT

        drhodz
        drhodT

        dSdp
        dSdz
        dSdT

        dMobdp
        dMobdz
        dMobdT

        dhdp
        dhdz
        dhdT
        
        rho_old
        rho
        P
        z_old
        z
        x

        h  % fluid enthalpy
        Tk % transmisibility of rock conductivity
        T  % temperature of rock, gas, liquid

        Thph
        Ghph
    end
    methods
        function obj = Overall_Composition_Energy_formulation(SimulationInput,n_phases,n_components)
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

            Start = End + 1;
            End = Start + Nm - 1;
            x(Start:End) = ProductionSystem.Reservoir.State.Properties('T').Value;

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

                    Start = End + 1;
                    End = Start + Nf(f) - 1;
                    x(Start:End) = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('T').Value;
                end
            end
        end
        %%
        function ComputePropertiesAndDerivatives(obj, ProductionSystem, FluidModel, ~)
            %% 1. reservoir
            % This is the bitchy part! (Indeed)
            [obj.dxdp, obj.dxdz, obj.dxdT] = FluidModel.DxDpDzDT(ProductionSystem.Reservoir.State, ProductionSystem.Reservoir.State.Properties('SinglePhase').Value);
            % %%%%%%%%%%%%%%%%%%%%%%%%
            obj.drhodp = FluidModel.ComputeDrhoDp(ProductionSystem.Reservoir.State, obj.dxdp);
            obj.drhodz = FluidModel.ComputeDrhoDz(ProductionSystem.Reservoir.State, obj.dxdz);
            obj.drhodT = FluidModel.ComputeDrhoDT(ProductionSystem.Reservoir.State, obj.dxdT);

            obj.dSdp = FluidModel.ComputeDSDp(ProductionSystem.Reservoir.State, obj.drhodp, obj.dxdp(:, end));
            obj.dSdz = FluidModel.ComputeDSDz(ProductionSystem.Reservoir.State, obj.drhodz, obj.dxdz(:, end, :));
            obj.dSdT = FluidModel.ComputeDSDT(ProductionSystem.Reservoir.State, obj.drhodT, obj.dxdT(:, end));
            
            obj.dhdp = FluidModel.ComputeDhDp(ProductionSystem.Reservoir.State, obj.dxdp);
            obj.dhdz = FluidModel.ComputeDhDz(ProductionSystem.Reservoir.State, obj.dxdz);
            obj.dhdT = FluidModel.ComputeDhDT(ProductionSystem.Reservoir.State, obj.dxdT);

            obj.drhoTdp = FluidModel.ComputeDrhotDp(ProductionSystem.Reservoir.State,obj.drhodp, obj.dSdp);
            obj.drhoTdz = FluidModel.ComputeDrhotDz(ProductionSystem.Reservoir.State,obj.drhodz, obj.dSdz);
            obj.drhoTdT = FluidModel.ComputeDrhotDT(ProductionSystem.Reservoir.State,obj.drhodT, obj.dSdT);

            if FluidModel.HysteresisActive
                obj.dMobdp = FluidModel.ComputeDMobDpHysteresis(ProductionSystem.Reservoir.State, [obj.dSdp, obj.dSdp], obj.dxdp);
                obj.dMobdz = FluidModel.ComputeDMobDzHysteresis(ProductionSystem.Reservoir.State, [obj.dSdz, obj.dSdz], obj.dxdz);
                obj.dMobdT = FluidModel.ComputeDMobDTHysteresis(ProductionSystem.Reservoir.State, [obj.dSdT, obj.dSdT], obj.dxdT);
            else
                obj.dMobdp = FluidModel.ComputeDMobDp(ProductionSystem.Reservoir.State, [obj.dSdp, obj.dSdp], obj.dxdp); % I am not using it !
                obj.dMobdz = FluidModel.ComputeDMobDz(ProductionSystem.Reservoir.State, [obj.dSdz, obj.dSdz], obj.dxdz);
                obj.dMobdT = FluidModel.ComputeDMobDT(ProductionSystem.Reservoir.State, [obj.dSdT, obj.dSdT], obj.dxdT);
            end
            obj.dPc = FluidModel.ComputeDPcDz(ProductionSystem.Reservoir.State, obj.dSdz, 1);
            %% 2. fracture
            for f = 1:ProductionSystem.FracturesNetwork.NumOfFrac

                [dxdp_temp, dxdz_temp, dxdT_temp] = FluidModel.DxDpDzDT(ProductionSystem.FracturesNetwork.Fractures(f).State, ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('SinglePhase').Value);

                % %%%%%%%%%%%%%%%%%%%%%%%%
                % Mengjie:
                % obj.dxdp > dxdp_temp   & obj.dxdz > dxdz_temp           
                drhodp_temp = FluidModel.ComputeDrhoDp(ProductionSystem.FracturesNetwork.Fractures(f).State, dxdp_temp);
                drhodz_temp = FluidModel.ComputeDrhoDz(ProductionSystem.FracturesNetwork.Fractures(f).State, dxdz_temp);
                drhodT_temp = FluidModel.ComputeDrhoDT(ProductionSystem.FracturesNetwork.Fractures(f).State, dxdT_temp);

                dSdp_temp = FluidModel.ComputeDSDp(ProductionSystem.FracturesNetwork.Fractures(f).State, drhodp_temp, dxdp_temp(:,end));
                dSdz_temp = FluidModel.ComputeDSDz(ProductionSystem.FracturesNetwork.Fractures(f).State, drhodz_temp, dxdz_temp(:,end, :));
                dSdT_temp = FluidModel.ComputeDSDT(ProductionSystem.FracturesNetwork.Fractures(f).State, drhodT_temp, dxdT_temp(:,end));
                
                obj.dhdp = [obj.dhdp; FluidModel.ComputeDhDp(ProductionSystem.FracturesNetwork.Fractures(f).State, dxdp_temp)];
                obj.dhdz = [obj.dhdz; FluidModel.ComputeDhDz(ProductionSystem.FracturesNetwork.Fractures(f).State, dxdz_temp)];
                obj.dhdT = [obj.dhdT; FluidModel.ComputeDhDT(ProductionSystem.FracturesNetwork.Fractures(f).State, dxdT_temp)];

                obj.drhoTdp = [obj.drhoTdp; FluidModel.ComputeDrhotDp(ProductionSystem.FracturesNetwork.Fractures(f).State, drhodp_temp, dSdp_temp)];
                obj.drhoTdz = [obj.drhoTdz; FluidModel.ComputeDrhotDz(ProductionSystem.FracturesNetwork.Fractures(f).State, drhodz_temp, dSdz_temp)];
                obj.drhoTdT = [obj.drhoTdT; FluidModel.ComputeDrhotDT(ProductionSystem.FracturesNetwork.Fractures(f).State, drhodT_temp, dSdT_temp)];
       
                obj.dMobdp  = [obj.dMobdp;  FluidModel.ComputeDMobDp(ProductionSystem.FracturesNetwork.Fractures(f).State, [dSdp_temp, dSdp_temp], dxdp_temp)];
                obj.dMobdz  = [obj.dMobdz;  FluidModel.ComputeDMobDz(ProductionSystem.FracturesNetwork.Fractures(f).State, [dSdz_temp, dSdz_temp], dxdz_temp)];
                obj.dMobdT  = [obj.dMobdT;  FluidModel.ComputeDMobDT(ProductionSystem.FracturesNetwork.Fractures(f).State, [dSdT_temp, dSdT_temp], dxdT_temp)];
                obj.dPc     = [obj.dPc;     FluidModel.ComputeDPcDz(ProductionSystem.FracturesNetwork.Fractures(f).State, dSdz_temp, 0)];

                obj.dxdp = [obj.dxdp; dxdp_temp];
                obj.dxdz = [obj.dxdz; dxdz_temp];
                obj.dxdT = [obj.dxdT; dxdT_temp];
                obj.drhodp = [obj.drhodp; drhodp_temp];
                obj.drhodz = [obj.drhodz; drhodz_temp];
                obj.drhodT = [obj.drhodT; drhodT_temp];
                obj.dSdp = [obj.dSdp; dSdp_temp];
                obj.dSdz = [obj.dSdz; dSdz_temp];
                obj.dSdT = [obj.dSdT; dSdT_temp];
            end

        end

        %%
        function [Residual, RHS] = BuildMassMediumResidual(obj, Medium, Grid, dt, State0, Index, qw, qf, f, c)
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

        function [Residual, RHS] = BuildEnergyMediumResidual(obj, Medium, Grid, dt, State0, Index, qw_h, qf_h, f)
            I = Index.Start:Index.End;
            phi = Medium.Por(I);
            N = Grid.N;
            V = Grid.Volume;

            A = sparse(1:N, 1:N, V/dt);

            n_1 = Medium.State.Properties('S_1').Value(I) .* ( Medium.State.Properties('rho_1').Value(I) .* Medium.State.Properties('h_1').Value(I) - Medium.State.Properties('P_1').Value(I));
            n_old_1 = State0.Properties('S_1').Value(I) .* ( State0.Properties('rho_1').Value(I) .* State0.Properties('h_1').Value(I) - State0.Properties('P_1').Value(I));

            n_2 = Medium.State.Properties('S_2').Value(I) .* (Medium.State.Properties('rho_2').Value(I) .* Medium.State.Properties('h_2').Value(I) - Medium.State.Properties('P_2').Value(I));
            n_old_2 = State0.Properties('S_2').Value(I) .* (State0.Properties('rho_2').Value(I) .* State0.Properties('h_2').Value(I) - State0.Properties('P_2').Value(I) );
            
            RHS = sum(qw_h(I,:), 2) + sum(qf_h(I,:), 2);

            Residual = A * ( (1 - phi) .* Medium.Rho .* Medium.State.Properties('hRock').Value(I) + phi .* ( n_1 + n_2 ) - ( ...
                (1 - phi) .* Medium.Rho .* State0.Properties('hRock').Value(I) + phi .* ( n_old_1 + n_old_2 ) ) ) ...
                + ( obj.Thph{1,1+f} * Medium.State.Properties('P_1').Value(I) - obj.Ghph{1,1+f} * Grid.Depth ) ...
                + ( obj.Thph{2,1+f} * Medium.State.Properties('P_2').Value(I) - obj.Ghph{2,1+f} * Grid.Depth ) ...
                + obj.Tk{1,1+f} * Medium.State.Properties('T').Value(I) ...
                - sum(qw_h(I,:), 2) ...
                - sum(qf_h(I,:), 2); 
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
                obj.h(Index.Start:Index.End, i) = ProductionSystem.Reservoir.State.Properties(['h_', num2str(i)]).Value;
                obj.S(Index.Start:Index.End, i) = ProductionSystem.Reservoir.State.Properties(['S_', num2str(i)]).Value;  
            end
            obj.T(Index.Start:Index.End, 1) = ProductionSystem.Reservoir.State.Properties('T').Value;

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
                    obj.h(Index.Start:Index.End, i) = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['h_', num2str(i)]).Value;
                    obj.S(Index.Start:Index.End, i) = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['S_', num2str(i)]).Value;
                end
                obj.T(Index.Start:Index.End, 1) = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('T').Value;
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
            Residual = zeros(Nt * (obj.NofComponents + 1), 1);
            RHS = zeros(Nt * (obj.NofComponents + 1), 1);
            Nm = DiscretizationModel.ReservoirGrid.N;

            % Copy values to local variables
            CopytoLocalVariables(obj, ProductionSystem, DiscretizationModel, State0);

            % Source terms
            % well
            qw = obj.ComputeSourceTerms(Nt, ProductionSystem.Wells);
            qw_h = obj.ComputeSourceTerms_heat(Nt, ProductionSystem.Wells);

            % fracture_matrix
            qf = zeros(Nt, obj.NofComponents);
            qf_h = zeros(Nt, obj.NofPhases);
            if ProductionSystem.FracturesNetwork.Active
                qf = ComputeSourceTerms_frac_mat(obj, ProductionSystem, DiscretizationModel);
                qf_h = ComputeSourceTerms_frac_mat_heat(obj, ProductionSystem, DiscretizationModel);
            end

            % unique(qf)
            % Build Mass Balance Residual
            for c = 1:obj.NofComponents
                % Reservoir
                % Index
                Index.Start = 1;
                Index.End = Nm;
                % Phase Transmissibilities
                obj.TransmissibilityMatrix(DiscretizationModel.ReservoirGrid, ProductionSystem.Reservoir, obj.rho(Index.Start:Index.End,:), obj.x(Index.Start:Index.End,(c-1)*2+1:(c-1)*2+2), Index, c, 0);
                % Residual
                [Residualm, RHSm] = BuildMassMediumResidual(obj, ProductionSystem.Reservoir, DiscretizationModel.ReservoirGrid, dt, State0, Index, qw, qf, 0, c);
                Residual((c-1)*Nt + Index.Start: (c-1)*Nt + Index.End) = Residualm;
                RHS(     (c-1)*Nt + Index.Start: (c-1)*Nt + Index.End) = RHSm;

                % Fractures
                for f = 1:ProductionSystem.FracturesNetwork.NumOfFrac
                    % Index
                    Index.Start = Index.End + 1;
                    Index.End = Index.Start + DiscretizationModel.FracturesGrid.N(f) - 1;
                    % Phase Transmissibilities
                    obj.TransmissibilityMatrix(DiscretizationModel.FracturesGrid.Grids(f), ProductionSystem.FracturesNetwork.Fractures(f), obj.rho(Index.Start:Index.End,:), obj.x(Index.Start:Index.End,(c-1)*2+1:(c-1)*2+2), Index, c, f);
                    % Residual
                    Residual_frac_f = BuildMassMediumResidual(obj, ProductionSystem.FracturesNetwork.Fractures(f), DiscretizationModel.FracturesGrid.Grids(f), dt, State0, Index, qw, qf, f, c);
                    Residual((c-1)*Nt + Index.Start: (c-1)*Nt + Index.End) = Residual_frac_f;
                end
            end

            off = obj.NofComponents * Nt;

            % Build Energy Balance Residual Reservoir
            Index.Start = 1;
            Index.End = Nm;
            for ph = 1 : obj.NofPhases
                % heat transmissibility
                obj.ConvectiveHeatTransmissibilityMatrix(DiscretizationModel.ReservoirGrid, Index, obj.rho(Index.Start:Index.End,:), obj.h(Index.Start:Index.End,:), ph, 0);
            end

            obj.ConductiveHeatTransmissibilityMatrix( DiscretizationModel.ReservoirGrid , ProductionSystem.Reservoir.State, 0);
            [Residual_EB, RHS_EB] = BuildEnergyMediumResidual(obj, ProductionSystem.Reservoir, DiscretizationModel.ReservoirGrid, dt, State0, Index, qw_h, qf_h, 0);
            Residual(off + Index.Start : off + Index.End ) = Residual_EB;
            RHS     (off + Index.Start : off + Index.End ) = RHS_EB;

            % Fractures
            for f = 1:ProductionSystem.FracturesNetwork.NumOfFrac
                Index.Start = Index.End + 1;
                Index.End = Index.Start + DiscretizationModel.FracturesGrid.N(f) - 1;
                for ph = 1 : obj.NofPhases
                    obj.ConvectiveHeatTransmissibilityMatrix(DiscretizationModel.FracturesGrid.Grids(f), Index, obj.rho(Index.Start:Index.End,:), obj.h(Index.Start:Index.End,:), ph, f);
                end

                obj.ConductiveHeatTransmissibilityMatrix(DiscretizationModel.FracturesGrid.Grids(f), ProductionSystem.Fracture(f).State, f)
                Residual_EB_f = BuildEnergyMediumResidual(obj, ProductionSystem.FracturesNetwork.Fractures(f), DiscretizationModel.FracturesGrid.Grids(f), dt, State0, Index, qw_h, qf_h, 0);
                Residual(off + Index.Start : off + Index.End ) = Residual_EB_f;
            end
        end

        %% corner_point_grid case is added to BuildMediumJacobian function
        function [Jp, Jz, JT] = BuildMassMediumJacobian(obj, Medium, Grid, dt, Index, f, c)
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
            Jz = sparse(N, N * (obj.NofComponents-1));
            JT = sparse(N, N);
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
                    % 2. Component c composition block
                    for j=1:obj.NofComponents-1
                        dMupPh1 = obj.UpWind{1,1+f} * ( obj.dMobdz(Index.Start:Index.End, 1, j) .* obj.x(Index.Start:Index.End,(c-1)*2+1) .* obj.rho(Index.Start:Index.End,1) + obj.Mob(Index.Start:Index.End, 1) .* obj.dxdz(Index.Start:Index.End,(c-1)*2+1, j) .* obj.rho(Index.Start:Index.End,1) + obj.Mob(Index.Start:Index.End, 1) .* obj.x(Index.Start:Index.End, (c-1)*2+1) .* obj.drhodz(Index.Start:Index.End, 1));
                        dMupPh2 = obj.UpWind{2,1+f} * ( obj.dMobdz(Index.Start:Index.End, 2, j) .* obj.x(Index.Start:Index.End,(c-1)*2+2) .* obj.rho(Index.Start:Index.End,2) + obj.Mob(Index.Start:Index.End, 2) .* obj.dxdz(Index.Start:Index.End,(c-1)*2+2, j) .* obj.rho(Index.Start:Index.End,2) + obj.Mob(Index.Start:Index.End, 2) .* obj.x(Index.Start:Index.End, (c-1)*2+2) .* obj.drhodz(Index.Start:Index.End, 2));

                        FluxDerivativePh1 = dMupPh1 .* abs(obj.U{1,1+f});
                        FluxDerivativePh2 = dMupPh2 .* abs(obj.U{2,1+f});

                        if c == obj.NofComponents
                            acc = - pv/dt .* rhoT + pv/dt .* obj.z(Index.Start:Index.End,c) .* obj.drhoTdz(Index.Start:Index.End, j);
                        elseif c == j
                            acc = pv/dt .* rhoT + pv/dt .* obj.z(Index.Start:Index.End,c) .* obj.drhoTdz(Index.Start:Index.End, j);
                        end

                        Jz(:, 1+N*(j-1):N*j) = sparse(1:N, 1:N, acc)...
                            + (UpwindPermutation1Ph1 * sparse(1:N_Face, 1:N_Face, FluxDerivativePh1) * UpwindPermutation2Ph1')'...
                            + (UpwindPermutation1Ph2 * sparse(1:N_Face, 1:N_Face, FluxDerivativePh2) * UpwindPermutation2Ph2')';

                        %                         Jz = spdiags(1:N, 1:N, acc)...
                        %                             + (UpwindPermutation1Ph1 * spdiags(FluxDerivativePh1, 0, N_Face, N_Face) * UpwindPermutation2Ph1')'...
                        %                             + (UpwindPermutation1Ph2 * spdiags(FluxDerivativePh2, 0, N_Face, N_Face) * UpwindPermutation2Ph2')';

                        % Capillarity
                        Jz(:, 1+N*(j-1):N*j) = Jz(:, 1+N*(j-1):N*j) - obj.Tph{2*c-1,1+f} * sparse(1:N, 1:N, obj.dPc(Index.Start:Index.End));
                        % Jz = Jz - obj.Tph{2*c-1,1+f} * spdiags(obj.dPc(Index.Start:Index.End), 0, N, N);
                    end

                    % 3. temperature part
                    dMupPh1 = obj.UpWind{1,1+f} * (obj.Mob(Index.Start:Index.End, 1) .* obj.x(Index.Start:Index.End,(c-1)*2+1) .* obj.drhodT(Index.Start:Index.End,1) + obj.Mob(Index.Start:Index.End, 1) .* obj.dxdT(Index.Start:Index.End,(c-1)*2+1) .* obj.rho(Index.Start:Index.End,1) + obj.dMobdT(Index.Start:Index.End, 1) .* obj.x(Index.Start:Index.End,(c-1)*2+1) .* obj.rho(Index.Start:Index.End,1));
                    dMupPh2 = obj.UpWind{2,1+f} * (obj.Mob(Index.Start:Index.End, 2) .* obj.x(Index.Start:Index.End,(c-1)*2+2) .* obj.drhodT(Index.Start:Index.End,2) + obj.Mob(Index.Start:Index.End, 2) .* obj.dxdT(Index.Start:Index.End,(c-1)*2+2) .* obj.rho(Index.Start:Index.End,2) + obj.dMobdT(Index.Start:Index.End, 2) .* obj.x(Index.Start:Index.End,(c-1)*2+2) .* obj.rho(Index.Start:Index.End,2));

                    FluxDerivativePh1 = dMupPh1 .* abs(obj.U{1,1+f}); % nf x 1
                    FluxDerivativePh2 = dMupPh2 .* abs(obj.U{2,1+f});

                    acc = pv/dt .* (obj.drhoTdT(Index.Start:Index.End) .* obj.z(Index.Start:Index.End,c));
                    JT = JT + sparse(1:N, 1:N, acc)...
                        + (UpwindPermutation1Ph1 * sparse(1:N_Face, 1:N_Face, FluxDerivativePh1) * UpwindPermutation2Ph1')'... % ncxnf * nfxnf * nfxnc = ncxnc
                        + (UpwindPermutation1Ph2 * sparse(1:N_Face, 1:N_Face, FluxDerivativePh2) * UpwindPermutation2Ph2')';

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
                    NeighborPermutation1 = sparse(repmat((1:nf)', 1, 2), [CellNeighbor1Index, CellNeighbor2Index], [ones(nf,1), zeros(nf,1)], nf, nc)'; % nc x nf
                    NeighborPermutation2 = sparse(repmat((1:nf)', 1, 2), [CellNeighbor1Index, CellNeighbor2Index], [zeros(nf,1), ones(nf,1)], nf, nc)'; % nc x nf

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
                    
                    vecL = obj.U_D{2*c-1,1+f} .* dDT_ph1_L .* dS_rho_D_ph1(CellNeighbor1Index) + obj.U_D{2*c,  1+f} .* dDT_ph2_L .* dS_rho_D_ph2(Neighbor1Index);
                    vecR = obj.U_D{2*c-1,1+f} .* dDT_ph1_R .* dS_rho_D_ph1(CellNeighbor2Index) + obj.U_D{2*c,  1+f} .* dDT_ph2_R .* dS_rho_D_ph2(Neighbor2Index);

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
                        dMupPh1 = obj.UpWind{1,1+f} * ( obj.dMobdz(Index.Start:Index.End, 1, j) .* obj.x(Index.Start:Index.End,(c-1)*2+1) .* obj.rho(Index.Start:Index.End,1) + obj.Mob(Index.Start:Index.End, 1) .* obj.dxdz(Index.Start:Index.End,(c-1)*2+1, j) .* obj.rho(Index.Start:Index.End,1) + obj.Mob(Index.Start:Index.End, 1) .* obj.x(Index.Start:Index.End, (c-1)*2+1) .* obj.drhodz(Index.Start:Index.End, 1));
                        dMupPh2 = obj.UpWind{2,1+f} * ( obj.dMobdz(Index.Start:Index.End, 2, j) .* obj.x(Index.Start:Index.End,(c-1)*2+2) .* obj.rho(Index.Start:Index.End,2) + obj.Mob(Index.Start:Index.End, 2) .* obj.dxdz(Index.Start:Index.End,(c-1)*2+2, j) .* obj.rho(Index.Start:Index.End,2) + obj.Mob(Index.Start:Index.End, 2) .* obj.x(Index.Start:Index.End, (c-1)*2+2) .* obj.drhodz(Index.Start:Index.End, 2));

                        FluxDerivativePh1 = dMupPh1 .* abs(obj.U{1,1+f});
                        FluxDerivativePh2 = dMupPh2 .* abs(obj.U{2,1+f});

                        if c == obj.NofComponents
                            acc = - pv/dt .* rhoT + pv/dt .* obj.z(Index.Start:Index.End,c) .* obj.drhoTdz(Index.Start:Index.End, j);
                        elseif c == j
                            acc = pv/dt .* rhoT + pv/dt .* obj.z(Index.Start:Index.End,c) .* obj.drhoTdz(Index.Start:Index.End, j);
                        end

                        Jz(:, 1+N*(j-1):N*j) = sparse(1:N, 1:N, acc)...
                            + (UpwindPermutation1Ph1 * sparse(1:N_Face, 1:N_Face, FluxDerivativePh1) * UpwindPermutation2Ph1')'...
                            + (UpwindPermutation1Ph2 * sparse(1:N_Face, 1:N_Face, FluxDerivativePh2) * UpwindPermutation2Ph2')';

                        % Diffusion
                        dS_rho_D_ph1 =  obj.dSdz(Index.Start:Index.End,1,j) .* Medium.Diff(:,1) .* Grid.Por;
                        dS_rho_D_ph2 = -obj.dSdz(Index.Start:Index.End,1,j) .* Medium.Diff(:,2) .* Grid.Por;
                        vecL = obj.U_D{2*c-1,1+f} .* dDT_ph1_L .* dS_rho_D_ph1(CellNeighbor1Index) + obj.U_D{2*c,  1+f} .* dDT_ph2_L .* dS_rho_D_ph2(CellNeighbor1Index);
                        vecR = obj.U_D{2*c-1,1+f} .* dDT_ph1_R .* dS_rho_D_ph1(CellNeighbor2Index) + obj.U_D{2*c,  1+f} .* dDT_ph2_R .* dS_rho_D_ph2(CellNeighbor2Index);

                        Jz(:,1+N*(j-1):N*j) = Jz(:,1+N*(j-1):N*j) + (NeighborPermutation1 * sparse(1:nf, 1:nf, vecL, nf, nf) * Grid.ConnectivityMatrix')' ...
                            + (NeighborPermutation2 * sparse(1:nf, 1:nf, vecR, nf, nf) * Grid.ConnectivityMatrix')' ...
                            + obj.Dph{2*c-1,1+f} * spdiags(obj.dxdz(Index.Start:Index.End,2*c-1,j), 0, N, N) ...
                            + obj.Dph{2*c,  1+f} * spdiags(obj.dxdz(Index.Start:Index.End,2*c,  j), 0, N, N);

                        drhodz1 = obj.drhodz(Index.Start:Index.End,1,j);
                        drhodz2 = obj.drhodz(Index.Start:Index.End,2,j);
                        vecL = obj.U_D{2*c-1,1+f} .* (drhodz1(CellNeighbor1Index) / 2) .* Grid.Diff(:,1) + obj.U_D{2*c,  1+f} .* (drhodz2(CellNeighbor1Index) / 2) .* Grid.Diff(:,2);
                        vecR = obj.U_D{2*c-1,1+f} .* (drhodz1(CellNeighbor2Index) / 2) .* Grid.Diff(:,1) + obj.U_D{2*c,  1+f} .* (drhodz2(CellNeighbor2Index) / 2) .* Grid.Diff(:,2);

                        Jz(:,col) = Jz(:,col) + (NeighborPermutation1 * sparse(1:nf, 1:nf, vecL, nf, nf) * Grid.ConnectivityMatrix')' ...
                            + (NeighborPermutation2 * sparse(1:nf, 1:nf, vecR, nf, nf) * Grid.ConnectivityMatrix')';

                        % Capillarity
                        Jz(:, 1+N*(j-1):N*j) = Jz(:, 1+N*(j-1):N*j) - obj.Tph{2*c-1,1+f} * sparse(1:N, 1:N, obj.dPc(Index.Start:Index.End, j));
                    end

                    % 3. temperature part
                    dMupPh1 = obj.UpWind{1,1+f} * (obj.Mob(Index.Start:Index.End, 1) .* obj.x(Index.Start:Index.End,(c-1)*2+1) .* obj.drhodT(Index.Start:Index.End,1) + obj.Mob(Index.Start:Index.End, 1) .* obj.dxdT(Index.Start:Index.End,(c-1)*2+1) .* obj.rho(Index.Start:Index.End,1) + obj.dMobdT(Index.Start:Index.End, 1) .* obj.x(Index.Start:Index.End,(c-1)*2+1) .* obj.rho(Index.Start:Index.End,1));
                    dMupPh2 = obj.UpWind{2,1+f} * (obj.Mob(Index.Start:Index.End, 2) .* obj.x(Index.Start:Index.End,(c-1)*2+2) .* obj.drhodT(Index.Start:Index.End,2) + obj.Mob(Index.Start:Index.End, 2) .* obj.dxdT(Index.Start:Index.End,(c-1)*2+2) .* obj.rho(Index.Start:Index.End,2) + obj.dMobdT(Index.Start:Index.End, 2) .* obj.x(Index.Start:Index.End,(c-1)*2+2) .* obj.rho(Index.Start:Index.End,2));

                    FluxDerivativePh1 = dMupPh1 .* abs(obj.U{1,1+f}); % nf x 1
                    FluxDerivativePh2 = dMupPh2 .* abs(obj.U{2,1+f});

                    acc = pv/dt .* (obj.drhoTdT(Index.Start:Index.End) .* obj.z(Index.Start:Index.End,c));
                    JT = JT + sparse(1:N, 1:N, acc)...
                        + (UpwindPermutation1Ph1 * sparse(1:N_Face, 1:N_Face, FluxDerivativePh1) * UpwindPermutation2Ph1')'... % ncxnf * nfxnf * nfxnc = ncxnc
                        + (UpwindPermutation1Ph2 * sparse(1:N_Face, 1:N_Face, FluxDerivativePh2) * UpwindPermutation2Ph2')';

                    dS_rho_D_ph1 =  obj.dSdT(Index.Start:Index.End,1) .* Medium.Diff(:,1) .* Grid.Por;
                    dS_rho_D_ph2 = -obj.dSdT(Index.Start:Index.End,1) .* Medium.Diff(:,2) .* Grid.Por;

                    vecL = obj.U_D{2*c-1,1+f} .* dDT_ph1_L .* dS_rho_D_ph1(CellNeighbor1Index) + obj.U_D{2*c,  1+f} .* dDT_ph2_L .* dS_rho_D_ph2(CellNeighbor1Index);
                    vecR = obj.U_D{2*c-1,1+f} .* dDT_ph1_R .* dS_rho_D_ph1(CellNeighbor2Index) + obj.U_D{2*c,  1+f} .* dDT_ph2_R .* dS_rho_D_ph2(CellNeighbor2Index);

                    JT = JT + (NeighborPermutation1 * sparse(1:nf, 1:nf, vecL, nf, nf) * Grid.ConnectivityMatrix')' ...
                        + (NeighborPermutation2 * sparse(1:nf, 1:nf, vecR, nf, nf) * Grid.ConnectivityMatrix')' ...
                        + obj.Dph{2*c-1,1+f} * spdiags(obj.dxdT(Index.Start:Index.End,2*c-1), 0, nc, nc) ...
                        + obj.Dph{2*c,  1+f} * spdiags(obj.dxdT(Index.Start:Index.End,2*c  ), 0, nc, nc);

                    drhodT1 = obj.drhodT(Index.Start:Index.End,1);
                    drhodT2 = obj.drhodT(Index.Start:Index.End,2);
                    vecL = obj.U_D{2*c-1,1+f} .* (drhodT1(CellNeighbor1Index) / 2) .* Grid.Diff(:,1) + obj.U_D{2*c,  1+f} .* (drhodT2(CellNeighbor1Index) / 2) .* Grid.Diff(:,2);
                    vecR = obj.U_D{2*c-1,1+f} .* (drhodT1(CellNeighbor2Index) / 2) .* Grid.Diff(:,1) + obj.U_D{2*c,  1+f} .* (drhodT2(CellNeighbor2Index) / 2) .* Grid.Diff(:,2);

                    JT = JT + (NeighborPermutation1 * sparse(1:nf, 1:nf, vecL, nf, nf) * Grid.ConnectivityMatrix')' ...
                        + (NeighborPermutation2 * sparse(1:nf, 1:nf, vecR, nf, nf) * Grid.ConnectivityMatrix')';

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
                    Jp = Jp + spdiags(DiagVecs, DiagIndx, Grid.N, Grid.N) + obj.Dph{2*c-1,1+f} * spdiags(obj.dxdp(:,2*c-1), 0, N, N) + obj.Dph{2*c,1+f} * spdiags(obj.dxdp(:,2*c), 0, N, N);
                    
                    vecX1 = reshape(obj.U_D{2*c-1,1+f}.x(1:Nx,:,:),N,1).* obj.drhodp(Index.Start:Index.End,1) / 2 .* reshape(Grid.Diffx(1:Nx,:,:,1),N,1) + reshape(obj.U_D{2*c,1+f}.x(1:Nx,:,:),N,1) .* obj.drhodp(Index.Start:Index.End,2) / 2  .* reshape(Grid.Diffx(1:Nx,:,:,2),N,1);
                    vecX2 = reshape(obj.U_D{2*c-1,1+f}.x(2:Nx+1,:,:),N,1).* obj.drhodp(Index.Start:Index.End,1) / 2 .* reshape(Grid.Diffx(2:Nx+1,:,:,1),N,1) + reshape(obj.U_D{2*c,1+f}.x(2:Nx+1,:,:),N,1) .* obj.drhodp(Index.Start:Index.End,2) / 2  .* reshape(Grid.Diffx(2:Nx+1,:,:,2),N,1);
                    vecY1 = reshape(obj.U_D{2*c-1,1+f}.y(:,1:Ny,:),N,1).* obj.drhodp(Index.Start:Index.End,1) / 2.* reshape(Grid.Diffy(:,1:Ny,:,1),N,1) + reshape(obj.U_D{2*c,1+f}.y(:,1:Ny,:),N,1) .* obj.drhodp(Index.Start:Index.End,2) / 2  .* reshape(Grid.Diffy(:,1:Ny,:,2),N,1);
                    vecY2 = reshape(obj.U_D{2*c-1,1+f}.y(:,2:Ny+1,:),N,1).* obj.drhodp(Index.Start:Index.End,1) / 2 .* reshape(Grid.Diffy(:,2:Ny+1,:,1),N,1) + reshape(obj.U_D{2*c,1+f}.y(:,2:Ny+1,:),N,1) .* obj.drhodp(Index.Start:Index.End,2) / 2  .* reshape(Grid.Diffy(:,2:Ny+1,:,2),N,1);
                    vecZ1 = reshape(obj.U_D{2*c-1,1+f}.z(:,:,1:Nz),N,1) .* obj.drhodp(Index.Start:Index.End,1) / 2 .* reshape(Grid.Diffz(:,:,1:Nz,1),N,1) + reshape(obj.U_D{2*c,1+f}.z(:,:,1:Nz),N,1) .* obj.drhodp(Index.Start:Index.End,2) / 2  .* reshape(Grid.Diffz(:,:,1:Nz,2),N,1);
                    vecZ2 = reshape(obj.U_D{2*c-1,1+f}.z(:,:,2:Nz+1),N,1) .* obj.drhodp(Index.Start:Index.End,1) / 2 .* reshape(Grid.Diffz(:,:,2:Nz+1,1),N,1) + reshape(obj.U_D{2*c,1+f}.z(:,:,2:Nz+1),N,1) .* obj.drhodp(Index.Start:Index.End,2) / 2  .* reshape(Grid.Diffz(:,:,2:Nz+1,2),N,1);
                    
                    DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2 + vecY2 + vecX2 - vecY1 - vecX1 - vecZ1, vecX1, vecY1, vecZ1];
                    DiagIndx = [-Grid.Nx*Grid.Ny, -Grid.Nx, -1, 0, 1, Grid.Nx, Grid.Nx*Grid.Ny];
                    Jp = Jp + spdiags(DiagVecs, DiagIndx, Grid.N, Grid.N);

                    % 2. Component c composition block
                    for j=1:obj.NofComponents-1
                        dMupxPh1 = reshape(obj.UpWind{1,1+f}.x*...
                            ( obj.dMobdz(Index.Start:Index.End, 1, j) .* obj.x(Index.Start:Index.End,(c-1)*2+1) .* obj.rho(Index.Start:Index.End,1)...
                            + obj.Mob(Index.Start:Index.End, 1) .* obj.dxdz(Index.Start:Index.End,(c-1)*2+1, j) .* obj.rho(Index.Start:Index.End,1)...
                            + obj.Mob(Index.Start:Index.End, 1) .* obj.x(Index.Start:Index.End, (c-1)*2+1) .* obj.drhodz(Index.Start:Index.End, 1, j)), Nx, Ny, Nz);
                        dMupyPh1 = reshape(obj.UpWind{1,1+f}.y*...
                            ( obj.dMobdz(Index.Start:Index.End, 1, j) .* obj.x(Index.Start:Index.End, (c-1)*2+1) .* obj.rho(Index.Start:Index.End,1)...
                            + obj.Mob(Index.Start:Index.End, 1) .* obj.dxdz(Index.Start:Index.End, (c-1)*2+1, j) .* obj.rho(Index.Start:Index.End,1)...
                            + obj.Mob(Index.Start:Index.End, 1) .* obj.x(Index.Start:Index.End,(c-1)*2+1) .* obj.drhodz(Index.Start:Index.End, 1, j)), Nx, Ny, Nz);
                        dMupzPh1 = reshape(obj.UpWind{1,1+f}.z*...
                            ( obj.dMobdz(Index.Start:Index.End, 1, j) .* obj.x(Index.Start:Index.End,(c-1)*2+1) .* obj.rho(Index.Start:Index.End,1)...
                            + obj.Mob(Index.Start:Index.End,1) .* obj.dxdz(Index.Start:Index.End,(c-1)*2+1, j) .* obj.rho(Index.Start:Index.End,1)...
                            + obj.Mob(Index.Start:Index.End,1) .* obj.x(Index.Start:Index.End,(c-1)*2+1) .* obj.drhodz(Index.Start:Index.End, 1, j)), Nx, Ny, Nz);
                        dMupxPh2 = reshape(obj.UpWind{2,1+f}.x*...
                            ( obj.dMobdz(Index.Start:Index.End, 2, j) .* obj.x(Index.Start:Index.End,(c-1)*2+2) .* obj.rho(Index.Start:Index.End,2)...
                            + obj.Mob(Index.Start:Index.End,2) .* obj.dxdz(Index.Start:Index.End,(c-1)*2+2, j) .* obj.rho(Index.Start:Index.End,2)...
                            + obj.Mob(Index.Start:Index.End,2) .* obj.x(Index.Start:Index.End,(c-1)*2+2) .* obj.drhodz(Index.Start:Index.End, 2, j)), Nx, Ny, Nz);
                        dMupyPh2 = reshape(obj.UpWind{2,1+f}.y*...
                            ( obj.dMobdz(Index.Start:Index.End, 2, j) .* obj.x(Index.Start:Index.End,(c-1)*2+2) .* obj.rho(Index.Start:Index.End,2)...
                            + obj.Mob(Index.Start:Index.End, 2) .* obj.dxdz(Index.Start:Index.End,(c-1)*2+2, j) .* obj.rho(Index.Start:Index.End,2)...
                            + obj.Mob(Index.Start:Index.End, 2) .* obj.x(Index.Start:Index.End, (c-1)*2+2) .* obj.drhodz(Index.Start:Index.End, 2, j)), Nx, Ny, Nz);
                        dMupzPh2 = reshape(obj.UpWind{2,1+f}.z*...
                            ( obj.dMobdz(Index.Start:Index.End, 2, j) .* obj.x(Index.Start:Index.End, (c-1)*2+2) .* obj.rho(Index.Start:Index.End,2)...
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
                        Jz(:, 1+N*(j-1):N*j) = spdiags(DiagVecs,DiagIndx, N, N);

                        dS_rho_D_ph1_x = (obj.dSdz(Index.Start:Index.End,1, j) .* Medium.Diff(:,1) .* Grid.Por) .* Grid.Ax / Grid.dx;
                        dS_rho_D_ph1_y = (obj.dSdz(Index.Start:Index.End,1, j) .* Medium.Diff(:,1) .* Grid.Por) .* Grid.Ay / Grid.dy;
                        dS_rho_D_ph1_z = (obj.dSdz(Index.Start:Index.End,1, j) .* Medium.Diff(:,1) .* Grid.Por) .* Grid.Az / Grid.dz;
                        dS_rho_D_ph2_x = (-obj.dSdz(Index.Start:Index.End,1, j) .* Medium.Diff(:,2) .* Grid.Por) .* Grid.Ax / Grid.dx;
                        dS_rho_D_ph2_y = (-obj.dSdz(Index.Start:Index.End,1, j) .* Medium.Diff(:,2) .* Grid.Por) .* Grid.Ay / Grid.dy;
                        dS_rho_D_ph2_z = (-obj.dSdz(Index.Start:Index.End,1, j) .* Medium.Diff(:,2) .* Grid.Por) .* Grid.Az / Grid.dz;

                        vecX1 = reshape(obj.U_D{2*c-1,1+f}.x(1:Nx,:,:),N,1).* dS_rho_D_ph1_x .* reshape(dDTx_ph1_L,N,1) + reshape(obj.U_D{2*c,1+f}.x(1:Nx,:,:),N,1) .* dS_rho_D_ph2_x .* reshape(dDTx_ph2_L,N,1);
                        vecX2 = reshape(obj.U_D{2*c-1,1+f}.x(2:Nx+1,:,:),N,1).* dS_rho_D_ph1_x .* reshape(dDTx_ph1_R,N,1) + reshape(obj.U_D{2*c,1+f}.x(2:Nx+1,:,:),N,1) .* dS_rho_D_ph2_x .* reshape(dDTx_ph2_R,N,1);
                        vecY1 = reshape(obj.U_D{2*c-1,1+f}.y(:,1:Ny,:),N,1).* dS_rho_D_ph1_y .* reshape(dDTy_ph1_L,N,1) + reshape(obj.U_D{2*c,1+f}.y(:,1:Ny,:),N,1) .* dS_rho_D_ph2_y .* reshape(dDTy_ph2_L,N,1);
                        vecY2 = reshape(obj.U_D{2*c-1,1+f}.y(:,2:Ny+1,:),N,1).* dS_rho_D_ph1_y .* reshape(dDTy_ph1_R,N,1) + reshape(obj.U_D{2*c,1+f}.y(:,2:Ny+1,:),N,1) .* dS_rho_D_ph2_y .* reshape(dDTy_ph2_R,N,1);
                        vecZ1 = reshape(obj.U_D{2*c-1,1+f}.z(:,:,1:Nz),N,1) .* dS_rho_D_ph1_z .* reshape(dDTz_ph1_L,N,1) + reshape(obj.U_D{2*c,1+f}.z(:,:,1:Nz),N,1) .* dS_rho_D_ph2_z .* reshape(dDTz_ph2_L,N,1);
                        vecZ2 = reshape(obj.U_D{2*c-1,1+f}.z(:,:,2:Nz+1),N,1) .* dS_rho_D_ph1_z .* reshape(dDTz_ph1_R,N,1) + reshape(obj.U_D{2*c,1+f}.z(:,:,2:Nz+1),N,1) .* dS_rho_D_ph2_z .* reshape(dDTz_ph2_R,N,1);
                        
                        DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2 + vecY2 + vecX2 - vecY1 - vecX1 - vecZ1, vecX1, vecY1, vecZ1];
                        DiagIndx = [-Grid.Nx*Grid.Ny, -Grid.Nx, -1, 0, 1, Grid.Nx, Grid.Nx*Grid.Ny];
                        Jz(:, 1+N*(j-1):N*j) = Jz(:, 1+N*(j-1):N*j) + spdiags(DiagVecs, DiagIndx, Grid.N, Grid.N) + obj.Dph{2*c-1,1+f} * spdiags(obj.dxdz(:,2*c-1,j), 0, N, N) + obj.Dph{2*c,1+f} * spdiags(obj.dxdz(:,2*c), 0, N, N);

                        vecX1 = reshape(obj.U_D{2*c-1,1+f}.x(1:Nx,:,:),N,1).* obj.drhodz(Index.Start:Index.End,1,j) / 2 .* reshape(Grid.Diffx(1:Nx,:,:,1),N,1) + reshape(obj.U_D{2*c,1+f}.x(1:Nx,:,:),N,1) .* obj.drhodz(Index.Start:Index.End,2,j) / 2  .* reshape(Grid.Diffx(1:Nx,:,:,2),N,1);
                        vecX2 = reshape(obj.U_D{2*c-1,1+f}.x(2:Nx+1,:,:),N,1).* obj.drhodz(Index.Start:Index.End,1,j) / 2 .* reshape(Grid.Diffx(2:Nx+1,:,:,1),N,1) + reshape(obj.U_D{2*c,1+f}.x(2:Nx+1,:,:),N,1) .* obj.drhodz(Index.Start:Index.End,2,j) / 2  .* reshape(Grid.Diffx(2:Nx+1,:,:,2),N,1);
                        vecY1 = reshape(obj.U_D{2*c-1,1+f}.y(:,1:Ny,:),N,1).* obj.drhodz(Index.Start:Index.End,1,j) / 2.* reshape(Grid.Diffy(:,1:Ny,:,1),N,1) + reshape(obj.U_D{2*c,1+f}.y(:,1:Ny,:),N,1) .* obj.drhodz(Index.Start:Index.End,2,j) / 2  .* reshape(Grid.Diffy(:,1:Ny,:,2),N,1);
                        vecY2 = reshape(obj.U_D{2*c-1,1+f}.y(:,2:Ny+1,:),N,1).* obj.drhodz(Index.Start:Index.End,1,j) / 2 .* reshape(Grid.Diffy(:,2:Ny+1,:,1),N,1) + reshape(obj.U_D{2*c,1+f}.y(:,2:Ny+1,:),N,1) .* obj.drhodz(Index.Start:Index.End,2,j) / 2  .* reshape(Grid.Diffy(:,2:Ny+1,:,2),N,1);
                        vecZ1 = reshape(obj.U_D{2*c-1,1+f}.z(:,:,1:Nz),N,1) .* obj.drhodz(Index.Start:Index.End,1,j) / 2 .* reshape(Grid.Diffz(:,:,1:Nz,1),N,1) + reshape(obj.U_D{2*c,1+f}.z(:,:,1:Nz),N,1) .* obj.drhodz(Index.Start:Index.End,2,j) / 2  .* reshape(Grid.Diffz(:,:,1:Nz,2),N,1);
                        vecZ2 = reshape(obj.U_D{2*c-1,1+f}.z(:,:,2:Nz+1),N,1) .* obj.drhodz(Index.Start:Index.End,1,j) / 2 .* reshape(Grid.Diffz(:,:,2:Nz+1,1),N,1) + reshape(obj.U_D{2*c,1+f}.z(:,:,2:Nz+1),N,1) .* obj.drhodz(Index.Start:Index.End,2,j) / 2  .* reshape(Grid.Diffz(:,:,2:Nz+1,2),N,1);

                        DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2 + vecY2 + vecX2 - vecY1 - vecX1 - vecZ1, vecX1, vecY1, vecZ1];
                        DiagIndx = [-Grid.Nx*Grid.Ny, -Grid.Nx, -1, 0, 1, Grid.Nx, Grid.Nx*Grid.Ny];
                        Jz(:, 1+N*(j-1):N*j) = Jz(:, 1+N*(j-1):N*j) + spdiags(DiagVecs, DiagIndx, Grid.N, Grid.N);

                        % Capillarity
                        % Jz{i,j} = Jz{i,j} - obj.Tph{i,1} * spdiags(obj.dPc, 0, N, N);
                        Jz(:, 1+N*(j-1):N*j) = Jz(:, 1+N*(j-1):N*j) - obj.Tph{2*c-1,1+f} * spdiags(obj.dPc(Index.Start:Index.End), 0, N, N);
                    end

                        % 3. temperature part
                        % JT
                        dMupxPh1 = reshape(obj.UpWind{1,1+f}.x * (obj.Mob(Index.Start:Index.End, 1) .* obj.x(Index.Start:Index.End,(c-1)*2+1) .* obj.drhodT(Index.Start:Index.End,1) + obj.Mob(Index.Start:Index.End, 1) .* obj.dxdT(Index.Start:Index.End,(c-1)*2+1) .* obj.rho(Index.Start:Index.End,1) + obj.dMobdT(Index.Start:Index.End, 1) .* obj.x(Index.Start:Index.End,(c-1)*2+1) .* obj.rho(Index.Start:Index.End,1)), Nx, Ny, Nz);
                        dMupyPh1 = reshape(obj.UpWind{1,1+f}.y * (obj.Mob(Index.Start:Index.End, 1) .* obj.x(Index.Start:Index.End,(c-1)*2+1) .* obj.drhodT(Index.Start:Index.End,1) + obj.Mob(Index.Start:Index.End, 1) .* obj.dxdT(Index.Start:Index.End,(c-1)*2+1) .* obj.rho(Index.Start:Index.End,1) + obj.dMobdT(Index.Start:Index.End, 1) .* obj.x(Index.Start:Index.End,(c-1)*2+1) .* obj.rho(Index.Start:Index.End,1)), Nx, Ny, Nz);
                        dMupzPh1 = reshape(obj.UpWind{1,1+f}.z * (obj.Mob(Index.Start:Index.End, 1) .* obj.x(Index.Start:Index.End,(c-1)*2+1) .* obj.drhodT(Index.Start:Index.End,1) + obj.Mob(Index.Start:Index.End, 1) .* obj.dxdT(Index.Start:Index.End,(c-1)*2+1) .* obj.rho(Index.Start:Index.End,1) + obj.dMobdT(Index.Start:Index.End, 1) .* obj.x(Index.Start:Index.End,(c-1)*2+1) .* obj.rho(Index.Start:Index.End,1)), Nx, Ny, Nz);
                        dMupxPh2 = reshape(obj.UpWind{2,1+f}.x * (obj.Mob(Index.Start:Index.End, 2) .* obj.x(Index.Start:Index.End,(c-1)*2+2) .* obj.drhodT(Index.Start:Index.End,2) + obj.Mob(Index.Start:Index.End, 2) .* obj.dxdT(Index.Start:Index.End,(c-1)*2+2) .* obj.rho(Index.Start:Index.End,2) + obj.dMobdT(Index.Start:Index.End, 2) .* obj.x(Index.Start:Index.End,(c-1)*2+2) .* obj.rho(Index.Start:Index.End,2)), Nx, Ny, Nz);
                        dMupyPh2 = reshape(obj.UpWind{2,1+f}.y * (obj.Mob(Index.Start:Index.End, 2) .* obj.x(Index.Start:Index.End,(c-1)*2+2) .* obj.drhodT(Index.Start:Index.End,2) + obj.Mob(Index.Start:Index.End, 2) .* obj.dxdT(Index.Start:Index.End,(c-1)*2+2) .* obj.rho(Index.Start:Index.End,2) + obj.dMobdT(Index.Start:Index.End, 2) .* obj.x(Index.Start:Index.End,(c-1)*2+2) .* obj.rho(Index.Start:Index.End,2)), Nx, Ny, Nz);
                        dMupzPh2 = reshape(obj.UpWind{2,1+f}.z * (obj.Mob(Index.Start:Index.End, 2) .* obj.x(Index.Start:Index.End,(c-1)*2+2) .* obj.drhodT(Index.Start:Index.End,2) + obj.Mob(Index.Start:Index.End, 2) .* obj.dxdT(Index.Start:Index.End,(c-1)*2+2) .* obj.rho(Index.Start:Index.End,2) + obj.dMobdT(Index.Start:Index.End, 2) .* obj.x(Index.Start:Index.End,(c-1)*2+2) .* obj.rho(Index.Start:Index.End,2)), Nx, Ny, Nz);

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
                        acc = pv/dt .* (obj.drhoTdT(Index.Start:Index.End) .* obj.z(Index.Start:Index.End,c));
                        DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2 + vecY2+vecX2-vecY1-vecX1-vecZ1+acc, vecX1, vecY1, vecZ1];
                        DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
                        JT = JT + spdiags(DiagVecs, DiagIndx, N, N);

                        dS_rho_D_ph1_x = (obj.dSdT(Index.Start:Index.End,1) .*  Medium.Diff(:,1) .* Grid.Por) .* Grid.Ax / Grid.dx;
                        dS_rho_D_ph1_y = (obj.dSdT(Index.Start:Index.End,1) .*  Medium.Diff(:,1) .* Grid.Por) .* Grid.Ay / Grid.dy;
                        dS_rho_D_ph1_z = (obj.dSdT(Index.Start:Index.End,1) .*  Medium.Diff(:,1) .* Grid.Por) .* Grid.Az / Grid.dz;
                        dS_rho_D_ph2_x = (-obj.dSdT(Index.Start:Index.End,1) .*  Medium.Diff(:,2) .* Grid.Por) .* Grid.Ax / Grid.dx;
                        dS_rho_D_ph2_y = (-obj.dSdT(Index.Start:Index.End,1) .*  Medium.Diff(:,2) .* Grid.Por) .* Grid.Ay / Grid.dy;
                        dS_rho_D_ph2_z = (-obj.dSdT(Index.Start:Index.End,1) .*  Medium.Diff(:,2) .* Grid.Por) .* Grid.Az / Grid.dz;

                        vecX1 = reshape(obj.U_D{2*c-1,1+f}.x(1:Nx,:,:),N,1).* dS_rho_D_ph1_x .* reshape(dDTx_ph1_L,N,1) + reshape(obj.U_D{2*c,1+f}.x(1:Nx,:,:),N,1) .* dS_rho_D_ph2_x .* reshape(dDTx_ph2_L,N,1);
                        vecX2 = reshape(obj.U_D{2*c-1,1+f}.x(2:Nx+1,:,:),N,1).* dS_rho_D_ph1_x .* reshape(dDTx_ph1_R,N,1) + reshape(obj.U_D{2*c,1+f}.x(2:Nx+1,:,:),N,1) .* dS_rho_D_ph2_x .* reshape(dDTx_ph2_R,N,1);
                        vecY1 = reshape(obj.U_D{2*c-1,1+f}.y(:,1:Ny,:),N,1).* dS_rho_D_ph1_y .* reshape(dDTy_ph1_L,N,1) + reshape(obj.U_D{2*c,1+f}.y(:,1:Ny,:),N,1) .* dS_rho_D_ph2_y .* reshape(dDTy_ph2_L,N,1);
                        vecY2 = reshape(obj.U_D{2*c-1,1+f}.y(:,2:Ny+1,:),N,1).* dS_rho_D_ph1_y .* reshape(dDTy_ph1_R,N,1) + reshape(obj.U_D{2*c,1+f}.y(:,2:Ny+1,:),N,1) .* dS_rho_D_ph2_y .* reshape(dDTy_ph2_R,N,1);
                        vecZ1 = reshape(obj.U_D{2*c-1,1+f}.z(:,:,1:Nz),N,1) .* dS_rho_D_ph1_z .* reshape(dDTz_ph1_L,N,1) + reshape(obj.U_D{2*c,1+f}.z(:,:,1:Nz),N,1) .* dS_rho_D_ph2_z .* reshape(dDTz_ph2_L,N,1);
                        vecZ2 = reshape(obj.U_D{2*c-1,1+f}.z(:,:,2:Nz+1),N,1) .* dS_rho_D_ph1_z .* reshape(dDTz_ph1_R,N,1) + reshape(obj.U_D{2*c,1+f}.z(:,:,2:Nz+1),N,1) .* dS_rho_D_ph2_z .* reshape(dDTz_ph2_R,N,1);

                        DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2 + vecY2 + vecX2 - vecY1 - vecX1 - vecZ1, vecX1, vecY1, vecZ1];
                        DiagIndx = [-Grid.Nx*Grid.Ny, -Grid.Nx, -1, 0, 1, Grid.Nx, Grid.Nx*Grid.Ny];
                        JT = JT + spdiags(DiagVecs, DiagIndx, Grid.N, Grid.N) + obj.Dph{2*c-1,1+f} * spdiags(obj.dxdT(:,2*c-1), 0, N, N) + obj.Dph{2*c,1+f} * spdiags(obj.dxdT(:,2*c), 0, N, N);

                        vecX1 = reshape(obj.U_D{2*c-1,1+f}.x(1:Nx,:,:),N,1).* obj.drhodT(Index.Start:Index.End,1) / 2 .* reshape(Grid.Diffx(1:Nx,:,:,1),N,1) + reshape(obj.U_D{2*c,1+f}.x(1:Nx,:,:),N,1) .* obj.drhodT(Index.Start:Index.End,2) / 2  .* reshape(Grid.Diffx(1:Nx,:,:,2),N,1);
                        vecX2 = reshape(obj.U_D{2*c-1,1+f}.x(2:Nx+1,:,:),N,1).* obj.drhodT(Index.Start:Index.End,1) / 2 .* reshape(Grid.Diffx(2:Nx+1,:,:,1),N,1) + reshape(obj.U_D{2*c,1+f}.x(2:Nx+1,:,:),N,1) .* obj.drhodT(Index.Start:Index.End,2) / 2  .* reshape(Grid.Diffx(2:Nx+1,:,:,2),N,1);
                        vecY1 = reshape(obj.U_D{2*c-1,1+f}.y(:,1:Ny,:),N,1).* obj.drhodT(Index.Start:Index.End,1) / 2.* reshape(Grid.Diffy(:,1:Ny,:,1),N,1) + reshape(obj.U_D{2*c,1+f}.y(:,1:Ny,:),N,1) .* obj.drhodT(Index.Start:Index.End,2) / 2  .* reshape(Grid.Diffy(:,1:Ny,:,2),N,1);
                        vecY2 = reshape(obj.U_D{2*c-1,1+f}.y(:,2:Ny+1,:),N,1).* obj.drhodT(Index.Start:Index.End,1) / 2 .* reshape(Grid.Diffy(:,2:Ny+1,:,1),N,1) + reshape(obj.U_D{2*c,1+f}.y(:,2:Ny+1,:),N,1) .* obj.drhodT(Index.Start:Index.End,2) / 2  .* reshape(Grid.Diffy(:,2:Ny+1,:,2),N,1);
                        vecZ1 = reshape(obj.U_D{2*c-1,1+f}.z(:,:,1:Nz),N,1) .* obj.drhodT(Index.Start:Index.End,1) / 2 .* reshape(Grid.Diffz(:,:,1:Nz,1),N,1) + reshape(obj.U_D{2*c,1+f}.z(:,:,1:Nz),N,1) .* obj.drhodT(Index.Start:Index.End,2) / 2  .* reshape(Grid.Diffz(:,:,1:Nz,2),N,1);
                        vecZ2 = reshape(obj.U_D{2*c-1,1+f}.z(:,:,2:Nz+1),N,1) .* obj.drhodT(Index.Start:Index.End,1) / 2 .* reshape(Grid.Diffz(:,:,2:Nz+1,1),N,1) + reshape(obj.U_D{2*c,1+f}.z(:,:,2:Nz+1),N,1) .* obj.drhodT(Index.Start:Index.End,2) / 2  .* reshape(Grid.Diffz(:,:,2:Nz+1,2),N,1);

                        DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2 + vecY2 + vecX2 - vecY1 - vecX1 - vecZ1, vecX1, vecY1, vecZ1];
                        DiagIndx = [-Grid.Nx*Grid.Ny, -Grid.Nx, -1, 0, 1, Grid.Nx, Grid.Nx*Grid.Ny];
                        JT = JT + spdiags(DiagVecs, DiagIndx, Grid.N, Grid.N);

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
                        dMupPh1 = obj.UpWind{1,1+f} * ( obj.dMobdz(Index.Start:Index.End, 1, j) .* obj.x(Index.Start:Index.End,(c-1)*2+1) .* obj.rho(Index.Start:Index.End,1) + obj.Mob(Index.Start:Index.End, 1) .* obj.dxdz(Index.Start:Index.End,(c-1)*2+1, j) .* obj.rho(Index.Start:Index.End,1) + obj.Mob(Index.Start:Index.End, 1) .* obj.x(Index.Start:Index.End, (c-1)*2+1) .* obj.drhodz(Index.Start:Index.End, 1));
                        dMupPh2 = obj.UpWind{2,1+f} * ( obj.dMobdz(Index.Start:Index.End, 2, j) .* obj.x(Index.Start:Index.End,(c-1)*2+2) .* obj.rho(Index.Start:Index.End,2) + obj.Mob(Index.Start:Index.End, 2) .* obj.dxdz(Index.Start:Index.End,(c-1)*2+2, j) .* obj.rho(Index.Start:Index.End,2) + obj.Mob(Index.Start:Index.End, 2) .* obj.x(Index.Start:Index.End, (c-1)*2+2) .* obj.drhodz(Index.Start:Index.End, 2));

                        FluxDerivativePh1 = dMupPh1 .* abs(obj.U{1,1+f});
                        FluxDerivativePh2 = dMupPh2 .* abs(obj.U{2,1+f});

                        if c == obj.NofComponents
                            acc = - pv/dt .* rhoT + pv/dt .* obj.z(Index.Start:Index.End,c) .* obj.drhoTdz(Index.Start:Index.End, j);
                        elseif c == j
                            acc = pv/dt .* rhoT + pv/dt .* obj.z(Index.Start:Index.End,c) .* obj.drhoTdz(Index.Start:Index.End, j);
                        end
                        Jz = sparse(1:N, 1:N, acc)...
                            + (UpwindPermutation1Ph1 * sparse(1:N_Face, 1:N_Face, FluxDerivativePh1) * UpwindPermutation2Ph1')'...
                            + (UpwindPermutation1Ph2 * sparse(1:N_Face, 1:N_Face, FluxDerivativePh2) * UpwindPermutation2Ph2')';
                        % Capillarity
                        Jz = Jz - obj.Tph{2*c-1,1+f} * sparse(1:N, 1:N, obj.dPc(Index.Start:Index.End));
                    end
                    % 3. temperature part
                    dMupPh1 = obj.UpWind{1,1+f} * (obj.Mob(Index.Start:Index.End, 1) .* obj.x(Index.Start:Index.End,(c-1)*2+1) .* obj.drhodT(Index.Start:Index.End,1) + obj.Mob(Index.Start:Index.End, 1) .* obj.dxdT(Index.Start:Index.End,(c-1)*2+1) .* obj.rho(Index.Start:Index.End,1) + obj.dMobdT(Index.Start:Index.End, 1) .* obj.x(Index.Start:Index.End,(c-1)*2+1) .* obj.rho(Index.Start:Index.End,1));
                    dMupPh2 = obj.UpWind{2,1+f} * (obj.Mob(Index.Start:Index.End, 2) .* obj.x(Index.Start:Index.End,(c-1)*2+2) .* obj.drhodT(Index.Start:Index.End,2) + obj.Mob(Index.Start:Index.End, 2) .* obj.dxdT(Index.Start:Index.End,(c-1)*2+2) .* obj.rho(Index.Start:Index.End,2) + obj.dMobdT(Index.Start:Index.End, 2) .* obj.x(Index.Start:Index.End,(c-1)*2+2) .* obj.rho(Index.Start:Index.End,2));

                    FluxDerivativePh1 = dMupPh1 .* abs(obj.U{1,1+f}); % nf x 1
                    FluxDerivativePh2 = dMupPh2 .* abs(obj.U{2,1+f});

                    acc = pv/dt .* (obj.drhoTdT(Index.Start:Index.End) .* obj.z(Index.Start:Index.End,c));
                    JT = JT + sparse(1:N, 1:N, acc)...
                        + (UpwindPermutation1Ph1 * sparse(1:N_Face, 1:N_Face, FluxDerivativePh1) * UpwindPermutation2Ph1')'... % ncxnf * nfxnf * nfxnc = ncxnc
                        + (UpwindPermutation1Ph2 * sparse(1:N_Face, 1:N_Face, FluxDerivativePh2) * UpwindPermutation2Ph2')';

                otherwise
                    error("Discretization not found.")
            end
        end

        %%
        function [Jp, Jz, JT] = BuildEnergyMediumJacobian(obj, Medium, Grid, dt, Index, f)

            N = Grid.N;
            N_Face = size(obj.UpWind{1,1+f}, 1);
            pv = Grid.Volume .* Medium.Por;
            rhoT = Medium.State.Properties('rhoT').Value;
            
            cond = zeros(N,obj.NofPhases);
            for i=1:obj.NofPhases
                lamda = Medium.State.Properties(['cond_',num2str(i)]);
                cond(:,i) = lamda.Value;
            end

            % BUILD FIM JACOBIAN BLOCK BY BLOCK
            Jz = sparse(N, N * (obj.NofComponents-1));
            JT = sparse(N, N);
            % 1.a: divergence
            Jp = obj.Thph{1,1+f}  + obj.Thph{2, 1+f};
            switch class(Grid)
                case 'unstructured_grid'
                    
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

                    %% 1.J_EB_P Block

                    dMupPh1 = obj.UpWind{1,1+f} * (obj.Mob(Index.Start:Index.End, 1) .* obj.h(Index.Start:Index.End,1) .* obj.drhodp(Index.Start:Index.End,1) + obj.Mob(Index.Start:Index.End, 1) .* obj.dhdp(Index.Start:Index.End,1) .* obj.rho(Index.Start:Index.End,1) + obj.dMobdp(Index.Start:Index.End, 1) .* obj.h(Index.Start:Index.End,1) .* obj.rho(Index.Start:Index.End,1));
                    dMupPh2 = obj.UpWind{2,1+f} * (obj.Mob(Index.Start:Index.End, 2) .* obj.h(Index.Start:Index.End,2) .* obj.drhodp(Index.Start:Index.End,2) + obj.Mob(Index.Start:Index.End, 2) .* obj.dhdp(Index.Start:Index.End,2) .* obj.rho(Index.Start:Index.End,2) + obj.dMobdp(Index.Start:Index.End, 2) .* obj.h(Index.Start:Index.End,2) .* obj.rho(Index.Start:Index.End,2));

                    FluxDerivativePh1 = dMupPh1 .* abs(obj.U{1,1+f}); % nf x 1
                    FluxDerivativePh2 = dMupPh2 .* abs(obj.U{2,1+f});

                    acc_1 = obj.rho(Index.Start:Index.End,1).*obj.h(Index.Start:Index.End,1) - obj.P(Index.Start:Index.End,1);
                    acc_2 = obj.rho(Index.Start:Index.End,2).*obj.h(Index.Start:Index.End,2) - obj.P(Index.Start:Index.End,2);
                    dacc_1 = obj.drhodp(Index.Start:Index.End,1).*obj.h(Index.Start:Index.End,1) + obj.rho(Index.Start:Index.End,1).*obj.dhdp(Index.Start:Index.End,1) - 1;
                    dacc_2 = obj.drhodp(Index.Start:Index.End,2).*obj.h(Index.Start:Index.End,2) + obj.rho(Index.Start:Index.End,2).*obj.dhdp(Index.Start:Index.End,2) - 1;

                    acc = pv /dt .* (obj.dSdp .* acc_1 + obj.S(Index.Start:Index.End,1) .* dacc_1 + obj.S(Index.Start:Index.End,2) .* dacc_2 - obj.dSdp .* acc_2);
                    Jp = Jp + sparse(1:N, 1:N, acc)...
                        + (UpwindPermutation1Ph1 * sparse(1:N_Face, 1:N_Face, FluxDerivativePh1) * UpwindPermutation2Ph1')'... % ncxnf * nfxnf * nfxnc = ncxnc
                        + (UpwindPermutation1Ph2 * sparse(1:N_Face, 1:N_Face, FluxDerivativePh2) * UpwindPermutation2Ph2')';

                    CellNeighbor1Index = Grid.CornerPointGridData.Internal_Faces.CellNeighbor1Index;
                    CellNeighbor2Index = Grid.CornerPointGridData.Internal_Faces.CellNeighbor2Index;
                    CellNeighbor1Vec   = Grid.CornerPointGridData.Internal_Faces.CellNeighbor1Vec;
                    CellNeighbor2Vec   = Grid.CornerPointGridData.Internal_Faces.CellNeighbor2Vec;
                    Nvec               = Grid.CornerPointGridData.Internal_Faces.Nvec;

                    NeighborPermutation1 = sparse(repmat((1:nf)',1,2), [CellNeighbor1Index, CellNeighbor2Index], [ones(nf,1), zeros(nf,1)], nf, nc)';   % nc x nf
                    NeighborPermutation2 = sparse(repmat((1:nf)',1,2), [CellNeighbor1Index, CellNeighbor2Index], [zeros(nf,1), ones(nf,1)], nf, nc)';   % nc x nf

                    S_D_ph1  = obj.S(Index.Start:Index.End,1) .* cond(:,1) .* Grid.Por;
                    S_D_ph2  = obj.S(Index.Start:Index.End,2) .* cond(:,2) .* Grid.Por;
                    dS_D_ph1 =  obj.dSdp(Index.Start:Index.End,1) .* cond(:,1) .* Grid.Por;
                    dS_D_ph2 = -obj.dSdp(Index.Start:Index.End,1) .* cond(:,2) .* Grid.Por;

                    Half1_ph1 = S_D_ph1(CellNeighbor1Index) .* abs(sum(CellNeighbor1Vec .* Nvec, 2) ./ (sum(CellNeighbor1Vec .* CellNeighbor1Vec, 2) + 1e-30));
                    Half2_ph1 = S_D_ph1(CellNeighbor2Index) .* abs(sum(CellNeighbor2Vec .* Nvec, 2) ./ (sum(CellNeighbor2Vec .* CellNeighbor2Vec, 2) + 1e-30));
                    Half1_ph2 = S_D_ph2(CellNeighbor1Index) .* abs(sum(CellNeighbor1Vec .* Nvec, 2) ./ (sum(CellNeighbor1Vec .* CellNeighbor1Vec, 2) + 1e-30));
                    Half2_ph2 = S_D_ph2(CellNeighbor2Index) .* abs(sum(CellNeighbor2Vec .* Nvec, 2) ./ (sum(CellNeighbor2Vec .* CellNeighbor2Vec, 2) + 1e-30));

                    % keep the same L/R convention as your Cartesian code
                    dDT_ph1_L = abs(sum(CellNeighbor1Vec .* Nvec, 2) ./ (sum(CellNeighbor1Vec .* CellNeighbor1Vec, 2) + 1e-30)) .* Half1_ph1.^2 ./ (Half1_ph1 + Half2_ph1 + 1e-30).^2;
                    dDT_ph1_R = abs(sum(CellNeighbor2Vec .* Nvec, 2) ./ (sum(CellNeighbor2Vec .* CellNeighbor2Vec, 2) + 1e-30)) .* Half2_ph1.^2 ./ (Half1_ph1 + Half2_ph1 + 1e-30).^2;
                    dDT_ph2_L = abs(sum(CellNeighbor1Vec .* Nvec, 2) ./ (sum(CellNeighbor1Vec .* CellNeighbor1Vec, 2) + 1e-30)) .* Half1_ph2.^2 ./ (Half1_ph2 + Half2_ph2 + 1e-30).^2;
                    dDT_ph2_R = abs(sum(CellNeighbor2Vec .* Nvec, 2) ./ (sum(CellNeighbor2Vec .* CellNeighbor2Vec, 2) + 1e-30)) .* Half2_ph2.^2 ./ (Half1_ph2 + Half2_ph2 + 1e-30).^2;

                    vecL = obj.U_T{1,1+f} .* dDT_ph1_L .* dS_D_ph1(CellNeighbor1Index) + obj.U_T{1,1+f} .* dDT_ph2_L .* dS_D_ph2(CellNeighbor1Index);
                    vecR = obj.U_T{1,1+f} .* dDT_ph1_R .* dS_D_ph1(CellNeighbor2Index) + obj.U_T{1,1+f} .* dDT_ph2_R .* dS_D_ph2(CellNeighbor2Index);

                    Jp = Jp  + (NeighborPermutation1 * sparse(1:nf, 1:nf, vecL, nf, nf) * Grid.ConnectivityMatrix')' ...
                        + (NeighborPermutation2 * sparse(1:nf, 1:nf, vecR, nf, nf) * Grid.ConnectivityMatrix')';
                    
                    %% 2. J_EB_Z Block

                    for j=1:obj.NofComponents-1
                        dMupPh1 = obj.UpWind{1,1+f} * ( obj.dMobdz(Index.Start:Index.End, 1, j) .* obj.h(Index.Start:Index.End,1) .* obj.rho(Index.Start:Index.End,1) + obj.Mob(Index.Start:Index.End, 1) .* obj.dhdz(Index.Start:Index.End, 1, j) .* obj.rho(Index.Start:Index.End,1) + obj.Mob(Index.Start:Index.End, 1) .* obj.h(Index.Start:Index.End, 1) .* obj.drhodz(Index.Start:Index.End, 1));
                        dMupPh2 = obj.UpWind{2,1+f} * ( obj.dMobdz(Index.Start:Index.End, 2, j) .* obj.h(Index.Start:Index.End,2) .* obj.rho(Index.Start:Index.End,2) + obj.Mob(Index.Start:Index.End, 2) .* obj.dhdz(Index.Start:Index.End, 2, j) .* obj.rho(Index.Start:Index.End,2) + obj.Mob(Index.Start:Index.End, 2) .* obj.h(Index.Start:Index.End, 2) .* obj.drhodz(Index.Start:Index.End, 2));

                        FluxDerivativePh1 = dMupPh1 .* abs(obj.U{1,1+f});
                        FluxDerivativePh2 = dMupPh2 .* abs(obj.U{2,1+f});

                        dacc_1 = obj.drhodz(Index.Start:Index.End,1,j).*obj.h(Index.Start:Index.End,1) + obj.rho(Index.Start:Index.End,1).*obj.dhdz(Index.Start:Index.End,1,j) - obj.dPc(Index.Start:Index.End,j);
                        dacc_2 = obj.drhodz(Index.Start:Index.End,2,j).*obj.h(Index.Start:Index.End,2) + obj.rho(Index.Start:Index.End,2).*obj.dhdz(Index.Start:Index.End,2,j);

                        acc = pv /dt .* (obj.dSdz(Index.Start:Index.End,1,j) .* acc_1 + obj.S(Index.Start:Index.End,1) .* dacc_1 + obj.S(Index.Start:Index.End,2) .* dacc_2 - obj.dSdz(Index.Start:Index.End,1,j) .* acc_2);

                        Jz(:, 1+N*(j-1):N*j) = sparse(1:N, 1:N, acc)...
                            + (UpwindPermutation1Ph1 * sparse(1:N_Face, 1:N_Face, FluxDerivativePh1) * UpwindPermutation2Ph1')'...
                            + (UpwindPermutation1Ph2 * sparse(1:N_Face, 1:N_Face, FluxDerivativePh2) * UpwindPermutation2Ph2')';

                        dS_D_ph1 =  obj.dSdz(Index.Start:Index.End,1,j) .* cond(:,1) .* Grid.Por;
                        dS_D_ph2 = -obj.dSdz(Index.Start:Index.End,1,j) .* cond(:,2) .* Grid.Por;

                        vecL = obj.U_T{1,1+f} .* dDT_ph1_L .* dS_D_ph1(CellNeighbor1Index) + obj.U_T{1,1+f} .* dDT_ph2_L .* dS_D_ph2(CellNeighbor1Index);
                        vecR = obj.U_T{1,1+f} .* dDT_ph1_R .* dS_D_ph1(CellNeighbor2Index) + obj.U_T{1,1+f} .* dDT_ph2_R .* dS_D_ph2(CellNeighbor2Index);
                        
                        Jz(:,1+N*(j-1):N*j) = Jz(:,1+N*(j-1):N*j) + (NeighborPermutation1 * sparse(1:nf, 1:nf, vecL, nf, nf) * Grid.ConnectivityMatrix')' + (NeighborPermutation2 * sparse(1:nf, 1:nf, vecR, nf, nf) * Grid.ConnectivityMatrix')';

                        % Capillarity
                        Jz(:, 1+N*(j-1):N*j) = Jz(:, 1+N*(j-1):N*j) - obj.Thph{1,1+f} * sparse(1:N, 1:N, obj.dPc(Index.Start:Index.End));
                    end

                    %% 3. J_EB_T Block

                    dMupPh1 = obj.UpWind{1,1+f} * (obj.Mob(Index.Start:Index.End, 1) .* obj.h(Index.Start:Index.End,1) .* obj.drhodT(Index.Start:Index.End,1) + obj.Mob(Index.Start:Index.End, 1) .* obj.dhdT(Index.Start:Index.End,1) .* obj.rho(Index.Start:Index.End,1) + obj.dMobdT(Index.Start:Index.End, 1) .* obj.h(Index.Start:Index.End,1) .* obj.rho(Index.Start:Index.End,1));
                    dMupPh2 = obj.UpWind{2,1+f} * (obj.Mob(Index.Start:Index.End, 2) .* obj.h(Index.Start:Index.End,2) .* obj.drhodT(Index.Start:Index.End,2) + obj.Mob(Index.Start:Index.End, 2) .* obj.dhdT(Index.Start:Index.End,2) .* obj.rho(Index.Start:Index.End,2) + obj.dMobdT(Index.Start:Index.End, 2) .* obj.h(Index.Start:Index.End,2) .* obj.rho(Index.Start:Index.End,2));

                    FluxDerivativePh1 = dMupPh1 .* abs(obj.U{1,1+f}); % nf x 1
                    FluxDerivativePh2 = dMupPh2 .* abs(obj.U{2,1+f});

                    dacc_1 = obj.drhodT(Index.Start:Index.End,1).*obj.h(Index.Start:Index.End,1) + obj.rho(Index.Start:Index.End,1).*obj.dhdT(Index.Start:Index.End,1);
                    dacc_2 = obj.drhodT(Index.Start:Index.End,2).*obj.h(Index.Start:Index.End,2) + obj.rho(Index.Start:Index.End,2).*obj.dhdT(Index.Start:Index.End,2);
                    acc = pv /dt .* (obj.dSdT .* acc_1 + obj.S(Index.Start:Index.End,1) .* dacc_1 + obj.S(Index.Start:Index.End,2) .* dacc_2 - obj.dSdT .* acc_2) + Grid.Volume .* (1 - Medium.Por) / dt .* Medium.Rho .* Medium.Cpr;

                    JT = JT + sparse(1:N, 1:N, acc)...
                        + (UpwindPermutation1Ph1 * sparse(1:N_Face, 1:N_Face, FluxDerivativePh1) * UpwindPermutation2Ph1')'... % ncxnf * nfxnf * nfxnc = ncxnc
                        + (UpwindPermutation1Ph2 * sparse(1:N_Face, 1:N_Face, FluxDerivativePh2) * UpwindPermutation2Ph2')' + obj.Tk{1,1+f};

                    dS_D_ph1 =  obj.dSdT(Index.Start:Index.End,1) .* cond(:,1) .* Grid.Por;
                    dS_D_ph2 = -obj.dSdT(Index.Start:Index.End,1) .* cond(:,2) .* Grid.Por;

                    vecL = obj.U_T{1,1+f} .* dDT_ph1_L .* dS_D_ph1(CellNeighbor1Index) + obj.U_T{1,1+f} .* dDT_ph2_L .* dS_D_ph2(CellNeighbor1Index);
                    vecR = obj.U_T{1,1+f} .* dDT_ph1_R .* dS_D_ph1(CellNeighbor2Index) + obj.U_T{1,1+f} .* dDT_ph2_R .* dS_D_ph2(CellNeighbor2Index);

                    JT = JT + (NeighborPermutation1 * sparse(1:nf, 1:nf, vecL, nf, nf) * Grid.ConnectivityMatrix')' ...
                        + (NeighborPermutation2 * sparse(1:nf, 1:nf, vecR, nf, nf) * Grid.ConnectivityMatrix')';

                case('cartesian_grid')
                    Nx = Grid.Nx;
                    Ny = Grid.Ny;
                    Nz = Grid.Nz;

                    Tx_ph1 = zeros(Nx+1, Ny, Nz);
                    Ty_ph1 = zeros(Nx, Ny+1, Nz);
                    Tz_ph1 = zeros(Nx, Ny, Nz+1);
                    Tx_ph2 = zeros(Nx+1, Ny, Nz);
                    Ty_ph2 = zeros(Nx, Ny+1, Nz);
                    Tz_ph2 = zeros(Nx, Ny, Nz+1);

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

                    %% 1. J_EB_P Block
                    dMupxPh1 = reshape(obj.UpWind{1,1+f}.x * (obj.Mob(Index.Start:Index.End, 1) .* obj.h(Index.Start:Index.End,1) .* obj.drhodp(Index.Start:Index.End,1) + obj.Mob(Index.Start:Index.End, 1) .* obj.dhdp(Index.Start:Index.End,1) .* obj.rho(Index.Start:Index.End,1) + obj.dMobdp(Index.Start:Index.End, 1) .* obj.h(Index.Start:Index.End,1) .* obj.rho(Index.Start:Index.End,1)), Nx, Ny, Nz);
                    dMupyPh1 = reshape(obj.UpWind{1,1+f}.y * (obj.Mob(Index.Start:Index.End, 1) .* obj.h(Index.Start:Index.End,1) .* obj.drhodp(Index.Start:Index.End,1) + obj.Mob(Index.Start:Index.End, 1) .* obj.dhdp(Index.Start:Index.End,1) .* obj.rho(Index.Start:Index.End,1) + obj.dMobdp(Index.Start:Index.End, 1) .* obj.h(Index.Start:Index.End,1) .* obj.rho(Index.Start:Index.End,1)), Nx, Ny, Nz);
                    dMupzPh1 = reshape(obj.UpWind{1,1+f}.z * (obj.Mob(Index.Start:Index.End, 1) .* obj.h(Index.Start:Index.End,1) .* obj.drhodp(Index.Start:Index.End,1) + obj.Mob(Index.Start:Index.End, 1) .* obj.dhdp(Index.Start:Index.End,1) .* obj.rho(Index.Start:Index.End,1) + obj.dMobdp(Index.Start:Index.End, 1) .* obj.h(Index.Start:Index.End,1) .* obj.rho(Index.Start:Index.End,1)), Nx, Ny, Nz);
                    dMupxPh2 = reshape(obj.UpWind{2,1+f}.x * (obj.Mob(Index.Start:Index.End, 2) .* obj.h(Index.Start:Index.End,2) .* obj.drhodp(Index.Start:Index.End,2) + obj.Mob(Index.Start:Index.End, 2) .* obj.dhdp(Index.Start:Index.End,2) .* obj.rho(Index.Start:Index.End,2) + obj.dMobdp(Index.Start:Index.End, 2) .* obj.h(Index.Start:Index.End,2) .* obj.rho(Index.Start:Index.End,2)), Nx, Ny, Nz);
                    dMupyPh2 = reshape(obj.UpWind{2,1+f}.y * (obj.Mob(Index.Start:Index.End, 2) .* obj.h(Index.Start:Index.End,2) .* obj.drhodp(Index.Start:Index.End,2) + obj.Mob(Index.Start:Index.End, 2) .* obj.dhdp(Index.Start:Index.End,2) .* obj.rho(Index.Start:Index.End,2) + obj.dMobdp(Index.Start:Index.End, 2) .* obj.h(Index.Start:Index.End,2) .* obj.rho(Index.Start:Index.End,2)), Nx, Ny, Nz);
                    dMupzPh2 = reshape(obj.UpWind{2,1+f}.z * (obj.Mob(Index.Start:Index.End, 2) .* obj.h(Index.Start:Index.End,2) .* obj.drhodp(Index.Start:Index.End,2) + obj.Mob(Index.Start:Index.End, 2) .* obj.dhdp(Index.Start:Index.End,2) .* obj.rho(Index.Start:Index.End,2) + obj.dMobdp(Index.Start:Index.End, 2) .* obj.h(Index.Start:Index.End,2) .* obj.rho(Index.Start:Index.End,2)), Nx, Ny, Nz);

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
                    
                    acc_1 = obj.rho(Index.Start:Index.End,1).*obj.h(Index.Start:Index.End,1) - obj.P(Index.Start:Index.End,1);
                    acc_2 = obj.rho(Index.Start:Index.End,2).*obj.h(Index.Start:Index.End,2) - obj.P(Index.Start:Index.End,2);
                    dacc_1 = obj.drhodp(Index.Start:Index.End,1).*obj.h(Index.Start:Index.End,1) + obj.rho(Index.Start:Index.End,1).*obj.dhdp(Index.Start:Index.End,1) - 1;
                    dacc_2 = obj.drhodp(Index.Start:Index.End,2).*obj.h(Index.Start:Index.End,2) + obj.rho(Index.Start:Index.End,2).*obj.dhdp(Index.Start:Index.End,2) - 1;

                    acc = pv /dt .* (obj.dSdp .* acc_1 + obj.S(Index.Start:Index.End,1) .* dacc_1 + obj.S(Index.Start:Index.End,2) .* dacc_2 - obj.dSdp .* acc_2);
                    DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2 + vecY2+vecX2-vecY1-vecX1-vecZ1+acc, vecX1, vecY1, vecZ1];
                    DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
                    Jp = Jp + spdiags(DiagVecs, DiagIndx, N, N);

                    S_D_ph1 = reshape(obj.S(Index.Start:Index.End,1) .* cond(:,1) .* Grid.Por, Grid.Nx, Grid.Ny, Grid.Nz);
                    S_D_ph2 = reshape(obj.S(Index.Start:Index.End,2) .* cond(:,2) .* Grid.Por, Grid.Nx, Grid.Ny, Grid.Nz);

                    dS_D_ph1_x = (obj.dSdp(Index.Start:Index.End,1) .* cond(:,1) .* Grid.Por) .* Grid.Ax / Grid.dx;
                    dS_D_ph1_y = (obj.dSdp(Index.Start:Index.End,1) .* cond(:,1) .* Grid.Por) .* Grid.Ay / Grid.dy;
                    dS_D_ph1_z = (obj.dSdp(Index.Start:Index.End,1) .* cond(:,1) .* Grid.Por) .* Grid.Az / Grid.dz;
                    dS_D_ph2_x = (-obj.dSdp(Index.Start:Index.End,1) .* cond(:,2) .* Grid.Por) .* Grid.Ax / Grid.dx;
                    dS_D_ph2_y = (-obj.dSdp(Index.Start:Index.End,1) .* cond(:,2).* Grid.Por) .* Grid.Ay / Grid.dy;
                    dS_D_ph2_z = (-obj.dSdp(Index.Start:Index.End,1) .* cond(:,2) .* Grid.Por) .* Grid.Az / Grid.dz;

                    dDTx_ph1_L(2:Nx,:,:)   = 2 .* S_D_ph1(1:Nx-1,:,:) .^ 2 ./ (S_D_ph1(1:Nx-1,:,:) + S_D_ph1(2:Nx,:,:) + 1e-10).^2;
                    dDTx_ph1_R(1:Nx-1,:,:) = 2 .* S_D_ph1(2:Nx,:,:)   .^ 2 ./ (S_D_ph1(1:Nx-1,:,:) + S_D_ph1(2:Nx,:,:) + 1e-10).^2;
                    dDTy_ph1_L(:,2:Ny,:)   = 2 .* S_D_ph1(:,1:Ny-1,:) .^ 2 ./ (S_D_ph1(:,1:Ny-1,:) + S_D_ph1(:,2:Ny,:) + 1e-10).^2;
                    dDTy_ph1_R(:,1:Ny-1,:) = 2 .* S_D_ph1(:,2:Ny,:)   .^ 2 ./ (S_D_ph1(:,1:Ny-1,:) + S_D_ph1(:,2:Ny,:) + 1e-10).^2;
                    dDTz_ph1_L(:,:,2:Nz)   = 2 .* S_D_ph1(:,:,1:Nz-1) .^ 2 ./ (S_D_ph1(:,:,1:Nz-1) + S_D_ph1(:,:,2:Nz) + 1e-10).^2;
                    dDTz_ph1_R(:,:,1:Nz-1) = 2 .* S_D_ph1(:,:,2:Nz)   .^ 2 ./ (S_D_ph1(:,:,1:Nz-1) + S_D_ph1(:,:,2:Nz) + 1e-10).^2;

                    dDTx_ph2_L(2:Nx,:,:)   = 2 .* S_D_ph2(1:Nx-1,:,:) .^ 2 ./ (S_D_ph2(1:Nx-1,:,:) + S_D_ph2(2:Nx,:,:) + 1e-10).^2;
                    dDTx_ph2_R(1:Nx-1,:,:) = 2 .* S_D_ph2(2:Nx,:,:)   .^ 2 ./ (S_D_ph2(1:Nx-1,:,:) + S_D_ph2(2:Nx,:,:) + 1e-10).^2;
                    dDTy_ph2_L(:,2:Ny,:)   = 2 .* S_D_ph2(:,1:Ny-1,:) .^ 2 ./ (S_D_ph2(:,1:Ny-1,:) + S_D_ph2(:,2:Ny,:) + 1e-10).^2;
                    dDTy_ph2_R(:,1:Ny-1,:) = 2 .* S_D_ph2(:,2:Ny,:)   .^ 2 ./ (S_D_ph2(:,1:Ny-1,:) + S_D_ph2(:,2:Ny,:) + 1e-10).^2;
                    dDTz_ph2_L(:,:,2:Nz)   = 2 .* S_D_ph2(:,:,1:Nz-1) .^ 2 ./ (S_D_ph2(:,:,1:Nz-1) + S_D_ph2(:,:,2:Nz) + 1e-10).^2;
                    dDTz_ph2_R(:,:,1:Nz-1) = 2 .* S_D_ph2(:,:,2:Nz)   .^ 2 ./ (S_D_ph2(:,:,1:Nz-1) + S_D_ph2(:,:,2:Nz) + 1e-10).^2;

                    vecX1 = reshape(obj.U_T{1,1+f}.x(1:Nx,:,:),N,1).* dS_D_ph1_x .* reshape(dDTx_ph1_L,N,1) + reshape(obj.U_T{1,1+f}.x(1:Nx,:,:),N,1) .* dS_D_ph2_x .* reshape(dDTx_ph2_L,N,1);
                    vecX2 = reshape(obj.U_T{1,1+f}.x(2:Nx+1,:,:),N,1).* dS_D_ph1_x .* reshape(dDTx_ph1_R,N,1) + reshape(obj.U_T{1,1+f}.x(2:Nx+1,:,:),N,1) .* dS_D_ph2_x .* reshape(dDTx_ph2_R,N,1);
                    vecY1 = reshape(obj.U_T{1,1+f}.y(:,1:Ny,:),N,1).* dS_D_ph1_y .* reshape(dDTy_ph1_L,N,1) + reshape(obj.U_T{1,1+f}.y(:,1:Ny,:),N,1) .* dS_D_ph2_y .* reshape(dDTy_ph2_L,N,1);
                    vecY2 = reshape(obj.U_T{1,1+f}.y(:,2:Ny+1,:),N,1).* dS_D_ph1_y .* reshape(dDTy_ph1_R,N,1) + reshape(obj.U_T{1,1+f}.y(:,2:Ny+1,:),N,1) .* dS_D_ph2_y .* reshape(dDTy_ph2_R,N,1);
                    vecZ1 = reshape(obj.U_T{1,1+f}.z(:,:,1:Nz),N,1) .* dS_D_ph1_z .* reshape(dDTz_ph1_L,N,1) + reshape(obj.U_T{1,1+f}.z(:,:,1:Nz),N,1) .* dS_D_ph2_z .* reshape(dDTz_ph2_L,N,1);
                    vecZ2 = reshape(obj.U_T{1,1+f}.z(:,:,2:Nz+1),N,1) .* dS_D_ph1_z .* reshape(dDTz_ph1_R,N,1) + reshape(obj.U_T{1,1+f}.z(:,:,2:Nz+1),N,1) .* dS_D_ph2_z .* reshape(dDTz_ph2_R,N,1);

                    DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2 + vecY2 + vecX2 - vecY1 - vecX1 - vecZ1, vecX1, vecY1, vecZ1];
                    DiagIndx = [-Grid.Nx*Grid.Ny, -Grid.Nx, -1, 0, 1, Grid.Nx, Grid.Nx*Grid.Ny];
                    Jp = Jp + spdiags(DiagVecs, DiagIndx, Grid.N, Grid.N);

                    %% 2. J_EB_Z Block
                    for j=1:obj.NofComponents-1
                        dMupxPh1 = reshape(obj.UpWind{1,1+f}.x*...
                            ( obj.dMobdz(Index.Start:Index.End, 1, j) .* obj.h(Index.Start:Index.End,1) .* obj.rho(Index.Start:Index.End,1)...
                            + obj.Mob(Index.Start:Index.End, 1) .* obj.dhdz(Index.Start:Index.End,1, j) .* obj.rho(Index.Start:Index.End,1)...
                            + obj.Mob(Index.Start:Index.End, 1) .* obj.h(Index.Start:Index.End, 1) .* obj.drhodz(Index.Start:Index.End, 1, j)), Nx, Ny, Nz);
                        dMupyPh1 = reshape(obj.UpWind{1,1+f}.y*...
                            ( obj.dMobdz(Index.Start:Index.End, 1, j) .* obj.h(Index.Start:Index.End, 1) .* obj.rho(Index.Start:Index.End,1)...
                            + obj.Mob(Index.Start:Index.End, 1) .* obj.dhdz(Index.Start:Index.End, 1, j) .* obj.rho(Index.Start:Index.End,1)...
                            + obj.Mob(Index.Start:Index.End, 1) .* obj.h(Index.Start:Index.End,1) .* obj.drhodz(Index.Start:Index.End, 1, j)), Nx, Ny, Nz);
                        dMupzPh1 = reshape(obj.UpWind{1,1+f}.z*...
                            ( obj.dMobdz(Index.Start:Index.End, 1, j) .* obj.h(Index.Start:Index.End,1) .* obj.rho(Index.Start:Index.End,1)...
                            + obj.Mob(Index.Start:Index.End,1) .* obj.dhdz(Index.Start:Index.End,1, j) .* obj.rho(Index.Start:Index.End,1)...
                            + obj.Mob(Index.Start:Index.End,1) .* obj.h(Index.Start:Index.End,1) .* obj.drhodz(Index.Start:Index.End, 1, j)), Nx, Ny, Nz);
                        dMupxPh2 = reshape(obj.UpWind{2,1+f}.x*...
                            ( obj.dMobdz(Index.Start:Index.End, 2, j) .* obj.h(Index.Start:Index.End, 2) .* obj.rho(Index.Start:Index.End,2)...
                            + obj.Mob(Index.Start:Index.End,2) .* obj.dhdz(Index.Start:Index.End,2, j) .* obj.rho(Index.Start:Index.End,2)...
                            + obj.Mob(Index.Start:Index.End,2) .* obj.h(Index.Start:Index.End, 2) .* obj.drhodz(Index.Start:Index.End, 2, j)), Nx, Ny, Nz);
                        dMupyPh2 = reshape(obj.UpWind{2,1+f}.y*...
                            ( obj.dMobdz(Index.Start:Index.End, 2, j) .* obj.h(Index.Start:Index.End,2) .* obj.rho(Index.Start:Index.End,2)...
                            + obj.Mob(Index.Start:Index.End, 2) .* obj.dhdz(Index.Start:Index.End,2, j) .* obj.rho(Index.Start:Index.End,2)...
                            + obj.Mob(Index.Start:Index.End, 2) .* obj.h(Index.Start:Index.End, 2) .* obj.drhodz(Index.Start:Index.End, 2, j)), Nx, Ny, Nz);
                        dMupzPh2 = reshape(obj.UpWind{2,1+f}.z*...
                            ( obj.dMobdz(Index.Start:Index.End, 2, j) .* obj.h(Index.Start:Index.End, 2) .* obj.rho(Index.Start:Index.End,2)...
                            + obj.Mob(Index.Start:Index.End,2) .* obj.dhdz(Index.Start:Index.End, 2, j) .* obj.rho(Index.Start:Index.End,2)...
                            + obj.Mob(Index.Start:Index.End,2) .* obj.h(Index.Start:Index.End, 2) .* obj.drhodz(Index.Start:Index.End, 2, j)), Nx, Ny, Nz);

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
                        
                        dacc_1 = obj.drhodz(Index.Start:Index.End,1,j).*obj.h(Index.Start:Index.End,1) + obj.rho(Index.Start:Index.End,1).*obj.dhdz(Index.Start:Index.End,1,j) - obj.dPc(Index.Start:Index.End,j);
                        dacc_2 = obj.drhodz(Index.Start:Index.End,2,j).*obj.h(Index.Start:Index.End,2) + obj.rho(Index.Start:Index.End,2).*obj.dhdz(Index.Start:Index.End,2,j);

                        acc = pv /dt .* (obj.dSdz(Index.Start:Index.End,1,j) .* acc_1 + obj.S(Index.Start:Index.End,1) .* dacc_1 + obj.S(Index.Start:Index.End,2) .* dacc_2 - obj.dSdz(Index.Start:Index.End,1,j) .* acc_2);

                        DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2+vecY2+vecX2-vecZ1-vecY1-vecX1+acc, vecX1, vecY1, vecZ1];
                        DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
                        Jz(:, 1+N*(j-1):N*j) = spdiags(DiagVecs,DiagIndx, N, N);
                        
                        dS_D_ph1_x = (obj.dSdz(Index.Start:Index.End,1,j) .* cond(:,1) .* Grid.Por) .* Grid.Ax / Grid.dx;
                        dS_D_ph1_y = (obj.dSdz(Index.Start:Index.End,1,j) .* cond(:,1) .* Grid.Por) .* Grid.Ay / Grid.dy;
                        dS_D_ph1_z = (obj.dSdz(Index.Start:Index.End,1,j) .* cond(:,1) .* Grid.Por) .* Grid.Az / Grid.dz;
                        dS_D_ph2_x = (-obj.dSdz(Index.Start:Index.End,1,j) .* cond(:,2) .* Grid.Por) .* Grid.Ax / Grid.dx;
                        dS_D_ph2_y = (-obj.dSdz(Index.Start:Index.End,1,j) .* cond(:,2) .* Grid.Por) .* Grid.Ay / Grid.dy;
                        dS_D_ph2_z = (-obj.dSdz(Index.Start:Index.End,1,j) .* cond(:,2) .* Grid.Por) .* Grid.Az / Grid.dz;

                        vecX1 = reshape(obj.U_T{1,1+f}.x(1:Nx,:,:),N,1).* dS_D_ph1_x .* reshape(dDTx_ph1_L,N,1) + reshape(obj.U_T{1,1+f}.x(1:Nx,:,:),N,1) .* dS_D_ph2_x .* reshape(dDTx_ph2_L,N,1);
                        vecX2 = reshape(obj.U_T{1,1+f}.x(2:Nx+1,:,:),N,1).* dS_D_ph1_x .* reshape(dDTx_ph1_R,N,1) + reshape(obj.U_T{1,1+f}.x(2:Nx+1,:,:),N,1) .* dS_D_ph2_x .* reshape(dDTx_ph2_R,N,1);
                        vecY1 = reshape(obj.U_T{1,1+f}.y(:,1:Ny,:),N,1).* dS_D_ph1_y .* reshape(dDTy_ph1_L,N,1) + reshape(obj.U_T{1,1+f}.y(:,1:Ny,:),N,1) .* dS_D_ph2_y .* reshape(dDTy_ph2_L,N,1);
                        vecY2 = reshape(obj.U_T{1,1+f}.y(:,2:Ny+1,:),N,1).* dS_D_ph1_y .* reshape(dDTy_ph1_R,N,1) + reshape(obj.U_T{1,1+f}.y(:,2:Ny+1,:),N,1) .* dS_D_ph2_y .* reshape(dDTy_ph2_R,N,1);
                        vecZ1 = reshape(obj.U_T{1,1+f}.z(:,:,1:Nz),N,1) .* dS_D_ph1_z .* reshape(dDTz_ph1_L,N,1) + reshape(obj.U_T{1,1+f}.z(:,:,1:Nz),N,1) .* dS_D_ph2_z .* reshape(dDTz_ph2_L,N,1);
                        vecZ2 = reshape(obj.U_T{1,1+f}.z(:,:,2:Nz+1),N,1) .* dS_D_ph1_z .* reshape(dDTz_ph1_R,N,1) + reshape(obj.U_T{1,1+f}.z(:,:,2:Nz+1),N,1) .* dS_D_ph2_z .* reshape(dDTz_ph2_R,N,1);

                        DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2 + vecY2 + vecX2 - vecY1 - vecX1 - vecZ1, vecX1, vecY1, vecZ1];
                        DiagIndx = [-Grid.Nx*Grid.Ny, -Grid.Nx, -1, 0, 1, Grid.Nx, Grid.Nx*Grid.Ny];
                        Jz(:, 1+N*(j-1):N*j) = Jz(:, 1+N*(j-1):N*j) + spdiags(DiagVecs, DiagIndx, Grid.N, Grid.N);

                        % Capillarity
                        % Jz{i,j} = Jz{i,j} - obj.Tph{i,1} * spdiags(obj.dPc, 0, N, N);
                        Jz(:, 1+N*(j-1):N*j) = Jz(:, 1+N*(j-1):N*j) - obj.Thph{1,1+f} * spdiags(obj.dPc(Index.Start:Index.End), 0, N, N);
                    end

                    %% 3. J_EB_T Block
                    dMupxPh1 = reshape(obj.UpWind{1,1+f}.x * (obj.Mob(Index.Start:Index.End, 1) .* obj.h(Index.Start:Index.End,1) .* obj.drhodT(Index.Start:Index.End,1) + obj.Mob(Index.Start:Index.End, 1) .* obj.dhdT(Index.Start:Index.End,1) .* obj.rho(Index.Start:Index.End,1) + obj.dMobdT(Index.Start:Index.End, 1) .* obj.h(Index.Start:Index.End,1) .* obj.rho(Index.Start:Index.End,1)), Nx, Ny, Nz);
                    dMupyPh1 = reshape(obj.UpWind{1,1+f}.y * (obj.Mob(Index.Start:Index.End, 1) .* obj.h(Index.Start:Index.End,1) .* obj.drhodT(Index.Start:Index.End,1) + obj.Mob(Index.Start:Index.End, 1) .* obj.dhdT(Index.Start:Index.End,1) .* obj.rho(Index.Start:Index.End,1) + obj.dMobdT(Index.Start:Index.End, 1) .* obj.h(Index.Start:Index.End,1) .* obj.rho(Index.Start:Index.End,1)), Nx, Ny, Nz);
                    dMupzPh1 = reshape(obj.UpWind{1,1+f}.z * (obj.Mob(Index.Start:Index.End, 1) .* obj.h(Index.Start:Index.End,1) .* obj.drhodT(Index.Start:Index.End,1) + obj.Mob(Index.Start:Index.End, 1) .* obj.dhdT(Index.Start:Index.End,1) .* obj.rho(Index.Start:Index.End,1) + obj.dMobdT(Index.Start:Index.End, 1) .* obj.h(Index.Start:Index.End,1) .* obj.rho(Index.Start:Index.End,1)), Nx, Ny, Nz);
                    dMupxPh2 = reshape(obj.UpWind{2,1+f}.x * (obj.Mob(Index.Start:Index.End, 2) .* obj.h(Index.Start:Index.End,2) .* obj.drhodT(Index.Start:Index.End,2) + obj.Mob(Index.Start:Index.End, 2) .* obj.dhdT(Index.Start:Index.End,2) .* obj.rho(Index.Start:Index.End,2) + obj.dMobdT(Index.Start:Index.End, 2) .* obj.h(Index.Start:Index.End,2) .* obj.rho(Index.Start:Index.End,2)), Nx, Ny, Nz);
                    dMupyPh2 = reshape(obj.UpWind{2,1+f}.y * (obj.Mob(Index.Start:Index.End, 2) .* obj.h(Index.Start:Index.End,2) .* obj.drhodT(Index.Start:Index.End,2) + obj.Mob(Index.Start:Index.End, 2) .* obj.dhdT(Index.Start:Index.End,2) .* obj.rho(Index.Start:Index.End,2) + obj.dMobdT(Index.Start:Index.End, 2) .* obj.h(Index.Start:Index.End,2) .* obj.rho(Index.Start:Index.End,2)), Nx, Ny, Nz);
                    dMupzPh2 = reshape(obj.UpWind{2,1+f}.z * (obj.Mob(Index.Start:Index.End, 2) .* obj.h(Index.Start:Index.End,2) .* obj.drhodT(Index.Start:Index.End,2) + obj.Mob(Index.Start:Index.End, 2) .* obj.dhdT(Index.Start:Index.End,2) .* obj.rho(Index.Start:Index.End,2) + obj.dMobdT(Index.Start:Index.End, 2) .* obj.h(Index.Start:Index.End,2) .* obj.rho(Index.Start:Index.End,2)), Nx, Ny, Nz);

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

                    dacc_1 = obj.drhodT(Index.Start:Index.End,1).*obj.h(Index.Start:Index.End,1) + obj.rho(Index.Start:Index.End,1).*obj.dhdT(Index.Start:Index.End,1);
                    dacc_2 = obj.drhodT(Index.Start:Index.End,2).*obj.h(Index.Start:Index.End,2) + obj.rho(Index.Start:Index.End,2).*obj.dhdT(Index.Start:Index.End,2);

                    acc = pv /dt .* (obj.dSdT .* acc_1 + obj.S(Index.Start:Index.End,1) .* dacc_1 + obj.S(Index.Start:Index.End,2) .* dacc_2 - obj.dSdT .* acc_2) + Grid.Volume .* (1 - Medium.Por) / dt .* Medium.Rho .* Medium.Cpr;
                    DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2 + vecY2+vecX2-vecY1-vecX1-vecZ1+acc, vecX1, vecY1, vecZ1];
                    DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
                    JT = JT + spdiags(DiagVecs, DiagIndx, N, N) + obj.Tk{1,1+f};

                    dS_D_ph1_x = (obj.dSdT(Index.Start:Index.End,1) .* cond(:,1) .* Grid.Por) .* Grid.Ax / Grid.dx;
                    dS_D_ph1_y = (obj.dSdT(Index.Start:Index.End,1) .* cond(:,1) .* Grid.Por) .* Grid.Ay / Grid.dy;
                    dS_D_ph1_z = (obj.dSdT(Index.Start:Index.End,1) .* cond(:,1) .* Grid.Por) .* Grid.Az / Grid.dz;
                    dS_D_ph2_x = (-obj.dSdT(Index.Start:Index.End,1) .* cond(:,2) .* Grid.Por) .* Grid.Ax / Grid.dx;
                    dS_D_ph2_y = (-obj.dSdT(Index.Start:Index.End,1) .* cond(:,2) .* Grid.Por) .* Grid.Ay / Grid.dy;
                    dS_D_ph2_z = (-obj.dSdT(Index.Start:Index.End,1) .* cond(:,2) .* Grid.Por) .* Grid.Az / Grid.dz;

                    vecX1 = reshape(obj.U_T{1,1+f}.x(1:Nx,:,:),N,1).* dS_D_ph1_x .* reshape(dDTx_ph1_L,N,1) + reshape(obj.U_T{1,1+f}.x(1:Nx,:,:),N,1) .* dS_D_ph2_x .* reshape(dDTx_ph2_L,N,1);
                    vecX2 = reshape(obj.U_T{1,1+f}.x(2:Nx+1,:,:),N,1).* dS_D_ph1_x .* reshape(dDTx_ph1_R,N,1) + reshape(obj.U_T{1,1+f}.x(2:Nx+1,:,:),N,1) .* dS_D_ph2_x .* reshape(dDTx_ph2_R,N,1);
                    vecY1 = reshape(obj.U_T{1,1+f}.y(:,1:Ny,:),N,1).* dS_D_ph1_y .* reshape(dDTy_ph1_L,N,1) + reshape(obj.U_T{1,1+f}.y(:,1:Ny,:),N,1) .* dS_D_ph2_y .* reshape(dDTy_ph2_L,N,1);
                    vecY2 = reshape(obj.U_T{1,1+f}.y(:,2:Ny+1,:),N,1).* dS_D_ph1_y .* reshape(dDTy_ph1_R,N,1) + reshape(obj.U_T{1,1+f}.y(:,2:Ny+1,:),N,1) .* dS_D_ph2_y .* reshape(dDTy_ph2_R,N,1);
                    vecZ1 = reshape(obj.U_T{1,1+f}.z(:,:,1:Nz),N,1) .* dS_D_ph1_z .* reshape(dDTz_ph1_L,N,1) + reshape(obj.U_T{1,1+f}.z(:,:,1:Nz),N,1) .* dS_D_ph2_z .* reshape(dDTz_ph2_L,N,1);
                    vecZ2 = reshape(obj.U_T{1,1+f}.z(:,:,2:Nz+1),N,1) .* dS_D_ph1_z .* reshape(dDTz_ph1_R,N,1) + reshape(obj.U_T{1,1+f}.z(:,:,2:Nz+1),N,1) .* dS_D_ph2_z .* reshape(dDTz_ph2_R,N,1);

                    DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2 + vecY2 + vecX2 - vecY1 - vecX1 - vecZ1, vecX1, vecY1, vecZ1];
                    DiagIndx = [-Grid.Nx*Grid.Ny, -Grid.Nx, -1, 0, 1, Grid.Nx, Grid.Nx*Grid.Ny];
                    JT = JT + spdiags(DiagVecs, DiagIndx, Grid.N, Grid.N);

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
                        dMupPh1 = obj.UpWind{1,1+f} * ( obj.dMobdz(Index.Start:Index.End, 1, j) .* obj.x(Index.Start:Index.End,(c-1)*2+1) .* obj.rho(Index.Start:Index.End,1) + obj.Mob(Index.Start:Index.End, 1) .* obj.dxdz(Index.Start:Index.End,(c-1)*2+1, j) .* obj.rho(Index.Start:Index.End,1) + obj.Mob(Index.Start:Index.End, 1) .* obj.x(Index.Start:Index.End, (c-1)*2+1) .* obj.drhodz(Index.Start:Index.End, 1));
                        dMupPh2 = obj.UpWind{2,1+f} * ( obj.dMobdz(Index.Start:Index.End, 2, j) .* obj.x(Index.Start:Index.End,(c-1)*2+2) .* obj.rho(Index.Start:Index.End,2) + obj.Mob(Index.Start:Index.End, 2) .* obj.dxdz(Index.Start:Index.End,(c-1)*2+2, j) .* obj.rho(Index.Start:Index.End,2) + obj.Mob(Index.Start:Index.End, 2) .* obj.x(Index.Start:Index.End, (c-1)*2+2) .* obj.drhodz(Index.Start:Index.End, 2));

                        FluxDerivativePh1 = dMupPh1 .* abs(obj.U{1,1+f});
                        FluxDerivativePh2 = dMupPh2 .* abs(obj.U{2,1+f});

                        if c == obj.NofComponents
                            acc = - pv/dt .* rhoT + pv/dt .* obj.z(Index.Start:Index.End,c) .* obj.drhoTdz(Index.Start:Index.End, j);
                        elseif c == j
                            acc = pv/dt .* rhoT + pv/dt .* obj.z(Index.Start:Index.End,c) .* obj.drhoTdz(Index.Start:Index.End, j);
                        end
                        Jz = sparse(1:N, 1:N, acc)...
                            + (UpwindPermutation1Ph1 * sparse(1:N_Face, 1:N_Face, FluxDerivativePh1) * UpwindPermutation2Ph1')'...
                            + (UpwindPermutation1Ph2 * sparse(1:N_Face, 1:N_Face, FluxDerivativePh2) * UpwindPermutation2Ph2')';
                        % Capillarity
                        Jz = Jz - obj.Tph{2*c-1,1+f} * sparse(1:N, 1:N, obj.dPc(Index.Start:Index.End));
                    end
                otherwise
                    error("Discretization not found.")
            end
        end

        function Jacobian = BuildFullJacobian(obj, ProductionSystem, DiscretizationModel, Time, dt)
            Jhorizontal = cell(obj.NofComponents + 1, 1);
            Nt = DiscretizationModel.N;
            Nm = DiscretizationModel.ReservoirGrid.N;

            for c = 1:obj.NofComponents
                Jph{c} = sparse(Nt, obj.NofComponents*Nt);
                % Reservoir
                Index.Start = 1;
                Index.End = Nm;
                [Jp_res, Jz_res, JT_res] = BuildMassMediumJacobian(obj, ProductionSystem.Reservoir, DiscretizationModel.ReservoirGrid, dt, Index, 0, c);

                % Wells
                for i_well = 1:ProductionSystem.Wells.NofProd
                    Prod = ProductionSystem.Wells.Prod{i_well};
                    % Check if the well is open/close
                    isOpen = mod(sum(Time >= Prod.Schedule), 2);
                    if isOpen
                        [Jp_res, Jz_res, JT_res] = obj.AddWellsToMassJacobian(Jp_res, Jz_res, JT_res, obj.P(:,obj.NofPhases), obj.rho, obj.x, Prod, ProductionSystem.Reservoir.K, c);
                    end
                end

                % Fractures
                Jp = [];  Jz = []; JT = [];
                for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                    Nf = DiscretizationModel.FracturesGrid.N;
                    Index.Start = DiscretizationModel.Index_Local_to_Global( Nm, f, 1     );
                    Index.End   = DiscretizationModel.Index_Local_to_Global( Nm, f, Nf(f) );
                    [Jp_frac, Jz_frac, JT_frac] = BuildMassMediumJacobian(obj, ProductionSystem.FracturesNetwork.Fractures(f), DiscretizationModel.FracturesGrid.Grids(f), dt, Index, f, c);
                    Jp  = blkdiag(Jp, Jp_frac);
                    Jz  = blkdiag(Jz, Jz_frac);
                    JT  = blkdiag(JT, JT_frac);
                end

                Jp  = blkdiag(Jp_res,Jp);
                Jz  = blkdiag(Jz_res,Jz);
                JT  = blkdiag(JT_res,JT);

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
                    dDeltaP1_dTm = - 0.5 * obj.drhodT(idx_m,1) .* obj.GravityModel.g .* depthDiff;
                    dDeltaP1_dTf = - 0.5 * obj.drhodT(idx_f,1) .* obj.GravityModel.g .* depthDiff;
                    dDeltaP1_dzm = - 0.5 * obj.drhodz(idx_m,1) .* obj.GravityModel.g .* depthDiff;
                    dDeltaP1_dzf = - 0.5 * obj.drhodz(idx_f,1) .* obj.GravityModel.g .* depthDiff;

                    % Phase 2
                    DeltaP2 = (obj.P(idx_m,2) - obj.P(idx_f,2)) - rho_avg_2 .* obj.GravityModel.g .* depthDiff;
                    dDeltaP2_dPm = 1 - 0.5 * obj.drhodp(idx_m,2) .* obj.GravityModel.g .* depthDiff;
                    dDeltaP2_dPf = 1 - 0.5 * obj.drhodp(idx_f,2) .* obj.GravityModel.g .* depthDiff;
                    dDeltaP2_dTm = - 0.5 * obj.drhodT(idx_m,2) .* obj.GravityModel.g .* depthDiff;
                    dDeltaP2_dTf = - 0.5 * obj.drhodT(idx_f,2) .* obj.GravityModel.g .* depthDiff;
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

                    % 2. Temperature part
                    JT_conn = - T_Geo .* ( ...
                        UpWind(:,1) .* ( ...
                        obj.Mob(idx_m,1)    .* obj.rho(idx_m,1).*obj.x(idx_m,(c-1)*2+1) .* dDeltaP1_dTm + ...
                        obj.Mob(idx_m,1)    .* obj.rho(idx_m,1).*obj.dxdT(idx_m,(c-1)*2+1).*DeltaP1 + ...
                        obj.Mob(idx_m,1)    .* obj.drhodT(idx_m,1).*obj.x(idx_m,(c-1)*2+1).*DeltaP1  + ...
                        obj.dMobdT(idx_m,1) .* obj.rho(idx_m,1).*obj.x(idx_m,(c-1)*2+1).*DeltaP1)  + ...
                        ~UpWind(:,1) .*  ( ...
                        obj.Mob(idx_f,1) .*obj.rho(idx_f,1).*obj.x(idx_f,(c-1)*2+1) .* dDeltaP1_dTf + ...
                        obj.Mob(idx_f,1) .*obj.rho(idx_f,1).*obj.dxdT(idx_f,(c-1)*2+1).*DeltaP1 + ...
                        obj.Mob(idx_f,1) .*obj.drhodT(idx_f,1).*obj.x(idx_f,(c-1)*2+1).*DeltaP1 + ...
                        obj.dMobdT(idx_f,1) .*obj.rho(idx_f,1).*obj.x(idx_f,(c-1)*2+1).*DeltaP1) + ...
                        UpWind(:,2) .* ( ...
                        obj.Mob(idx_m,2) .* obj.rho(idx_m,2).*obj.x(idx_m,(c-1)*2+2) .* dDeltaP2_dTm + ...
                        obj.Mob(idx_m,2) .* obj.rho(idx_m,2).*obj.dxdT(idx_m,(c-1)*2+2).*DeltaP2 + ...
                        obj.Mob(idx_m,2) .* obj.drhodT(idx_m,2).*obj.x(idx_m,(c-1)*2+2).*DeltaP2  + ...
                        obj.dMobdT(idx_m,2) .* obj.rho(idx_m,2).*obj.x(idx_m,(c-1)*2+2).*DeltaP2)  + ...
                        ~UpWind(:,2) .*  ( ...
                        obj.Mob(idx_f,2) .*obj.rho(idx_f,2).*obj.x(idx_f,(c-1)*2+2) .* dDeltaP2_dTf + ...
                        obj.Mob(idx_f,2) .*obj.rho(idx_f,2).*obj.dxdT(idx_f,(c-1)*2+2).*DeltaP2 + ...
                        obj.Mob(idx_f,2) .*obj.drhodT(idx_f,2).*obj.x(idx_f,(c-1)*2+2).*DeltaP2 + ...
                        obj.dMobdT(idx_f,2) .*obj.rho(idx_f,2).*obj.x(idx_f,(c-1)*2+2).*DeltaP2 ) );

                    % 3. Composition part
                    Jz1_conn = T_Geo .* ( ...
                        UpWind(:,1) .* (obj.drhodz(idx_m,1) .* obj.Mob(idx_m,1) .* obj.x(idx_m,(c-1)*2+1) .* (-DeltaP1) + obj.rho(idx_m,1) .* obj.dMobdz(idx_m,1) .* obj.x(idx_m,(c-1)*2+1) .* (-DeltaP1) + obj.rho(idx_m,1) .* obj.Mob(idx_m,1) .* obj.dxdz(idx_m,(c-1)*2+1) .* (-DeltaP1) + obj.rho(idx_m,1) .* obj.Mob(idx_m,1) .* obj.x(idx_m,(c-1)*2+1) .* (-dDeltaP1_dzm))  + ...
                        UpWind(:,2) .* (obj.drhodz(idx_m,2) .* obj.Mob(idx_m,2) .* obj.x(idx_m,(c-1)*2+2) .* (-DeltaP2) + obj.rho(idx_m,2) .* obj.dMobdz(idx_m,2) .* obj.x(idx_m,(c-1)*2+2) .* (-DeltaP2) + obj.rho(idx_m,2) .* obj.Mob(idx_m,2) .* obj.dxdz(idx_m,(c-1)*2+2) .* (-DeltaP2) + obj.rho(idx_m,2) .* obj.Mob(idx_m,2) .* obj.x(idx_m,(c-1)*2+2) .* (-dDeltaP2_dzm)) );

                    Jz2_conn = T_Geo .* ( ...
                        ~UpWind(:,1) .* (obj.drhodz(idx_f,1) .* obj.Mob(idx_f,1).* obj.x(idx_f,(c-1)*2+1) .* (-DeltaP1) + obj.rho(idx_f,1) .* obj.dMobdz(idx_f,1).* obj.x(idx_f,(c-1)*2+1) .*(-DeltaP1) + obj.rho(idx_f,1) .* obj.Mob(idx_f,1) .* obj.dxdz(idx_f,(c-1)*2+1) .* (-DeltaP1) + obj.rho(idx_f,1) .* obj.Mob(idx_f,1) .* obj.x(idx_f,(c-1)*2+1) .* (-dDeltaP1_dzf)) + ...
                        ~UpWind(:,2) .* (obj.drhodz(idx_f,2) .* obj.Mob(idx_f,2).* obj.x(idx_f,(c-1)*2+2) .* (-DeltaP2) + obj.rho(idx_f,2) .* obj.dMobdz(idx_f,2).* obj.x(idx_f,(c-1)*2+2) .*(-DeltaP2) + obj.rho(idx_f,2) .* obj.Mob(idx_f,2) .* obj.dxdz(idx_f,(c-1)*2+2) .* (-DeltaP2) + obj.rho(idx_f,2) .* obj.Mob(idx_f,2) .* obj.x(idx_f,(c-1)*2+2) .* (-dDeltaP2_dzf)));

                    Jp(idx_f, idx_m) = Jp_conn;
                    JT(idx_f, idx_m) = JT_conn;
                    Jp(idx_f, idx_f) = Jp(idx_f, idx_f) - sum(Jp_conn);
                    JT(idx_f, idx_f) = JT(idx_f, idx_f) - sum(JT_conn);
                    % mat-frac or frac2 - frac1
                    Jp(idx_m, idx_f) = Jp_conn';
                    JT(idx_m, idx_f) = JT_conn';
                    Jp(sub2ind([Nt, Nt], idx_m, idx_m)) = Jp(sub2ind([Nt, Nt], idx_m, idx_m)) - Jp_conn;
                    JT(sub2ind([Nt, Nt], idx_m, idx_m)) = JT(sub2ind([Nt, Nt], idx_m, idx_m)) - JT_conn;

                    Jz(idx_f, idx_m) = Jz1_conn;
                    Jz(idx_f, idx_f) = Jz(idx_f, idx_f) + sum(Jz2_conn);
                    % mat-frac or frac2 - frac1
                    Jz(idx_m, idx_f) = -Jz2_conn;
                    % diag of mat or frac2
                    Jz(sub2ind([Nt, Nt], idx_m, idx_m)) = Jz(sub2ind([Nt, Nt], idx_m, idx_m)) - Jz1_conn;
                end

                % Combine the Jacobian blocks
                Jhorizontal{c} = horzcat(Jp, Jz, JT);
            end

            [J_EB_P, J_EB_Z, J_EB_T] = BuildEnergyMediumJacobian(obj, ProductionSystem.Reservoir, DiscretizationModel.ReservoirGrid, dt, Index, 0);
            
            % Wells
            for i_well = 1:ProductionSystem.Wells.NofProd
                Prod = ProductionSystem.Wells.Prod{i_well};
                % Check if the well is open/close
                isOpen = mod(sum(Time >= Prod.Schedule), 2);
                if isOpen
                    [J_EB_P, J_EB_Z, J_EB_T] = obj.AddWellsToEnergyJacobian(J_EB_P, J_EB_Z, J_EB_T, obj.P(:,obj.NofPhases), obj.rho, obj.h, Prod, ProductionSystem.Reservoir.K);
                end
            end
            
            Jp = [];  Jz = []; JT = [];         
            Jp  = blkdiag(J_EB_P, Jp);
            Jz  = blkdiag(J_EB_Z, Jz);
            JT  = blkdiag(J_EB_T, JT);

            Jhorizontal{obj.NofComponents+1} = horzcat(Jp, Jz, JT);
            Jacobian = Jhorizontal{1};
            %Jacobian_noconn = Jphnoconn{1};
            for c=2:obj.NofComponents+1
                Jacobian = vertcat(Jacobian, Jhorizontal{c});
            end
              
        end

        %%
        function [Jp, Jz, JT] = AddWellsToMassJacobian(obj, Jp, Jz, JT, p, rho, x, Prod, K, c)
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
                        ( obj.dMobdz(b(j), 1, zI) * rho(b(j), 1) * x(b(j), (c-1)*2+1) ...
                        + obj.dMobdz(b(j), 2, zI) * rho(b(j), 2) * x(b(j), (c-1)*2+2)...
                        + obj.Mob(b(j), 1) * rho(b(j), 1) * obj.dxdz(b(j), (c-1)*2+1, zI) ...
                        + obj.Mob(b(j), 2) * rho(b(j), 2) * obj.dxdz(b(j), (c-1)*2+2, zI)...
                        + obj.Mob(b(j), 1) * obj.drhodz(b(j), 1, zI) * x(b(j), (c-1)*2+1) ...
                        + obj.Mob(b(j), 2) * obj.drhodz(b(j), 2, zI) * x(b(j), (c-1)*2+2)...
                        ) * (Prod.p(j) - p(b(j)));
                end
                JT(b(j),b(j)) = JT(b(j),b(j)) ...
                            - Prod.PI(j) * K(b(j)) * ...
                            ( obj.dMobdT(b(j), 1) * rho(b(j), 1) * x(b(j), (c-1)*2+1) ...
                            + obj.dMobdT(b(j), 2) * rho(b(j), 2) * x(b(j), (c-1)*2+2)...
                            + obj.Mob(b(j), 1) * obj.drhodT(b(j), 1) * x(b(j), (c-1)*2+1) ...
                            + obj.Mob(b(j), 2) * obj.drhodT(b(j), 2) * x(b(j), (c-1)*2+2) ...
                            + obj.Mob(b(j), 1) * rho(b(j), 1) * obj.dxdT(b(j), (c-1)*2+1) ...
                            + obj.Mob(b(j), 2) * rho(b(j), 2) * obj.dxdT(b(j), (c-1)*2+2)) * (Prod.p(j) - p(b(j)));
            end
        end

        function [Jp, Jz, JT] = AddWellsToEnergyJacobian(obj, Jp, Jz, JT, p, rho, h, Prod, K)
            N = length(p);
            b = Prod.Cells;
            for j=1:length(b)
                %Pressure blocks
                Jp(b(j),b(j)) = Jp(b(j),b(j))...
                    + Prod.PI(j) * K(b(j)) *...
                    ( obj.Mob(b(j), 1) * rho(b(j),1) * h(b(j), 1) ...
                    + obj.Mob(b(j), 2) * rho(b(j),2) * h(b(j), 2)...
                    )...
                    - Prod.PI(j) * K(b(j)) * obj.Mob(b(j), 1) * h(b(j), 1) * obj.drhodp(b(j),1) * (Prod.p(j) - p(b(j))) ...
                    - Prod.PI(j) * K(b(j)) * obj.Mob(b(j), 2) * h(b(j), 2) * obj.drhodp(b(j),2) * (Prod.p(j) - p(b(j))) ...
                    - Prod.PI(j) * K(b(j)) * obj.Mob(b(j), 1) * obj.dhdp(b(j), 1) * rho(b(j),1) * (Prod.p(j) - p(b(j))) ...
                    - Prod.PI(j) * K(b(j)) * obj.Mob(b(j), 2) * obj.dhdp(b(j), 2) * rho(b(j),2) * (Prod.p(j) - p(b(j))) ...
                    - Prod.PI(j) * K(b(j)) * obj.dMobdp(b(j), 1) * h(b(j), 1) * rho(b(j),1) * (Prod.p(j) - p(b(j))) ...
                    - Prod.PI(j) * K(b(j)) * obj.dMobdp(b(j), 2) * h(b(j), 2) * rho(b(j),2) * (Prod.p(j) - p(b(j)));
                for zI = 1:obj.NofComponents-1
                    Jz(b(j),(zI-1)*N+b(j)) = Jz(b(j),(zI-1)*N+b(j))...
                        - Prod.PI(j) * K(b(j)) *...
                        ( obj.dMobdz(b(j), 1, zI) * rho(b(j), 1) * h(b(j), 1) ...
                        + obj.dMobdz(b(j), 2, zI) * rho(b(j), 2) * h(b(j), 2)...
                        + obj.Mob(b(j), 1) * rho(b(j), 1) * obj.dhdz(b(j), 1, zI) ...
                        + obj.Mob(b(j), 2) * rho(b(j), 2) * obj.dhdz(b(j), 2, zI)...
                        + obj.Mob(b(j), 1) * obj.drhodz(b(j), 1, zI) * h(b(j), 1) ...
                        + obj.Mob(b(j), 2) * obj.drhodz(b(j), 2, zI) * h(b(j), 2)...
                        ) * (Prod.p(j) - p(b(j)));
                end
                JT(b(j),b(j)) = JT(b(j),b(j)) ...
                            - Prod.PI(j) * K(b(j)) * ...
                            ( obj.dMobdT(b(j), 1) * rho(b(j), 1) * h(b(j), 1) ...
                            + obj.dMobdT(b(j), 2) * rho(b(j), 2) * h(b(j), 2)...
                            + obj.Mob(b(j), 1) * obj.drhodT(b(j), 1) * h(b(j), 1) ...
                            + obj.Mob(b(j), 2) * obj.drhodT(b(j), 2) * h(b(j), 2) ...
                            + obj.Mob(b(j), 1) * rho(b(j), 1) * obj.dhdT(b(j), 1) ...
                            + obj.Mob(b(j), 2) * rho(b(j), 2) * obj.dhdT(b(j), 2)) * (Prod.p(j) - p(b(j)));
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
                deltaz = delta(Nt+1: end-Nt);
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
                
                % 1.c Update temperature T
                deltaT = delta((obj.NofComponents * Nt + 1):(obj.NofComponents * Nt + Nt));
                Tm = ProductionSystem.Reservoir.State.Properties('T');
                Tm.update(deltaT(1:Nm));

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
                            %IZ = Nt*(c-1) + IP;
                            %EZ = Nt*(c-1) + EP;
                            Zf = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['z_', num2str(c)]);
                            % Zf.update(delta(c*Nf(f) + 1:(c+1)*Nf(f)));
                            Deltaz = delta(c*Nt + IP : c*Nt + EP); 

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

                        % 2.c Update T
                        Tf = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('T');
                        Tf.update(deltaT(IP:EP));

                        % Update Remaining properties (x, S, rho, rhoT, Pc)
                        obj.UpdatePhaseCompositions(ProductionSystem.FracturesNetwork.Fractures(f), DiscretizationModel, FluidModel, 0);
                    end
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

    end
end