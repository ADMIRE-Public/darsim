% Compositional fluid model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini and Barnaby Fryer
%TU Delft
%Created: 14 July 2016
%Last modified: 7 March 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Compositional_fluid_model < fluid_model
    properties
        FlashCalculator
    end
    methods
        function obj = Compositional_fluid_model(n_phases, n_comp)
            obj@fluid_model(n_phases, n_comp);
            obj.name = 'Compositional';
        end
        %% Initialize Injectors
        function InitializeInjectors(obj, Inj)  
            nc = length(obj.Components);
            nph = length(obj.Phases);
            % Loop over all injectors
            for i=1:length(Inj)
                N = length(Inj{i}.Cells);
                c = Inj{i}.Composition;
                Inj{i}.z = nan(N, nc);
                Inj{i}.x = zeros(N, nph * nc); 
                Inj{i}.x(:, nph*(c-1)+1:nph*c) = 1;
                for iph = 1:nph
                    Inj{i}.rho(:, iph) = obj.Phases(iph).ComputeDensity(Inj{i}.p, obj.Components, Inj{i}.x);
                end
                Inj{i}.S = 1;
                % Inj{i}.Mob = obj.ComputePhaseMobilities(Inj{i}.S);
                if c < nc
                    mu = obj.Phases(1).mu;
                    if isempty(mu) || all(isnan(mu(:)))
                        mu = obj.Phases(1).ComputeViscosity(Inj{i}.p, Inj{i}.T, Inj{i}.x, obj.Thermo);
                    end
                    Inj{i}.Mob = [1/mu 0];
                else
                    mu = obj.Phases(2).mu;
                    if isempty(mu) || all(isnan(mu(:)))
                        mu = obj.Phases(2).ComputeViscosity(Inj{i}.p, Inj{i}.T, Inj{i}.x, obj.Thermo);
                    end
                    Inj{i}.Mob = [0 1/mu];
                end
            end        
        end
        %% Call Flash
        function SinglePhase = Flash(obj, Status)
           SinglePhase = obj.FlashCalculator.Flash(Status, obj.Components, obj.Phases, obj.Thermo);
        end
        %% Phase density and Phase saturation
        function ComputePhaseSaturation(obj, Status, SinglePhase)
            N = length(SinglePhase);
            rho = zeros(N, obj.NofPhases);
            ni = zeros(N, obj.NofPhases);
            % Store in local variables
            NiRhoT = zeros(N, 1);
            for i=1:obj.NofPhases
                rho(:,i) = Status.Properties(['rho_',num2str(i)]).Value;
                ni(:,i) = Status.Properties(['ni_',num2str(i)]).Value;
                NiRhoT = NiRhoT + ni(:,i)./rho(:,i);
            end
            % Compute saturations
            for i=1:obj.NofPhases
                S = Status.Properties(['S_',num2str(i)]);
                S.Value = (ni(:,i)./rho(:,i)) ./ NiRhoT;
                % S.Value = max(S.Value, 0.0001); % safe-guard
            end
        end
        function ComputePhaseDensities(obj, Status, SinglePhase)
            % [obj.Rs, obj.dRs] = obj.FlashCalculator.KvaluesCalculator.ComputeRs(Status, obj.Phases);

            for i=1:obj.NofPhases
                % Correct Rs and dRs for undersaturated cells
                % [obj.Rs(:,i), obj.dRs(:,i)] = obj.Phases(i).RsOfUnderSaturatedPhase(Status.Properties('z_1').Value, obj.Components, obj.Rs(:,i), obj.dRs(:,i), SinglePhase);
                rho = Status.Properties(['rho_', num2str(i)]);
                % if i == 1
                %     if all(Status.Properties('P_1').Value == 0)
                %         rho.Value = obj.Phases(i).ComputeDensity(Status.Properties('P_2').Value, obj.Components, Status);
                %     else
                %         rho.Value = obj.Phases(i).ComputeDensity(Status.Properties('P_1').Value, obj.Components, Status);
                %     end
                % else
                %     rho.Value = obj.Phases(i).ComputeDensity(Status.Properties('P_2').Value, obj.Components, Status.Properties('x_1ph2').Value);
                % end
                rho.Value = obj.Phases(i).ComputeDensity(Status.Properties('P_2').Value, obj.Components, Status);
            end
        end
        function z = ComputeTotalFractions(obj, Status, N)
            %Two phase, two component total mole fraction updater
            %Based on mass balance equation z_1 * rho_t = x11*rho1*s1 + x12*rho2*s2
            S = zeros(N, obj.NofPhases);
            rho = zeros(N, obj.NofPhases);
            x = zeros(N, obj.NofComp*obj.NofPhases);
            
            for ph=1:obj.NofPhases
                rho(:, ph) = Status.Properties(['rho_', num2str(ph)]).Value;
                S(:, ph) = Status.Properties(['S_', num2str(ph)]).Value;
                for c=1:obj.NofComp
                    x(:,(c-1)*obj.NofPhases + ph) = Status.Properties(['x_', num2str(c),'ph',num2str(ph)]).Value;
                end
            end
            rhoT = Status.Properties('rhoT').Value;
            for i=1:obj.NofComp
                z = Status.Properties(['z_', num2str(i)]);
                Num = zeros(N,1);
                for j=1:obj.NofPhases
                    Num = Num + x(:,(i-1)*obj.NofPhases + j) .* S(:, j) .* rho(:,j); 
                end
                z.Value = Num ./ rhoT; 
            end
        end
        function SinglePhase = CheckNumberOfPhases(obj, SinglePhase, PreviousSinglePhase, Status, k)
            N = length(SinglePhase);
            z = zeros(N, obj.NofComp);
            for i=1:obj.NofComp
                z(:, i) = Status.Properties(['z_', num2str(i)]).Value;
                for j=1:obj.NofPhases
                    x = Status.Properties(['x_', num2str(i),'ph',num2str(j)]);
                    x.Value(z(:,i) == 1) = i==j;
                    x.Value(z(:,i) == 0) = i==j;
                end
            end
            
            % Check if single component
            SinglePhase(z(:,1) == 1) = 1;
            SinglePhase(z(:,1) == 0) = 2;
            
            % Transform mass fractions into mole fractions to check phase
            % state
            BubCheck = zeros(length(z), 2);
            BubCheck(PreviousSinglePhase == 2, :) = z(PreviousSinglePhase == 2,:) .* k (PreviousSinglePhase == 2, :);
            BubCheck = sum(BubCheck, 2);
            SinglePhase(BubCheck > 1) = 0;
            
            % 2.b: checking if it 's all vapor: checks if mix is above dew
            % point
            DewCheck = zeros(length(z), 2);
            DewCheck(PreviousSinglePhase == 1,:) = z(PreviousSinglePhase == 1,:) ./ k(PreviousSinglePhase == 1,:);
            DewCheck = sum(DewCheck, 2);
            SinglePhase(DewCheck(PreviousSinglePhase == 1) > 1) = 0;
        end
        %% K-value
        function k = ComputeKvalues(obj, Status)
            k = obj.FlashCalculator.KvaluesCalculator.Compute(Status, obj.Thermo);
        end
        function dkdp = DKvalDp(obj, Status)
            dkdp = obj.FlashCalculator.KvaluesCalculator.DKvalDp(Status, obj.Thermo);
        end
        %% d(x)/d(p)   &&   d(x)/d(z)
        function [dxdp, dxdz] = DxDpDz_backup(obj, Status, SinglePhase)
            % x = x1v, x1l, x2v, x2l, ni
            % z = p, z1
            N = length(SinglePhase);
            Ni = zeros(N, obj.NofPhases);
            x = zeros(N, obj.NofPhases*obj.NofComponents);
            for j=1:obj.NofPhases
                Ni(:, j) = Status.Properties(['ni_', num2str(j)]).Value;
                for i=1:obj.NofComponents
                    x(:,(i-1)*obj.NofPhases + j) = Status.Properties(['x_', num2str(i),'ph',num2str(j)]).Value;
                end
            end
            ni = Ni(:,1);
            
            k = obj.ComputeKvalues(Status);
            dk = obj.DKvalDp(Status);
            % Loop over all cells and do local inversion
            dxdp = zeros(N, 2*obj.NofComponents+1);
            dxdz = zeros(N, 2*obj.NofComponents+1, obj.NofComponents-1);
            for i = 1:N
                dFdx = zeros(2*obj.NofComponents + 1);
                dFdz = zeros(2*obj.NofComponents + 1, obj.NofComponents);
                for c=1:obj.NofComponents
                    % x_cv - k_c x_cl = 0
                    dFdx (c, (c-1)*2+1:(c-1)*2+2) = [1, -k(i,c)];
                    dFdz (c, 1) = -dk(i,c) * x(i, (c-1)*2+2);
                    % z_c - ni*x_cv - (1-ni)*x_cl = 0
                    dFdx (obj.NofComponents + c, (c-1)*2+1:(c-1)*2+2) = [-ni(i), ni(i)-1];
                    dFdx (obj.NofComponents + c, end) = x(i, (c-1)*2+2) - x(i, (c-1)*2+1);
                    if c<obj.NofComponents
                        dFdz (obj.NofComponents + c, c+1) = 1;
                    end
                    % Sum x_cv-x_cl = 0
                    dFdx (2*obj.NofComponents + 1, (c-1)*2+1:(c-1)*2+2) = [1, -1]; 
                end
                dFdz (2*obj.NofComponents, 2:obj.NofComponents) = -1;
                dxdp(i,:) = (dFdx\dFdz(:,1))';
                for j=1:obj.NofComponents-1
                    dxdz(i,:,j) = (dFdx\dFdz(:,j+1))';
                end
            end
            dxdp(SinglePhase == 1, :) = 0;
            dxdp(SinglePhase == 2, :) = 0;
            
            dxdz(SinglePhase == 1, 1) = 1; %1
            dxdz(SinglePhase == 1, 2) = 0;
            dxdz(SinglePhase == 1, 3) = -1; %-1
            dxdz(SinglePhase == 1, 4) = 0;
            dxdz(SinglePhase == 1, 5) = 0;
            
            % 
            dxdz(SinglePhase == 2, 1) = 0;
            dxdz(SinglePhase == 2, 2) = 1; %1
            dxdz(SinglePhase == 2, 3) = 0;
            dxdz(SinglePhase == 2, 4) = -1; %-1
            dxdz(SinglePhase == 2, 5) = 0;
            
        end
        function [dxdp, dxdz] = DxDpDz_backup2(obj, Status, SinglePhase)
            % Generated by Mengjie
            % x = x1v, x1l, x2v, x2l, ni
            % z = p, z1
            N = length(SinglePhase);
            Ni = zeros(N, obj.NofPhases);
            x = zeros(N, obj.NofPhases*obj.NofComponents);
            for j=1:obj.NofPhases
                Ni(:, j) = Status.Properties(['ni_', num2str(j)]).Value;
                for i=1:obj.NofComponents
                    x(:,(i-1)*obj.NofPhases + j) = Status.Properties(['x_', num2str(i),'ph',num2str(j)]).Value;
                end
            end
            z = Status.Properties('z_1').Value;
            ni = Ni(:,1);

            k = obj.ComputeKvalues(Status);
            
            dk = obj.DKvalDp(Status);
            dxdp = zeros(N, 2*obj.NofComponents+1);
            dxdz = zeros(N, 2*obj.NofComponents+1, obj.NofComponents-1);

            num_ni = (k(:,2)-1)+ z.* (k(:,1)-k(:,2));
            den_ni = (k(:,1)-1) .*(k(:,2)-1);
            dnum_ni = (1 - z) .* dk(:,2) + z .* dk(:,1);
            dden_di = (k(:,2)-1) .* dk(:,1) + (k(:,1)-1) .* dk(:,2);
            dxdp(:, 5) = - (dnum_ni .* den_ni - dden_di .* num_ni) ./ den_ni.^2;
            dxdp(:, 2) =  - ( z .* ( dxdp(:,5) .* (k(:,1) - 1) + ni .* dk(:,1) ) ) ./ ( ( ni .* (k(:,1) - 1) + 1 ).^2 );
            dxdp(:, 1) = k(:,1) .* dxdp(:,2) + x(:,2) .* dk(:,1);
            dxdp(:, 3) = - dxdp(:,1);
            dxdp(:, 4) = - dxdp(:,2);
            
            dxdz(:, 5) = (k(:,2) - k(:,1) )./ ((k(:,1)-1) .* (k(:,2)-1));
            dxdz(:, 2) = (ni .* (k(:,1) - 1) + 1 - z .* ( dxdz(:,5) .* (k(:,1) - 1) ) ) ./ ( ( ni .* (k(:,1) - 1) + 1 ).^2 );
            dxdz(:, 1) = k(:,1) .* dxdz(:,2);
            dxdz(:, 3) = - dxdz(:, 1);
            dxdz(:, 4) = - dxdz(:, 2);
                        
            dxdp(SinglePhase == 1, :) = 0;
            dxdp(SinglePhase == 2, :) = 0;
            
            dxdz(SinglePhase == 1, 1) = 1; %1
            dxdz(SinglePhase == 1, 2) = 0;
            dxdz(SinglePhase == 1, 3) = -1; %-1
            dxdz(SinglePhase == 1, 4) = 0;
            dxdz(SinglePhase == 1, 5) = 0;
            
            % 
            dxdz(SinglePhase == 2, 1) = 0;
            dxdz(SinglePhase == 2, 2) = 1; %1
            dxdz(SinglePhase == 2, 3) = 0;
            dxdz(SinglePhase == 2, 4) = -1; %-1
            dxdz(SinglePhase == 2, 5) = 0;
        end
        function [dxdp, dxdz] = DxDpDz(obj, Status, SinglePhase)
            % x = x1v, x1l, x2v, x2l, ..., ni
            % z = p, z1, z2, ...
            N = length(SinglePhase);
            nc = obj.NofComponents;
            Ni = zeros(N, obj.NofPhases);
            x = zeros(N, obj.NofPhases*obj.NofComponents);
            for j=1:obj.NofPhases
                Ni(:, j) = Status.Properties(['ni_', num2str(j)]).Value;
                for i=1:obj.NofComponents
                    x(:,(i-1)*obj.NofPhases + j) = Status.Properties(['x_', num2str(i),'ph',num2str(j)]).Value;
                end
            end
            p = Status.Properties('P_2').Value;
            z = Status.Properties('z_1').Value;
            V = Ni(:, 1);
            L = Ni(:, 2);
            
            switch nc
            case 2
                k = obj.ComputeKvalues(Status);
                dkdp = obj.DKvalDp(Status);
            case 3
                k = zeros(N, nc);
                dkdp = zeros(N, nc);
                dkdy1 = zeros(N, nc);

                Mw_components = [obj.Thermo.comp{1}.Mw, obj.Thermo.comp{2}.Mw, obj.Thermo.comp{3}.Mw];
                % As thermo. functions are a functions of f(p, y1_mole)
                y1_molar = x(:, 1) / Mw_components(1) ./ sum(x(:, 1:2:end) ./ Mw_components, 2);
                y1_molar(isnan(y1_molar)) = 0; % Single liquid phase
                
                for i = 1:nc
                    % Derivatives of K_mass-values
                    k(:, i) = obj.Thermo.comp{i}.K_mass(p, y1_molar);
                    dkdp(:, i) = obj.Thermo.comp{i}.dKdp(p, y1_molar);
                    dkdy1(:, i) = obj.Thermo.comp{i}.dKdy1(p, y1_molar);
                end
            end
            dxdp = zeros(N, 2*obj.NofComponents+1);
            dxdz = zeros(N, 2*obj.NofComponents+1, obj.NofComponents-1);

            TwoPhase_index = find(SinglePhase == 0);
            for i = TwoPhase_index'
                A = zeros(2*nc+2);
                b = zeros(2*nc+2, 1);
                for ic = 1:nc
                    A(ic, 2*ic-1) = 1;
                    A(ic, 2*ic  ) = -k(i, ic);
                    A(nc + ic, 2*ic-1) = -V(i);
                    A(nc + ic, 2*ic  ) = -L(i);
                    A(nc + ic, end-1 ) = -x(i, 2*ic-1);
                    A(nc + ic, end   ) = -x(i, 2*ic  );
                    b(ic) = dkdp(i, ic) .* x(i, 2*ic);
                end
                A(end-1, 1:2:end-3) = 1;
                A(end, 2:2:end-2) = 1;

                ANS = A\b;
                dxdp(i, :) = ANS(1:2*nc+1);

                if nc == 3
                    A(1:nc, 1) = A(1:nc, 1) - (dkdy1(i, :) .* x(i, 2:2:end))';
                end

                for ic = 1:obj.NofComponents-1
                    b = zeros(2*nc+2, 1);
                    b(nc + ic) = -1;
                    b(2*nc) = +1;

                    ANS = A\b;
                    dxdz(i, :, ic) = ANS(1:2*nc+1);
                end
            end
                        
            dxdp(SinglePhase == 1, :) = 0;
            dxdp(SinglePhase == 2, :) = 0;
            for ic = 1:nc-1
            dxdz(SinglePhase == 1, ic*2-1, ic) = 1;
            dxdz(SinglePhase == 1, ic*2  , ic) = 0;
            dxdz(SinglePhase == 1, nc*2-1, ic) = -1;
            dxdz(SinglePhase == 1, nc*2  , ic) = 0;
            dxdz(SinglePhase == 1, end   , ic) = 0;
            % 
            dxdz(SinglePhase == 2, ic*2-1, ic) = 0;
            dxdz(SinglePhase == 2, ic*2  , ic) = 1;
            dxdz(SinglePhase == 2, nc*2-1, ic) = 0;
            dxdz(SinglePhase == 2, nc*2  , ic) = -1;
            dxdz(SinglePhase == 2, end   , :) = 0;
            end
        end
        %% d(rho)/d(p)   &&   d(rho)/d(z)
        function drhodp = ComputeDrhoDp(obj, Status, SinglePhase, dxdp)
            drhodp = zeros(length(SinglePhase), obj.NofPhases);
            for i=1:obj.NofPhases
                drhodp(:, i) = obj.Phases(i).ComputeDrhoDp(Status.Properties('P_2').Value, obj.Components, dxdp, Status);
            end
        end
        function drhodz = ComputeDrhoDz(obj, Status, SinglePhase, dxdz)
            N = length(Status.Properties('P_2').Value);
            nph = obj.NofPhases;
            nc = obj.NofComponents;
            drhodz = zeros(N, nph, nc-1);
            % if nc > 2 % d(x)/d(z) = 0 in 2-component system
            %     for iph = 1:nph
            %         for ic = 1:nc-1
            %             drhodz(:, iph, ic) = obj.Phases(iph).ComputeDrhoDz(Status.Properties(['P_', num2str(iph)]).Value, obj.Components, dxdz, Status, ic);
            %         end
            %     end
            % end
            for iph = 1:nph
                for ic = 1:nc-1
                    drhodz(:, iph, ic) = obj.Phases(iph).ComputeDrhoDz(Status.Properties('P_2').Value, obj.Components, dxdz, Status, ic);
                end
            end
        end
        %% d(Saturation)/d(p)   &&   d(Saturation)/d(z)
        % Saturation of Vapor phase: S_V = (V/rho_V) / (V/rho_V + L/rho_L)
        %                                = (V*rho_L) / (V*rho_L + L*rho_V)
        function dSdp = ComputeDSDp(obj, Status, drhodp, dni)
            ni = Status.Properties('ni_1').Value;
            rhov = Status.Properties('rho_1').Value;
            rhol = Status.Properties('rho_2').Value;
            drhov = drhodp(:,1);
            drhol = drhodp(:,2);
            Num = rhol .* ni;
            Den = rhol .* ni + (1 - ni) .* rhov ;
            dNum = drhol .* ni + dni .* rhol;
            dDen = drhol .* ni + dni .* rhol + drhov - ni .* drhov - rhov .* dni;
            % final derivative
            dSdp = (dNum .* Den - dDen .* Num)./Den.^2;
            dSdp(Num == 0) = 0;
            % dSdp(Status.Properties('immobile_label').Value == 1) = 0;
        end
        function dSdz = ComputeDSDz(obj, Status, drhodz, dni)
            % Derivative of S with respect to overall composition
            rhov = Status.Properties('rho_1').Value;
            rhol = Status.Properties('rho_2').Value;
            ni = Status.Properties('ni_1').Value;
            dSdz = zeros(length(ni), obj.NofComponents - 1);
            for j = 1:obj.NofComponents-1
                Num1 = rhol .* ni;
                dNum1 = rhol .* dni(:, :, j) + drhodz(:, 2, j) .* ni;
                Den1 = rhol .* ni + (1 - ni) .* rhov;
                dDen1 = rhol .* dni(:, :, j) + drhodz(:, 2, j) .* ni + drhodz(:, 1, j) - dni(:, :, j) .* rhov - ni .* drhodz(:, 1, j);
                dSdz(:, j) =(Den1 .* dNum1 - Num1 .* dDen1) ./ Den1.^2;
            end
            dSdz(Status.Properties('S_1').Value == 0, :) = 0;
            dSdz(Status.Properties('S_1').Value == 1, :) = 0;
            dSdz(Status.Properties('immobile_label').Value == 1, :) = 0;
            % dSdz(Status.Properties('S_1').Value == 0, 2) = 0;
            % dSdz(Status.Properties('S_1').Value == 0, 1) = 0;
            % dSdz(Status.Properties('S_1').Value == 1, 2) = 0;
            % dSdz(Status.Properties('S_1').Value == 1, 1) = 0;
            % 
            % dSdz(Status.Properties('immobile_label').Value == 1, 2) = 0;
            % dSdz(Status.Properties('immobile_label').Value == 1, 1) = 0;
        end
        %% d(Mob)/d(p)   &&   d(Mob)/d(z)
        function dMobdp = ComputeDMobDp(obj, Medium, dSdp)
            % Compute d(kr)/d(p)
            Status = Medium.State;
            s = Status.Properties('S_1').Value;
            % % kr = obj.RelPermModel.kr;
            % kr = obj.RelPermModel.ComputeRelPerm(obj.Phases, Status);

            if isa(Medium, 'reservoir')
                if obj.HysteresisActive
                    kr = obj.RelPermModel.ComputeRelPermHysteresis(obj.Phases, Medium);
                    dkrdS = obj.RelPermModel.ComputeDerivativeHysteresis();
                else
                    kr = obj.RelPermModel.ComputeRelPerm(obj.Phases, Medium);
                    dkrdS = obj.RelPermModel.ComputeDerivative(obj.Phases, Medium);
                end
            elseif isa(Medium, 'fracture')
                kr = obj.RelPermModelFractures.ComputeRelPerm(obj.Phases, Medium, obj.RelPermModel);
                dkrdS = obj.RelPermModelFractures.ComputeDerivative(obj.Phases, Medium, obj.RelPermModel);
            else
                error('Medium is not recognized.')
            end

            dkrdp = dkrdS .* dSdp;
            dMobdp = zeros(length(kr), obj.NofPhases);
             % Compute d(mu)/d(p)
            for i=1:obj.NofPhases
                mu = obj.Phases(i).mu;
                if isempty(mu) || all(isnan(mu(:)))
                    mu = Status.Properties(['mu_', num2str(i)]).Value;
                    dmudp = obj.Phases(i).ComputeDmuDp(Status.Properties('P_2').Value, Status.T, Status, obj.Thermo);
                else
                    dmudp = 0;
                end
                 % Compute d(Mob)/d(p)
                dMobdp(:, i) = (dkrdp(:, i) .* mu - dmudp .* kr(:, i)) ./ mu.^2;
            end
        end
        % function dMobdp = ComputeDMobDpHysteresis(obj, Status, dSdp)
        %     % Compute d(kr)/d(p)
        %     kr = obj.RelPermModel.kr;
        %     dkrdS = obj.RelPermModel.ComputeDerivativeHysteresis();
        %     dkrdp = dkrdS .* dSdp;
        %     dMobdp = zeros(length(kr), obj.NofPhases);
        %      % Compute d(mu)/d(p)
        %     for i=1:obj.NofPhases
        %         mu = obj.Phases(i).mu;
        %         if isempty(mu) || all(isnan(mu(:)))
        %             mu = Status.Properties(['mu_', num2str(i)]).Value;
        %             dmudp = obj.Phases(i).ComputeDmuDp(Status.Properties(['P_', num2str(i)]).Value, Status.T, Status, obj.Thermo);
        %         else
        %             dmudp = 0;
        %         end
        %          % Compute d(Mob)/d(p)
        %         dMobdp(:, i) = (dkrdp(:, i) .* mu - dmudp .* kr(:, i)) ./ mu.^2;
        %     end
        % end
        % function dMobdp = ComputeDMobDp(obj, Status, dSdp)
        %     % Compute d(kr)/d(p)
        %     s = Status.Properties('S_1').Value;
        %     % kr = obj.RelPermModel.kr;
        %     kr = obj.RelPermModel.ComputeRelPerm(obj.Phases, Status);
        %     dkrdS = obj.RelPermModel.ComputeDerivative(obj.Phases, s);
        %     dkrdp = dkrdS .* dSdp;
        %     dMobdp = zeros(length(kr), obj.NofPhases);
        %      % Compute d(mu)/d(p)
        %     for i=1:obj.NofPhases
        %         mu = obj.Phases(i).mu;
        %         if isempty(mu) || all(isnan(mu(:)))
        %             mu = Status.Properties(['mu_', num2str(i)]).Value;
        %             dmudp = obj.Phases(i).ComputeDmuDp(Status.Properties(['P_', num2str(i)]).Value, Status.T, Status, obj.Thermo);
        %         else
        %             dmudp = 0;
        %         end
        %          % Compute d(Mob)/d(p)
        %         dMobdp(:, i) = (dkrdp(:, i) .* mu - dmudp .* kr(:, i)) ./ mu.^2;
        %     end
        % end
        function dMobdz = ComputeDMobDz(obj, Medium, dSdz, dxdz)
            Status = Medium.State;
            N = length(dSdz);
            dMobdz = zeros(N, obj.NofPhases, obj.NofComponents-1);
            % Compute d(kr)/d(z)
            dkrdz = zeros(N, obj.NofPhases, obj.NofComponents-1);
            s = Status.Properties('S_1').Value;

            if isa(Medium, 'reservoir')
                if obj.HysteresisActive
                    kr = obj.RelPermModel.ComputeRelPermHysteresis(obj.Phases, Medium);
                    dkrdS = obj.RelPermModel.ComputeDerivativeHysteresis();
                else
                    kr = obj.RelPermModel.ComputeRelPerm(obj.Phases, Medium);
                    dkrdS = obj.RelPermModel.ComputeDerivative(obj.Phases, Medium);
                end
            elseif isa(Medium, 'fracture')
                kr = obj.RelPermModelFractures.ComputeRelPerm(obj.Phases, Medium, obj.RelPermModel);
                dkrdS = obj.RelPermModelFractures.ComputeDerivative(obj.Phases, Medium, obj.RelPermModel);
            else
                error('Medium is not recognized.')
            end
            
            for i=1:obj.NofPhases
                for j=1:obj.NofComponents-1
                    dkrdz(:, i, j) = dkrdS(:, i) .* dSdz(:, j);
                end
            end
             % Compute d(mu)/d(z)
            for i=1:obj.NofPhases
                mu = obj.Phases(i).mu;
                dmudz = zeros(N, obj.NofComponents);
                if isempty(mu) || all(isnan(mu(:)))
                    mu = Status.Properties(['mu_', num2str(i)]).Value;
                    for j=1:obj.NofComponents-1
                        dmudz(:, j) = obj.Phases(i).ComputeDmuDz(Status.Properties('P_2').Value, Status.T, Status, obj.Thermo, dxdz, j);
                    end
                end
                % Compute d(Mob)/d(z)
                for j=1:obj.NofComponents-1
                    dMobdz(:, i, j) = (dkrdz(:, i, j) .* mu - dmudz(:, j) .* kr(:, i)) ./ mu.^2;
                end
            end
        end
        % function dMobdz = ComputeDMobDzHysteresis(obj, Status, dSdz, dxdz)
        %     N = length(dSdz);
        %     dMobdz = zeros(N, obj.NofPhases, obj.NofComponents-1);
        %     % Compute d(kr)/d(z)
        %     dkrdz = zeros(N, obj.NofPhases, obj.NofComponents-1);
        %     kr = obj.RelPermModel.kr;
        %     dkrdS = obj.RelPermModel.ComputeDerivativeHysteresis();
        %     for i=1:obj.NofPhases
        %         for j=1:obj.NofComponents-1
        %             dkrdz(:, i, j) = dkrdS(:, i) .* dSdz(:, j);
        %         end
        %     end
        %      % Compute d(mu)/d(z)
        %     for i=1:obj.NofPhases
        %         mu = obj.Phases(i).mu;
        %         dmudz = zeros(N, obj.NofComponents);
        %         if isempty(mu) || all(isnan(mu(:)))
        %             mu = Status.Properties(['mu_', num2str(i)]).Value;
        %             for j=1:obj.NofComponents-1
        %                 dmudz(:, j) = obj.Phases(i).ComputeDmuDz(Status.Properties(['P_', num2str(i)]).Value, Status.T, Status, obj.Thermo, dxdz, j);
        %             end
        %         end
        %         % Compute d(Mob)/d(z)
        %         for j=1:obj.NofComponents-1
        %             dMobdz(:, i, j) = (dkrdz(:, i, j) .* mu - dmudz(:, j) .* kr(:, i)) ./ mu.^2;
        %         end
        %     end
        % end
        % function dMobdz = ComputeDMobDz(obj, Status, dSdz, dxdz)
        %     N = length(dSdz);
        %     dMobdz = zeros(N, obj.NofPhases, obj.NofComponents-1);
        %     % Compute d(kr)/d(z)
        %     dkrdz = zeros(N, obj.NofPhases, obj.NofComponents-1);
        %     s = Status.Properties('S_1').Value;
        %     % kr = obj.RelPermModel.kr;
        %     kr = obj.RelPermModel.ComputeRelPerm(obj.Phases, Status);
        %     dkrdS = obj.RelPermModel.ComputeDerivative(obj.Phases, s);
        %     for i=1:obj.NofPhases
        %         for j=1:obj.NofComponents-1
        %             dkrdz(:, i, j) = dkrdS(:, i) .* dSdz(:, j);
        %         end
        %     end
        %      % Compute d(mu)/d(z)
        %     for i=1:obj.NofPhases
        %         mu = obj.Phases(i).mu;
        %         dmudz = zeros(N, obj.NofComponents);
        %         if isempty(mu) || all(isnan(mu(:)))
        %             mu = Status.Properties(['mu_', num2str(i)]).Value;
        %             for j=1:obj.NofComponents-1
        %                 dmudz(:, j) = obj.Phases(i).ComputeDmuDz(Status.Properties(['P_', num2str(i)]).Value, Status.T, Status, obj.Thermo, dxdz, j);
        %             end
        %         end
        %         % Compute d(Mob)/d(z)
        %         for j=1:obj.NofComponents-1
        %             dMobdz(:, i, j) = (dkrdz(:, i, j) .* mu - dmudz(:, j) .* kr(:, i)) ./ mu.^2;
        %         end
        %     end
        % end
        %% d(rho_T)/d(p)   &&   d(rho_T)/d(z)
        % Total density: rho_T = S_V * rho_V + S_L * rho_L
        function drhotdp = ComputeDrhotDp(obj, Status, drho, dS)
            N = length(Status.Properties('rho_1').Value);
            rho = zeros(N, obj.NofPhases);
            S = zeros(N, obj.NofPhases);
            for i = 1:obj.NofPhases
                rho(:,i) = Status.Properties(['rho_',num2str(i)]).Value;
                S(:,i) = Status.Properties(['S_',num2str(i)]).Value;
            end
            drhotdp = drho(:, 1) .* S(:, 1) + rho(:, 1) .* dS + drho(:, 2) .* S(:, 2) - rho(:, 2) .* dS;

            % When it's one phase derivative is the one of the existing phase
            drhotdp(S(:,1) == 1) = drho(S(:,1) == 1, 1);
            drhotdp(S(:,1) == 0) = drho(S(:,1) == 0, 2);
        end
        function drhotdz = ComputeDrhotDz(obj, Status, drho, dS)
            N = length(Status.Properties('rho_1').Value);
            rho = zeros(N, obj.NofPhases);
            S = zeros(N, obj.NofPhases);
            for i = 1:obj.NofPhases
                rho(:,i) = Status.Properties(['rho_',num2str(i)]).Value;
                S(:,i) = Status.Properties(['S_',num2str(i)]).Value;
            end
            drhotdz = zeros(length(S), obj.NofComponents-1);
            for j = 1:obj.NofComponents-1
                drhotdz(:,j) = drho(:, 1, j) .* S(:, 1) + rho(:, 1) .* dS(:, j) + drho(:, 2, j) .* S(:, 2) - rho(:, 2) .* dS(:, j);
            end
            % When it s one phase
            for j = 1:obj.NofComponents-1
                drhotdz(S(:,1) == 1, j) = drho(S(:,1) == 1, 1, j);
                drhotdz(S(:,1) == 0, j) = drho(S(:,1) == 0, 2, j);
            end
        end
        %% d(Pc)
        % function dPc = ComputeDPcDz(obj, Status, dSdz)
        %     dPcdS = obj.ComputeDPcDS(Status.Properties('S_1').Value);
        %     % dPcdS = obj.RelPermModel.dPc;
        % 
        %     dPc = zeros(length(dPcdS), obj.NofComponents-1);
        %     for j=1:obj.NofComponents-1
        %         % Use chain rule
        %         dPc(:,j) = dPcdS(:,1) .* dSdz(:,j);
        %     end
        % end
        function dPc = ComputeDPcDz(obj, Status, dSdz, label)
            dPcdS = obj.ComputeDPcDS(Status, label);
            % dPcdS = obj.RelPermModel.dPc;
            dPc = zeros(length(dPcdS), obj.NofComponents-1);
            for j=1:obj.NofComponents-1
                % Use chain rule
                dPc(:,j) = dPcdS(:,1) .* dSdz(:,j);
            end
        end
        function dPc = ComputeDPcDp(obj, Status, dSdp, label)
            dPcdS = obj.ComputeDPcDS(Status, label);
            % dPcdS = obj.RelPermModel.dPc;
            % Use chain rule
            dPc = dPcdS .* dSdp;
        end
    end
end