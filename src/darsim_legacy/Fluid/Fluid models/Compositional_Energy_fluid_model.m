% Compositional Energy fluid model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DARSim 2 Reservoir Simulator
% Author: Mengjie
%TU Delft
%Created: 02 March 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Compositional_Energy_fluid_model < Compositional_fluid_model
    properties

    end
    methods
        function obj = Compositional_Energy_fluid_model(n_phases, n_comp)
            obj@Compositional_fluid_model(n_phases, n_comp);
            obj.name = 'Compositional_Energy';
        end
        %% Initialize Injectors
        function InitializeInjectors(obj, Inj)  
            nc = length(obj.Components);
            % nph = length(obj.Phases);
            % Loop over all injectors
            for i=1:length(Inj)
                N = length(Inj{i}.Cells);
                c = Inj{i}.Composition;
                Inj{i}.z = zeros(N, nc);
                Inj{i}.z(:, c) = 1;
                %Inj{i}.x = zeros(N, nph * nc); 
                Inj{i}.x = [1 0 0 1];
                Inj{i}.rho(:,1) = obj.Thermo.rhoV_pTz_table(Inj{i}.p,Inj{i}.T, Inj{i}.z(:,1));
                Inj{i}.rho(:,2) = obj.Thermo.rhoAq_pT_table(Inj{i}.p,Inj{i}.T);
                Inj{i}.S = 1;

                Inj{i}.Mob = [1./obj.Thermo.muV_pT_table(Inj{i}.p,Inj{i}.T) zeros(size(Inj{i}.p))];
            end        
        end

        function ComputePhaseDensities(obj, Status, ~)
            for i=1:obj.NofPhases
                rho = Status.Properties(['rho_', num2str(i)]);
                rho.Value = obj.Phases(i).ComputeDensity(Status.Properties(['P_', num2str(i)]).Value, Status.Properties('T').Value, Status, obj.Thermo);
            end
        end

        function ComputePhaseEnthalpies(obj, Status)
            for i=1:obj.NofPhases
                PhaseEnthalpy = Status.Properties(['h_',num2str(i)]);
                PhaseEnthalpy.Value = obj.Phases(i).ComputeEnthalpy(Status.Properties(['P_', num2str(i)]).Value, Status.Properties('T').Value, Status, obj.Thermo);    
            end
        end
        
        function ComputeRockEnthalpy(obj, Medium)
            Hr = Medium.State.Properties('hRock');
            Hr.Value = Medium.Cpr .* (Medium.State.Properties('T').Value - obj.Thermo.T_ref); 
        end

        function Mob = ComputePhaseMobilitiesHysteresis(obj, Reservoir)
            Mob = zeros(length(Reservoir.State.Properties('z_1').Value), obj.NofPhases);
            kr = obj.RelPermModel.ComputeRelPermHysteresis(obj.Phases, Reservoir);
            for i=1:obj.NofPhases
                mu = Reservoir.State.Properties(['mu_', num2str(i)]);
                mu.Value = obj.Phases(i).ComputeViscosity(Reservoir.State.Properties(sprintf('P_%d', i)).Value, Reservoir.State.Properties('T').Value, Reservoir.State, obj.Thermo);
                Mob(:,i) = kr(:,i)./mu.Value;
            end
        end

        function Mob = ComputePhaseMobilities(obj, Medium)
            Mob = zeros(length(Medium.State.Properties('S_1').Value), obj.NofPhases);       
            kr = obj.RelPermModel.ComputeRelPerm(obj.Phases, Medium);
            for i=1:obj.NofPhases
                mu = obj.Phases(i).mu;
                if isempty(mu) || all(isnan(mu(:)))
                    mu = Medium.State.Properties(['mu_', num2str(i)]);
                    mu.Value = obj.Phases(i).ComputeViscosity(Medium.State.Properties(sprintf('P_%d', i)).Value, Medium.State.Properties('T').Value, Medium.State, obj.Thermo);
                    Mob(:,i) = kr(:,i)./mu.Value;
                else
                    Mob(:,i) = kr(:,i)./mu;
                end       
            end
        end

        function ComputeThermalConductivity(obj, Medium)
            % (1-phi)*lamda_r + phi*Sw*lamda_w + phi*Sg*lamda_g : vector
            CondEff = Medium.State.Properties('CondEff');
            CondEff.Value = (1 - Medium.Por) .* Medium.K_Cond_rock;
            for i=1:obj.NofPhases
                cond = obj.Phases(i).cond;
                lamda = Medium.State.Properties(['cond_',num2str(i)]);     
                S = Medium.State.Properties(['S_',num2str(i)]);
                lamda.Value = cond .* ones(size(S.Value));
                CondEff.Value = CondEff.Value + Medium.Por .* cond .* S.Value;
            end
        end

        %% K-value
        function dkdT = DKvalDT(obj, Status)
            dkdT = obj.FlashCalculator.KvaluesCalculator.DKvalDT(Status, obj.Thermo);
        end

        %% d(x)/d(p)   &&   d(x)/d(z)    &&   d(x)/d(T)
        function [dxdp, dxdz, dxdT] = DxDpDzDT(obj, Status, SinglePhase)
            % x = x1v, x1l, x2v, x2l, ..., ni
            % z = p, z1, z2, T
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
            dkdp = obj.DKvalDp(Status);
            dkdT = obj.DKvalDT(Status);
            
            % Allocating memories for dxdp, dxdz, dxdT
            dxdp = zeros(N, 2*obj.NofComponents+1);
            dxdz = zeros(N, 2*obj.NofComponents+1, obj.NofComponents-1);
            dxdT = zeros(N, 2*obj.NofComponents+1);

            % TwoPhase case
            num_ni = (k(:,2)-1)+ z.* (k(:,1)-k(:,2));
            den_ni = (k(:,1)-1) .*(k(:,2)-1);
            dnum_ni = (1 - z) .* dkdp(:,2) + z .* dkdp(:,1);
            dden_di = (k(:,2)-1) .* dkdp(:,1) + (k(:,1)-1) .* dkdp(:,2);
            dxdp(:, 5) = - (dnum_ni .* den_ni - dden_di .* num_ni) ./ den_ni.^2;
            dxdp(:, 2) =  - ( z .* ( dxdp(:,5) .* (k(:,1) - 1) + ni .* dkdp(:,1) ) ) ./ ( ( ni .* (k(:,1) - 1) + 1 ).^2 );
            dxdp(:, 1) = k(:,1) .* dxdp(:,2) + x(:,2) .* dkdp(:,1);
            dxdp(:, 3) = - dxdp(:,1);
            dxdp(:, 4) = - dxdp(:,2);

            dnum_dT = (1 - z) .* dkdT(:,2) + z .* dkdT(:,1);
            dden_dT = (k(:,2)-1) .* dkdT(:,1) + (k(:,1)-1) .* dkdT(:,2);
            dxdT(:, 5) = - (dnum_dT .* den_ni - dden_dT .* num_ni) ./ den_ni.^2;
            dxdT(:, 2) =  - ( z .* ( dxdT(:,5) .* (k(:,1) - 1) + ni .* dkdT(:,1) ) ) ./ ( ( ni .* (k(:,1) - 1) + 1 ).^2 );
            dxdT(:, 1) = k(:,1) .* dxdT(:,2) + x(:,2) .* dkdT(:,1);
            dxdT(:, 3) = - dxdT(:,1);
            dxdT(:, 4) = - dxdT(:,2);
    
            dxdz(:, 5) = (k(:,2) - k(:,1) )./ ((k(:,1)-1) .* (k(:,2)-1));
            dxdz(:, 2) = (ni .* (k(:,1) - 1) + 1 - z .* ( dxdz(:,5) .* (k(:,1) - 1) ) ) ./ ( ( ni .* (k(:,1) - 1) + 1 ).^2 );
            dxdz(:, 1) = k(:,1) .* dxdz(:,2);
            dxdz(:, 3) = - dxdz(:, 1);
            dxdz(:, 4) = - dxdz(:, 2);

            % SinglePhase case
            dxdp(SinglePhase == 1, :) = 0;
            dxdp(SinglePhase == 2, :) = 0;
            dxdT(SinglePhase == 1, :) = 0;
            dxdT(SinglePhase == 2, :) = 0;
            
            dxdz(SinglePhase == 1, 1) = 1; %1
            dxdz(SinglePhase == 1, 2) = 0;
            dxdz(SinglePhase == 1, 3) = -1; %-1
            dxdz(SinglePhase == 1, 4) = 0;
            dxdz(SinglePhase == 1, 5) = 0;
            
            dxdz(SinglePhase == 2, 1) = 0;
            dxdz(SinglePhase == 2, 2) = 1; %1
            dxdz(SinglePhase == 2, 3) = 0;
            dxdz(SinglePhase == 2, 4) = -1; %-1
            dxdz(SinglePhase == 2, 5) = 0;
        end


        %% d(rho)/d(p)   &&   d(rho)/d(z)   &&   d(rho)/d(T) 
        function drhodp = ComputeDrhoDp(obj, Status, dxdp)
            drhodp = zeros(length(dxdp), obj.NofPhases);
            for i=1:obj.NofPhases
                if i == 1
                    drhodp(:, i) = obj.Phases(i).ComputeDrhoDp(Status.Properties(['P_', num2str(i)]).Value, Status.Properties('T').Value, Status, obj.Thermo);
                else
                    drhodp(:, i) = obj.Phases(i).ComputeDrhoDp(Status.Properties(['P_', num2str(i)]).Value, Status.Properties('T').Value, Status, obj.Thermo, dxdp);
                end
            end     
        end
        function drhodz = ComputeDrhoDz(obj, Status, dxdz)
            drhodz = zeros(length(dxdz), obj.NofPhases);
            for i=1:obj.NofPhases
                if i == 1
                    drhodz(:, i) = obj.Phases(i).ComputeDrhoDz(Status.Properties(['P_', num2str(i)]).Value, Status.Properties('T').Value, Status, obj.Thermo);
                else
                    drhodz(:, i) = obj.Phases(i).ComputeDrhoDz(Status.Properties(['P_', num2str(i)]).Value, Status.Properties('T').Value, Status, obj.Thermo, dxdz);
                end
            end 
        end
        function drhodT = ComputeDrhoDT(obj, Status, dxdT)
            drhodT = zeros(length(dxdT), obj.NofPhases);
            for i=1:obj.NofPhases
                if i == 1
                    drhodT(:, i) = obj.Phases(i).ComputeDrhoDT(Status.Properties(['P_', num2str(i)]).Value, Status.Properties('T').Value, Status, obj.Thermo);
                else
                    drhodT(:, i) = obj.Phases(i).ComputeDrhoDT(Status.Properties(['P_', num2str(i)]).Value, Status.Properties('T').Value, Status, obj.Thermo, dxdT);
                end
            end 
        end


        %% d(H)/d(p)   &&   d(H)/d(z)   &&   d(H)/d(T) 
        function dhdp = ComputeDhDp(obj, Status, dxdp)
            dhdp = zeros(length(dxdp), obj.NofPhases);
            for i=1:obj.NofPhases
                if i == 1
                    dhdp(:, i) = obj.Phases(i).ComputeDhDp(Status.Properties(['P_', num2str(i)]).Value, Status.Properties('T').Value, Status, obj.Thermo);
                else
                    dhdp(:, i) = obj.Phases(i).ComputeDhDp(Status.Properties(['P_', num2str(i)]).Value, Status.Properties('T').Value, Status, obj.Thermo, dxdp);
                end
            end
        end

        function dhdz = ComputeDhDz(obj, Status, dxdz)
            dhdz = zeros(length(dxdz), obj.NofPhases);
            for i=1:obj.NofPhases
                if i == 1
                    dhdz(:, i) = obj.Phases(i).ComputeDhDz(Status.Properties(['P_', num2str(i)]).Value, Status.Properties('T').Value, Status, obj.Thermo);
                else
                    dhdz(:, i) = obj.Phases(i).ComputeDhDz(Status.Properties(['P_', num2str(i)]).Value, Status.Properties('T').Value, Status, obj.Thermo, dxdz);
                end
            end
        end

        function dhdT = ComputeDhDT(obj, Status, dxdT)
            dhdT = zeros(length(dxdT), obj.NofPhases);
            for i=1:obj.NofPhases
                if i == 1
                    dhdT(:, i) = obj.Phases(i).ComputeDhDT(Status.Properties(['P_', num2str(i)]).Value, Status.Properties('T').Value, Status, obj.Thermo);
                else
                    dhdT(:, i) = obj.Phases(i).ComputeDhDT(Status.Properties(['P_', num2str(i)]).Value, Status.Properties('T').Value, Status, obj.Thermo, dxdT);
                end
            end
        end


        %% d(Saturation)/d(p)   &&   d(Saturation)/d(z)   &&   d(Saturation)/d(T)
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
            drhov = drhodz(:,1);
            drhol = drhodz(:,2);
            for j = 1:obj.NofComponents-1
                Num1 = rhol .* ni;
                dNum1 = rhol .* dni(:, :, j) + drhol .* ni;
                Den1 = rhol .* ni + (1 - ni) .* rhov;
                dDen1 = rhol .* dni(:, :, j) + drhol .* ni + drhov - dni(:, :, j) .* rhov - ni .* drhov;
                dSdz(:, j) =(Den1 .* dNum1 - Num1 .* dDen1) ./ Den1.^2;
            end
            dSdz(Status.Properties('S_1').Value == 0, :) = 0;
            dSdz(Status.Properties('S_1').Value == 1, :) = 0;
            dSdz(Status.Properties('immobile_label').Value == 1, :) = 0;
        end

        function dSdT = ComputeDSDT(obj, Status, drhodT, dni)
            % Derivative of S with respect to overall composition
            rhov = Status.Properties('rho_1').Value;
            rhol = Status.Properties('rho_2').Value;
            ni = Status.Properties('ni_1').Value;
            Num = rhol .* ni;
            Den = rhol .* ni + (1 - ni) .* rhov;
            drhov = drhodT(:,1);
            drhol = drhodT(:,2);
            dNum = drhol .* ni + rhol .* dni;
            dDen = rhol .* dni + drhol .* ni + drhov - dni .* rhov - ni .* drhov;
            dSdT = (Den .* dNum - Num .* dDen) ./ Den.^2;
        end


        %% d(Mob)/dp   &&   d(Mob)/d(z)   &&   d(Mob)/d(T)
        function dMobdp = ComputeDMobDpHysteresis(obj, Status, dSdp, dxdp)
            kr = obj.RelPermModel.kr;
            dkrdS = obj.RelPermModel.ComputeDerivativeHysteresis(); %%%%%%%
            dkr = dkrdS .* dSdp;
            dMobdp = zeros(length(kr), obj.NofPhases);
            for i=1:obj.NofPhases
                mu = Status.Properties(['mu_', num2str(i)]).Value;
                if i == 1            
                    dmudp = obj.Phases(i).ComputeDmuDp(Status.Properties(['P_', num2str(i)]).Value, Status.Properties('T').Value, obj.Thermo);                    
                else
                    dmudp = obj.Phases(i).ComputeDmuDp(Status.Properties(['P_', num2str(i)]).Value, Status.Properties('T').Value, Status, obj.Thermo, dxdp);                 
                end
                dMobdp(:,i) = (dkr(:,i) .* mu - dmudp .* kr(:,i)) ./ mu.^2;
            end
        end
        function dMobdp = ComputeDMobDp(obj, Status, dSdp, dxdp)
            kr = obj.RelPermModel.kr;
            s = Status.Properties('S_1').Value;
            dkrdS = obj.RelPermModel.ComputeDerivative(obj.Phases, s); %%%%%%%
            dkr = dkrdS .* dSdp;
            dMobdp = zeros(length(kr), obj.NofPhases);
            for i=1:obj.NofPhases
                mu = Status.Properties(['mu_', num2str(i)]).Value;
                if i == 1            
                    dmudp = obj.Phases(i).ComputeDmuDp(Status.Properties(['P_', num2str(i)]).Value, Status.Properties('T').Value, obj.Thermo);                    
                else
                    dmudp = obj.Phases(i).ComputeDmuDp(Status.Properties(['P_', num2str(i)]).Value, Status.Properties('T').Value, Status, obj.Thermo, dxdp);                 
                end
                dMobdp(:,i) = (dkr(:,i) .* mu - dmudp .* kr(:,i)) ./ mu.^2;
            end
        end
        
        %% d(Mob)/d(z)
        function dMobdz = ComputeDMobDzHysteresis(obj, Status, dSdz, dxdz)
            kr = obj.RelPermModel.kr;
            dkrdS = obj.RelPermModel.ComputeDerivativeHysteresis(); %%%%%%%
            dkr = dkrdS .* dSdz;
            dMobdz = zeros(length(kr), obj.NofPhases);
            for i=1:obj.NofPhases
                mu = Status.Properties(['mu_', num2str(i)]).Value;
                if i == 1            
                    dmudz = obj.Phases(i).ComputeDmuDz(Status.Properties(['P_', num2str(i)]).Value, Status.Properties('T').Value, obj.Thermo);                    
                else
                    dmudz = obj.Phases(i).ComputeDmuDz(Status.Properties(['P_', num2str(i)]).Value, Status.Properties('T').Value, Status, obj.Thermo, dxdz);                 
                end
                dMobdz(:,i) = (dkr(:,i) .* mu - dmudz .* kr(:,i)) ./ mu.^2;
            end
        end
        function dMobdz = ComputeDMobDz(obj, Status, dSdz, dxdz)
            kr = obj.RelPermModel.kr;
            s = Status.Properties('S_1').Value;
            dkrdS = obj.RelPermModel.ComputeDerivative(obj.Phases, s); %%%%%%%
            dkr = dkrdS .* dSdz;
            dMobdz = zeros(length(kr), obj.NofPhases);
            for i=1:obj.NofPhases
                mu = Status.Properties(['mu_', num2str(i)]).Value;
                if i == 1            
                    dmudz = obj.Phases(i).ComputeDmuDz(Status.Properties(['P_', num2str(i)]).Value, Status.Properties('T').Value, obj.Thermo);                    
                else
                    dmudz = obj.Phases(i).ComputeDmuDz(Status.Properties(['P_', num2str(i)]).Value, Status.Properties('T').Value, Status, obj.Thermo, dxdz);                 
                end
                dMobdz(:,i) = (dkr(:,i) .* mu - dmudz .* kr(:,i)) ./ mu.^2;
            end
        end

        %% d(Mob)/d(T)
        function dMobdT = ComputeDMobDTHysteresis(obj, Status, dSdT, dxdT)
            kr = obj.RelPermModel.kr;
            dkrdS = obj.RelPermModel.ComputeDerivativeHysteresis(); %%%%%%%
            dkr = dkrdS .* dSdT;
            dMobdT = zeros(length(kr), obj.NofPhases);
            for i=1:obj.NofPhases
                mu = Status.Properties(['mu_', num2str(i)]).Value;
                if i == 1            
                    dmudT = obj.Phases(i).ComputeDmuDT(Status.Properties(['P_', num2str(i)]).Value, Status.Properties('T').Value,  obj.Thermo);                    
                else
                    dmudT = obj.Phases(i).ComputeDmuDT(Status.Properties(['P_', num2str(i)]).Value, Status.Properties('T').Value, Status, obj.Thermo, dxdT);                 
                end
                dMobdT(:,i) = (dkr(:,i) .* mu - dmudT .* kr(:,i)) ./ mu.^2;
            end
        end
        function dMobdT = ComputeDMobDT(obj, Status, dSdT, dxdT)
            kr = obj.RelPermModel.kr;
            s = Status.Properties('S_1').Value;
            dkrdS = obj.RelPermModel.ComputeDerivative(obj.Phases, s); %%%%%%%
            dkr = dkrdS .* dSdT;
            dMobdT = zeros(length(kr), obj.NofPhases);
            for i=1:obj.NofPhases
                mu = Status.Properties(['mu_', num2str(i)]).Value;
                if i == 1            
                    dmudT = obj.Phases(i).ComputeDmuDT(Status.Properties(['P_', num2str(i)]).Value, Status.Properties('T').Value, obj.Thermo);                    
                else
                    dmudT = obj.Phases(i).ComputeDmuDT(Status.Properties(['P_', num2str(i)]).Value, Status.Properties('T').Value, Status, obj.Thermo, dxdT);                 
                end
                dMobdT(:,i) = (dkr(:,i) .* mu - dmudT .* kr(:,i)) ./ mu.^2;
            end
        end


        %% d(rho_T)/d(p)   &&   d(rho_T)/d(z) (in Compositional_fluid_model.m)   &&   d(rho_T)/d(T)
        function drhotdT = ComputeDrhotDT(obj, Status, drho, dS)
            N = length(Status.Properties('rho_1').Value);
            rho = zeros(N, obj.NofPhases);
            S = zeros(N, obj.NofPhases);
            for i = 1:obj.NofPhases
                rho(:,i) = Status.Properties(['rho_',num2str(i)]).Value;
                S(:,i) = Status.Properties(['S_',num2str(i)]).Value;
            end
            drhotdT = drho(:, 1) .* S(:, 1) + rho(:, 1) .* dS + drho(:, 2) .* S(:, 2) - rho(:, 2) .* dS;

            % When it's one phase derivative is the one of the existing phase
            drhotdT(S(:,1) == 1) = drho(S(:,1) == 1, 1);
            drhotdT(S(:,1) == 0) = drho(S(:,1) == 0, 2);
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
    end
end