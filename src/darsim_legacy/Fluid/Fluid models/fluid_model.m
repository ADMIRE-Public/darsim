 % Fluid model base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini, Yuhang Wang
%TU Delft
%Created: 14 July 2016
%Last modified: 28 Jan 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef fluid_model < handle
    properties
        name
        NofPhases
        NofComponents
        Phases
        Components
        RelPermModel
        CapillaryModel
        WettingPhaseIndex
        HysteresisActive

        RelPermModelFractures
        
        Mob_fluid_model

        Thermo
    end
    methods
        function obj = fluid_model(n_phases, n_comp)
            obj.NofPhases = n_phases;
            obj.NofComponents = n_comp;
            obj.Phases = phase.empty;
            obj.Components = component.empty;
            obj.WettingPhaseIndex = 1;
        end
        function AddPhase(obj, Phase, index)
            obj.Phases(index) = Phase;
        end
        function AddComponent(obj, Comp, index)
            obj.Components(index) = Comp;
        end
        function Mob = ComputePhaseMobilitiesWells(obj, S)
            Mob = zeros(length(S), obj.NofPhases);
            kr = [1 0];
            for i=1:obj.NofPhases
                Mob(:,i) = kr(:,i)./obj.Phases(i).mu;
            end
        end
        function Mob = ComputePhaseMobilities(obj, Medium)
            Mob = zeros(length(Medium.State.Properties('S_1').Value), obj.NofPhases);   
            if isa(Medium, 'reservoir')
                if obj.HysteresisActive
                    kr = obj.RelPermModel.ComputeRelPermHysteresis(obj.Phases, Medium);
                else
                    kr = obj.RelPermModel.ComputeRelPerm(obj.Phases, Medium);
                end
            elseif isa(Medium, 'fracture')
                kr = obj.RelPermModelFractures.ComputeRelPerm(obj.Phases, Medium, obj.RelPermModel);
            else
                error('Medium is not recognized.')
            end

            for i=1:obj.NofPhases
                mu = obj.Phases(i).mu;
                if isempty(mu) || all(isnan(mu(:)))
                    mu = Medium.State.Properties(['mu_', num2str(i)]);
                    mu.Value = obj.Phases(i).ComputeViscosity(Medium.State.Properties('P_2').Value, Medium.State.T, Medium.State, obj.Thermo);
                    Mob(:,i) = kr(:,i)./mu.Value;
                else
                    Mob(:,i) = kr(:,i)./mu;
                end       
            end
        end
        % function Mob = ComputePhaseMobilitiesHysteresis(obj, Reservoir)
        %     Mob = zeros(length(Reservoir.State.Properties('z_1').Value), obj.NofPhases);
        %     kr = obj.RelPermModel.ComputeRelPermHysteresis(obj.Phases, Reservoir);
        %     for i=1:obj.NofPhases
        %         mu = obj.Phases(i).mu;
        %         if isempty(mu) || all(isnan(mu(:)))
        %             mu = Reservoir.State.Properties(['mu_', num2str(i)]);
        %             mu.Value = obj.Phases(i).ComputeViscosity(Reservoir.State.Properties(sprintf('P_%d', i)).Value, Reservoir.State.T, Reservoir.State, obj.Thermo);
        %             Mob(:,i) = kr(:,i)./mu.Value;
        %         else
        %             Mob(:,i) = kr(:,i)./mu;
        %         end       
        %     end
        % end
        % function Mob = ComputePhaseMobilitiesFractures(obj, Medium)
        %     Mob = zeros(length(Medium.State.Properties('S_1').Value), obj.NofPhases);       
        %     kr = obj.RelPermModelFractures.ComputeRelPerm(obj.Phases, Medium.State);
        %     for i=1:obj.NofPhases
        %         mu = obj.Phases(i).mu;
        %         if isempty(mu) || all(isnan(mu(:)))
        %             mu = Medium.State.Properties(['mu_', num2str(i)]);
        %             mu.Value = obj.Phases(i).ComputeViscosity(Medium.State.Properties(sprintf('P_%d', i)).Value, Medium.State.T, Medium.State, obj.Thermo);
        %             Mob(:,i) = kr(:,i)./mu.Value;
        %         else
        %             Mob(:,i) = kr(:,i)./mu;
        %         end       
        %     end
        % end
        function ComputeTotalDensity(obj, State)
            % Compute the total density
            s = State.Properties('S_1');
            rho = State.Properties('rho_1');
            rhoT = rho.Value .* s.Value;
            for i=2:obj.NofPhases
                s = State.Properties(['S_', num2str(i)]);
                rho = State.Properties(['rho_', num2str(i)]);
                rhoT = rhoT + rho.Value .* s.Value;
            end
            RhoT = State.Properties('rhoT');
            RhoT.Value = rhoT;
        end
        function dMobdS = ComputeDMobDS(obj, s)
            dMobdS = zeros(length(s), obj.NofPhases);
            dkr = obj.RelPermModel.ComputeDerivative(obj.Phases, s); %%%%%%
            for i=1:obj.NofPhases
                dMobdS(:,i) = dkr(:,i)/obj.Phases(i).mu;
            end
        end
        function dMobdS = ComputeDMobDSHysteresis(obj, Status)
            s = Status.Properties('S_1').Value;
            kr = obj.RelPermModel.kr;
            dMobdS = zeros(length(s), obj.NofPhases);      
            dkr = obj.RelPermModel.ComputeDerivativeHysteresis(); %%%%%%%
            for i=1:obj.NofPhases
                mu = obj.Phases(i).mu;
                if isempty(mu) || all(isnan(mu(:)))
                    mu = Status.Properties(['mu_', num2str(i)]).Value;
                    dMobdS(:,i) = dkr(:,i) ./ mu;
                else
                    dMobdS(:,i) = dkr(:,i) ./ mu;
                end
                
            end
        end
        function dMobdSdS = ComputeDMobDSDS(obj, s)
            dMobdSdS = zeros(length(s), obj.NofPhases);
            ddkr = obj.RelPermModel.ComputeSecondDerivative(obj.Phases, s);
            for i=1:obj.NofPhases
                dMobdSdS(:,i) = ddkr(:,i)/obj.Phases(i).mu;
            end
        end
        function dMobdz = ComputeDMobDzHysteresis(obj, Status, dSdz)
            dMobdS = obj.ComputeDMobDSHysteresis(Status);
            dMobdz = zeros(length(dMobdS), obj.NofComponents-1);
            for j=1:obj.NofComponents-1
                % Use chain rule
                dMobdz(:,1,j) = dMobdS(:,1) .* dSdz(:,j);
                dMobdz(:,2,j) = dMobdS(:,2) .* dSdz(:,j);
            end
        end
        function dMobdz = ComputeDMobDz(obj, Status, dSdz)
            dMobdS = obj.ComputeDMobDS(Status.Properties('S_1').Value);
            dMobdz = zeros(length(dMobdS), obj.NofComponents-1);
            for j=1:obj.NofComponents-1
                % Use chain rule
                dMobdz(:,1,j) = dMobdS(:,1) .* dSdz(:,j);
                dMobdz(:,2,j) = dMobdS(:,2) .* dSdz(:,j);
            end
        end
        function drhodp = ComputeDrhoDp(obj, Status, SinglePhase, dxdp)
            drhodp = zeros(length(Status.Properties('P_1').Value), obj.NofPhases);
            for i=1:obj.NofPhases
                if i == 1
                    drhodp(:, i) = obj.Phases(i).ComputeDrhoDp(Status.Properties('P_1').Value, obj.Components);
                    % drhodp(:, i) = obj.Phases(i).ComputeDrhoDp(Status.Properties(['P_', num2str(i)]).Value, obj.Components, obj.Rs(:,i), obj.dRs(:,i));
                else
                    drhodp(:, i) = obj.Phases(i).ComputeDrhoDp(Status.Properties('P_2').Value, obj.Components, dxdp, Status.Properties('x_1ph2').Value);
                end
            end
        end
        %% original ComputePc function
        function ComputePc_old(obj, Status, DiscretizationModel, label)
                % Define local handles
                S = Status.Properties('S_1').Value;
                P1 = Status.Properties('P_1');
                P2 = Status.Properties('P_2');
                Pc = Status.Properties('Pc');             
            if (label==1) % Matrix                
                if obj.HysteresisActive
                    obj.CapillaryModel.ComputePcHysteresis(obj.RelPermModel, DiscretizationModel, Status);
                else
                    % Compute Pc as a function of S_1
                    switch(obj.WettingPhaseIndex)
                        case(1)
                            s = (S - obj.Phases(1).sr)./(1 - obj.Phases(1).sr);
                            s = max(s, 0.05);
                            Pc.Value = obj.CapillaryModel.ComputePc(s);
                            Pc.Value(S < obj.Phases(1).sr) = 0;
                        case(2)
                            S = 1 - S;
                            s = (S - obj.Phases(2).sr)./(1 - obj.Phases(2).sr - obj.Phases(1).sr); % s = (S - obj.Phases(2).sr)./(1 - obj.Phases(2).sr)
                            s = max(s, 0.05); % 0.05
                            s = min(s, 0.99); % 0.99
                            Pc.Value = -obj.CapillaryModel.ComputePc(s, label);
                            Pc.Value(S < obj.Phases(2).sr) = 0;
                    end
                    % Update P_1
                    P1.Value = P2.Value - Pc.Value;
                end  
            else % Fractures
                % no pillarity
                Pc.Value = zeros(length(S),1);
                % Update P_1
                P1.Value = P2.Value;
                
                % obj.CapillaryModel.ComputePcnoHysteresis(obj.RelPermModel, Status);
            end
        end  
        %
        function dPcdS = ComputeDPcDS_old(obj, Status, label) 
            if (label==1) % Matrix
                if obj.HysteresisActive
                    dPcdS = obj.CapillaryModel.dPcdSHysteresis();
                else
                    S = Status.Properties('S_1').Value;
                    switch(obj.WettingPhaseIndex)                        
                        case(1)
                            s = (S - obj.Phases(1).sr)./(1 - obj.Phases(1).sr);
                            s = max(s, 0.05);
                            dPcdS = obj.CapillaryModel.dPcdS(s);
                            dPcdS (S < obj.Phases(1).sr) = 0.0;
                        case(2)
                            S = 1 - S;
                            s = (S - obj.Phases(2).sr)./(1 - obj.Phases(2).sr - obj.Phases(1).sr); % s = (S - obj.Phases(2).sr)./(1 - obj.Phases(2).sr);
                            s = max(s, 0.05);
                            s = min(s, 0.99);
                            % There is a double negative sign. P1 = P2 - Pc. Pc = -Pc so dPc = -dPc and dPcdS1 = -dPcdS2
                            % dPcdS = obj.CapillaryModel.dPcdS(s);
                            dPcdS = (1 - obj.Phases(2).sr)^(-1)*obj.CapillaryModel.dPcdS(s, label);
                            dPcdS (S < obj.Phases(2).sr) = 0.0;
                    end
                end
             
            else % Fractures
                % no pillarity
               dPcdS = zeros(length(Status.Properties('S_1').Value),1);   
                             
               % dPcdS = obj.CapillaryModel.dPcdSnoHysteresis(obj.RelPermModel, Status);
            end
        end
        %% Call this for no hysteresis
        function ComputePc(obj, Medium, DiscretizationModel)
                % Define local handles
                S = Medium.State.Properties('S_1').Value;
                P1 = Medium.State.Properties('P_1');
                P2 = Medium.State.Properties('P_2');
                Pc = Medium.State.Properties('Pc');   
                if isa(Medium, 'reservoir')
                    if obj.HysteresisActive
                        obj.CapillaryModel.ComputePcHysteresis(obj.RelPermModel, Medium, DiscretizationModel);
                    else
                        switch(obj.name)
                            case{'BlackOil','Compositional','Compositional_Energy'}
                                obj.CapillaryModel.ComputePcnoHysteresis(obj.RelPermModel, Medium, DiscretizationModel);
                            otherwise
                                % Compute Pc as a function of S_1
                                switch(obj.WettingPhaseIndex)
                                    case(1)
                                        s = (S - obj.Phases(1).sr)./(1 - obj.Phases(1).sr);
                                        s = max(s, 0.05);
                                        Pc.Value = obj.CapillaryModel.ComputePc(s);
                                        Pc.Value(S < obj.Phases(1).sr) = 0;
                                    case(2)
                                        S = 1 - S;
                                        s = (S - obj.Phases(2).sr)./(1 - obj.Phases(2).sr - obj.Phases(1).sr); % s = (S - obj.Phases(2).sr)./(1 - obj.Phases(2).sr)
                                        s = max(s, 0.05); % 0.05
                                        s = min(s, 0.99); % 0.99
                                        Pc.Value = -obj.CapillaryModel.ComputePc(s, label);
                                        Pc.Value(S < obj.Phases(2).sr) = 0;
                                end
                                % Update P_1
                                P1.Value = P2.Value - Pc.Value;
                        end
                    end
                elseif isa(Medium, 'fracture')
                    % no pillarity
                    Pc.Value = zeros(length(S),1);
                    % Update P_1
                    P1.Value = P2.Value;
                else
                    error('Medium is not recognized.')
                end
        end 

        function dPcdS = ComputeDPcDS(obj, Status, label) 
            if (label==1) % Matrix
                if obj.HysteresisActive
                    dPcdS = obj.CapillaryModel.dPcdSHysteresis();
                else
                    dPcdS = obj.CapillaryModel.dPcdSnoHysteresis(obj.RelPermModel, Status);
                end
             
            else % Fractures
                % no pillarity
               dPcdS = zeros(length(Status.Properties('S_1').Value),1);   
                             
               % dPcdS = obj.CapillaryModel.dPcdSnoHysteresis(obj.RelPermModel, Status);
            end
        end
        %%
        function dPcdS = ComputeDPcDSDS(obj, S)
            switch(obj.WettingPhaseIndex)
                case(1)
                    s = (S - obj.Phases(1).sr)./(1 - obj.Phases(1).sr);
                    s = max(s, 0.05);
                    dPcdS = obj.CapillaryModel.dPcdSdS(s);
                    dPcdS (S < obj.Phases(1).sr) = 0.0;
                case(2)
                    S = 1 - S;
                    s = (S - obj.Phases(2).sr)./(1 - obj.Phases(2).sr);
                    s = max(s, 0.05);
                    % There is a double negative sign. P1 = P2 - Pc. Pc = -Pc so dPc = -dPc and dPcdS1 = -dPcdS2
                    dPcdS = obj.CapillaryModel.dPcdSdS(s);  
                    dPcdS (S < obj.Phases(2).sr) = 0.0;
            end
            
        end
        function drhotdp = ComputeDrhotDp(obj, Status, drho, dS)
            N = length(Status.Properties('rho_1').Value);
            rho = zeros(N, obj.NofPhases);
            S = zeros(N, obj.NofPhases);
            for i=1:obj.NofPhases
                rho(:,i) = Status.Properties(['rho_',num2str(i)]).Value;
                S(:,i) = Status.Properties(['S_',num2str(i)]).Value;
            end
            drhotdp = drho(:,1) .* S(:,1) + rho(:,1) .* dS + drho(:,2) .* S(:,2) - rho(:,2) .* dS;
            % When it s one phase derivative is the one of the existing
            % phase
            drhotdp(S(:,1) == 1) = drho(S(:,1) == 1, 1);
            drhotdp(S(:,1) == 0) = drho(S(:,1) == 0, 2);
        end
    end
    methods (Abstract)
        obj = Flash(obj);
        obj = InitializeInjectors(obj);
    end
end
