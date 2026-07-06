% Immiscible Fluid model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini, Barnaby Fryer, Yuhang Wang
%TU Delft
%Created: 14 July 2016
%Last modified: 1 Feb 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef BO_fluid_model < Compositional_fluid_model
    properties
        Rs
        dRs
    end
    methods
        function obj = BO_fluid_model(n_phases, n_comp)
            obj@Compositional_fluid_model(n_phases, n_comp);
            obj.name = 'BlackOil';
        end
        function SinglePhase = Flash(obj, Status)
           SinglePhase = obj.FlashCalculator.Flash(Status, obj.Components, obj.Phases);
        end
        function InitializeInjectors(obj, Inj)
            % Loop over all injectors
            for i=1:length(Inj)
                Inj{i}.z = [1 0];
                Inj{i}.x = [1 0 0 1];        
                for ph=1:obj.NofPhases
                    Inj{i}.rho(:, ph) = obj.Phases(ph).ComputeDensity(Inj{i}.p, obj.Components, zeros(length(Inj{i}.Cells), 1));
                end
                Inj{i}.S = 1;
                % Inj{i}.Mob = obj.ComputePhaseMobilities(Inj{i}.S);
                                Inj{i}.Mob = [1/obj.Phases(1).mu 0];
            end            
        end
        function ComputePhaseDensities(obj, Status, SinglePhase)
            [obj.Rs, obj.dRs] = obj.FlashCalculator.KvaluesCalculator.ComputeRs(Status, obj.Phases);
            for i=1:obj.NofPhases
                % Correct Rs and dRs for undersaturated cells
                [obj.Rs(:,i), obj.dRs(:,i)] = obj.Phases(i).RsOfUnderSaturatedPhase(Status.Properties('z_1').Value, obj.Components, obj.Rs(:,i), obj.dRs(:,i), SinglePhase);
                rho = Status.Properties(['rho_', num2str(i)]);
                if i == 1
                    if all(Status.Properties('P_1').Value == 0)
                        rho.Value = obj.Phases(i).ComputeDensity(Status.Properties('P_2').Value, obj.Components);
                    else
                        rho.Value = obj.Phases(i).ComputeDensity(Status.Properties('P_1').Value, obj.Components);
                    end
                else
                    rho.Value = obj.Phases(i).ComputeDensity(Status.Properties('P_2').Value, obj.Components, obj.Rs(:,i));

                    % rho.Value = obj.Phases(i).ComputeDensity(Status.Properties(['P_', num2str(i)]).Value, obj.Components, obj.Rs(:,i));
                end
                % rho.Value = obj.Phases(i).ComputeDensity(Status.Properties(['P_', num2str(i)]).Value, obj.Components, obj.Rs(:,i));
            end
            status_Rs = Status.Properties('Rs'); % update Rs value
            status_Rs.Value = obj.Rs(:,2);
        end
        function k = ComputeKvalues(obj, Status)            
            k = obj.FlashCalculator.KvaluesCalculator.Compute(Status, obj.Components, obj.Phases);
        end
        function dkdp = DKvalDp(obj, Status)
            dkdp = obj.FlashCalculator.KvaluesCalculator.DKvalDp(Status, obj.Components, obj.Phases);
        end
        function drhodp = ComputeDrhoDp(obj, Status, SinglePhase, dxdp)
            drhodp = zeros(length(SinglePhase), obj.NofPhases);
            [obj.Rs, obj.dRs] = obj.FlashCalculator.KvaluesCalculator.ComputeRs(Status, obj.Phases);
            for i=1:obj.NofPhases
                [obj.Rs(:,i), obj.dRs(:,i)] = obj.Phases(i).RsOfUnderSaturatedPhase(Status.Properties('z_1').Value, obj.Components, obj.Rs(:,i), obj.dRs(:,i), SinglePhase);
                if i == 1
                drhodp(:, i) = obj.Phases(i).ComputeDrhoDp(Status.Properties('P_2').Value, obj.Components);
                    
                else
                drhodp(:, i) = obj.Phases(i).ComputeDrhoDp(Status.Properties('P_2').Value, obj.Components, obj.Rs(:,i), obj.dRs(:,i));
                end
                % drhodp(:, i) = obj.Phases(i).ComputeDrhoDp(Status.Properties(['P_', num2str(i)]).Value, obj.Components, obj.Rs(:,i), obj.dRs(:,i));                  
            end
        end
        function drhodz = ComputeDrhoDz(obj, Status, SinglePhase, dxdz)
            N = length(Status.Properties('P_2').Value);
            p = Status.Properties('P_2').Value;
            z = Status.Properties('z_1').Value;
            drhodz = zeros(N, obj.NofPhases); % currently rho does not depend on z;
            % drhodz(:, 2) = obj.Phases(2).ComputeDrhoDz(p, z, obj.Components, SinglePhase);
        end
    end
end