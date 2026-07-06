% Table relative permeability model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Yuhang Wang
%TU Delft
%Created: 6 Jan 2021
%Last modified: 30 June 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef relperm_model_table < relperm_model
    properties
        krg_Sg_table
        nonwetting_imbibation_table
        krw_Sg_table
        dkrg_Sg_table
        dkrw_Sg_table

        Pc_drainage_Sg_table
        dPc_drainage_Sg_table
        Pc_imbibation_Sg_table
        dPc_imbibation_Sg_table

        scanning_curve_kr_table
        scanning_curve_dkr_table
        scanning_curve_Pc_table
        scanning_curve_dPc_table

        % hysteresis parameters
        % safeguard_S_label
        % transport_label
        kr
        dkr
        % Pc
        % dPc
        scanning_curve_Sg
        scanning_curve_kr
        scanning_curve_dkr
        scanning_curve_Pc
        scanning_curve_dPc
        scanning_curve_Sgi

        reverse_point
        process_label % 0: primary drainage; 1: scanning curve
        immobile_label % 0; mobile; 1: immobile
        Sg_r % residual saturation of gas phase

        C % Land trapping coefficient
        S_gcr % critical gas phase saturation
        S_gimax % maximum gas phase saturation
        S_gtmax % maximum residual saturation of gas phase

        Sgr_Facies % Nonwetting phase residual saturation of each facies
        Swr_Facies % Wetting phase residual saturation of each facies

    end
    properties (Constant)
        NScanningCurves = 1000
        NScanningPoints = 35

        fr = 0.5 % Residual Gas Fraction: S_gt = S_gi * fr
        PcScanningEpsilon = 0.1

        MinReversalSaturation = 2.0e-3
        ReversalTolerance = 1.0e-4
    end
    methods

        %% Initialize relperm_model_table
        function obj = relperm_model_table(SimulationInput)

            [N, FaciesIndex, N_Facies, rawTables] = obj.GetModelInfo(SimulationInput);

            obj.AllocateArrays(N, obj.NScanningCurves, obj.NScanningPoints);

            for faciesID = 1:N_Facies

                [wetting_kr, nonwetting_drainage_kr, nonwetting_imbibation_kr, Pc_drainage, Pc_imbibation] = obj.InitialiseFaciesProperties(rawTables, FaciesIndex, faciesID);

                obj.BuildRelPermInterpolations(faciesID, wetting_kr, nonwetting_drainage_kr, nonwetting_imbibation_kr);

                obj.BuildCapillaryPressureInterpolations(faciesID, Pc_drainage, Pc_imbibation);

                obj.BuildScanningCurveInterpolations(faciesID, obj.NScanningCurves, obj.NScanningPoints);

            end
        end

        %% Read lookup tables
        function [N, FaciesIndex, N_Facies, relpermTables] = GetModelInfo(obj, SimulationInput)

            N = SimulationInput.ReservoirProperties.Grid.N_ActiveCells;

            FaciesIndex = SimulationInput.ReservoirProperties.FaciesIndex;

            N_Facies = max(FaciesIndex);

            if N_Facies ~= size(SimulationInput.FluidProperties.RelPerm.Files, 2)
                error("Number of the facies doesn't match with the number of the lookup tables.")
            end

            relpermTables = cell(N_Facies, 1);
            for faciesID = 1:N_Facies
                file = fullfile(SimulationInput.Directory, SimulationInput.FluidProperties.RelPerm.Files{faciesID});
                relpermTables{faciesID} = read_table(file);
            end

        end

        %% Allocate Arrays
        function AllocateArrays(obj, N, N_scanning_curves, N_scanning_points)
            obj.kr = zeros(N, 2);
            obj.dkr = zeros(N, 2);

            obj.scanning_curve_Sg  = zeros(N_scanning_curves, N_scanning_points);
            obj.scanning_curve_kr  = zeros(N_scanning_curves, N_scanning_points);
            obj.scanning_curve_dkr = zeros(N_scanning_curves, N_scanning_points);
            obj.scanning_curve_Pc  = zeros(N_scanning_curves, N_scanning_points);
            obj.scanning_curve_dPc = zeros(N_scanning_curves, N_scanning_points);

            obj.reverse_point = zeros(N, 1);
            obj.process_label = zeros(N, 1);
            obj.immobile_label = zeros(N, 1);
            obj.Sg_r = zeros(N, 1);
        end

        %% Initialise Facies Properties
        function [wetting_kr, nonwetting_drainage_kr, nonwetting_imbibation_kr, Pc_drainage, Pc_imbibation] = InitialiseFaciesProperties(obj, rawTables, FaciesIndex, faciesID)

            % relations are functions of resident wetting phase;
            wetting_kr = rawTables{faciesID}.wetting_kr;
            nonwetting_drainage_kr = rawTables{faciesID}.nonwetting_drainage_kr;
            nonwetting_imbibation_kr = rawTables{faciesID}.nonwetting_imbibation_kr;
            Pc_drainage = rawTables{faciesID}.Pc_drainage;
            Pc_imbibation = rawTables{faciesID}.Pc_imbibation;

            obj.Sgr_Facies{faciesID} = 1 - nonwetting_imbibation_kr(end, 1); % Gas phase residual saturation based on table #i
            obj.Swr_Facies{faciesID} = wetting_kr(1, 1); % Liquid phase residual saturation based on table #i
            Index = faciesID == FaciesIndex;
            obj.reverse_point(Index) = 1 - obj.Swr_Facies{faciesID};
            % obj.Sg_r(Index) = obj.Sgr_Facies{faciesID};
            obj.S_gimax{faciesID} = 1 - obj.Swr_Facies{faciesID};
            obj.S_gtmax{faciesID} = obj.S_gimax{faciesID} * obj.fr;

            % relative permeability
            % convert relations to functions of S_g;
            wetting_kr(:,1) = 1.0 - wetting_kr(:,1);
            nonwetting_drainage_kr(:,1) = 1.0 - nonwetting_drainage_kr(:,1);
            nonwetting_imbibation_kr(:,1) = 1.0 - nonwetting_imbibation_kr(:,1);

            wetting_kr = flip(wetting_kr); % stored in ascending order
            nonwetting_imbibation_kr = flip(nonwetting_imbibation_kr); % stored in ascending order
            nonwetting_drainage_kr = flip(nonwetting_drainage_kr); % stored in ascending order
        end

        %% Build RelPerm Interpolations
        function BuildRelPermInterpolations(obj, faciesID, wetting_kr, nonwetting_drainage_kr, nonwetting_imbibation_kr)

            obj.nonwetting_imbibation_table{faciesID} = griddedInterpolant(nonwetting_imbibation_kr(:,1),nonwetting_imbibation_kr(:,2));

            % Making interpolations
            obj.krg_Sg_table{faciesID} = griddedInterpolant(nonwetting_drainage_kr(:,1),nonwetting_drainage_kr(:,2));

            obj.krw_Sg_table{faciesID}  = griddedInterpolant(wetting_kr(:,1),wetting_kr(:,2));

            dkrg_Sg_value = (nonwetting_drainage_kr(2:end,2)-nonwetting_drainage_kr(1:end-1,2))./(nonwetting_drainage_kr(2:end,1)-nonwetting_drainage_kr(1:end-1,1));
            dkrg_Sg_value = [dkrg_Sg_value(1);dkrg_Sg_value];
            obj.dkrg_Sg_table{faciesID} = griddedInterpolant(nonwetting_drainage_kr(:,1),dkrg_Sg_value);

            dkro_Sg_value = (wetting_kr(2:end,2)-wetting_kr(1:end-1,2))./(wetting_kr(2:end,1)-wetting_kr(1:end-1,1));
            dkro_Sg_value = [dkro_Sg_value;dkro_Sg_value(end)];
            obj.dkrw_Sg_table{faciesID} = griddedInterpolant(wetting_kr(:,1),dkro_Sg_value);
        end

        %% Build CapillaryPressure Interpolations
        function BuildCapillaryPressureInterpolations(obj, faciesID, Pc_drainage, Pc_imbibation)

            % capillary pressure
            % convert relations to functions of S_g;
            Pc_drainage(:,1) = 1.0 - Pc_drainage(:,1);
            Pc_imbibation(:,1) = 1.0 - Pc_imbibation(:,1);
            Pc_imbibation = flip(Pc_imbibation); % stored in ascending order
            Pc_drainage = flip(Pc_drainage); % stored in ascending order

            drainage_dPc_value = (Pc_drainage(2:end,2)-Pc_drainage(1:end-1,2))./(Pc_drainage(2:end,1)-Pc_drainage(1:end-1,1));
            drainage_dPc_value = [drainage_dPc_value(1);drainage_dPc_value];

            obj.Pc_drainage_Sg_table{faciesID} = griddedInterpolant(Pc_drainage(:,1),Pc_drainage(:,2));
            obj.dPc_drainage_Sg_table{faciesID} = griddedInterpolant(Pc_drainage(:,1),drainage_dPc_value);
            obj.Pc_imbibation_Sg_table{faciesID} = griddedInterpolant(Pc_imbibation(:,1),Pc_imbibation(:,2));

        end

        %% Build Scanning Curve Interpolations
        function BuildScanningCurveInterpolations(obj, faciesID, N_scanning_curves, N_scanning_points)

            % Compute scanning curves mapping kr(Sg, Sgi), dkr(Sg, Sgi)
            % relative permeability
            % S_gi = linspace(2e-3, 1-Sr(2), N_scanning_curves);
            S_gi = linspace(2e-3, 1-obj.Swr_Facies{faciesID}, N_scanning_curves);

            for i = 1:length(S_gi)
                % compute residual saturation corresponding to S_gi
                % S_gt = obj.S_gcr + (S_gi-obj.S_gcr)/(1+obj.C*(S_gi-obj.S_gcr)); % Land's model
                S_gt = S_gi(i) * obj.fr; % linear model

                % Saturations of scanning curves range from Sgt to Sgi
                % S_g = linspace(S_gt, S_gi(i), N_scanning_points)';

                % Extending the saturations of scanning curves (from 0 to Sgi)
                S_g = linspace(S_gt, S_gi(i), N_scanning_points - 2);
                S_g0 = linspace(0, S_gt, 3);
                S_g = [S_g0(1:end-1), S_g]';

                S_gstar = obj.S_gtmax{faciesID} + (S_g-S_gt).*(obj.S_gimax{faciesID}-obj.S_gtmax{faciesID})/(S_gi(i)-S_gt);
                k_rgi = obj.nonwetting_imbibation_table{faciesID}(S_gstar)*obj.krg_Sg_table{faciesID}(S_gi(i))/(obj.krg_Sg_table{faciesID}(obj.S_gimax{faciesID})); % Killough's model
                k_rgi = max(0, k_rgi);

                dk_rgi = (k_rgi(2:end)-k_rgi(1:end-1))./(S_g(2:end)-S_g(1:end-1));
                dk_rgi = [dk_rgi;dk_rgi(end)];

                obj.scanning_curve_Sg(i, :) = S_g';
                obj.scanning_curve_kr(i, :) = k_rgi';
                obj.scanning_curve_dkr(i, :) = dk_rgi';


                % capillary pressure
                epsilon = 0.1;
                F = (1./(1.0-S_g-(1.0-S_gi(i))+epsilon)-1.0/epsilon)./(1./((1.0-S_gt)-(1.0-S_gi(i))+epsilon)-1.0/epsilon); % the fraction of difference between Pc_im(Sg) and Pc_dr(Sg), 0<F<1
                P_c = obj.Pc_drainage_Sg_table{faciesID}(S_g) + F.*(obj.Pc_imbibation_Sg_table{faciesID}(S_g)-obj.Pc_drainage_Sg_table{faciesID}(S_g)); % Killough's model
                P_c = max(0, P_c);

                dP_c = (P_c(2:end)-P_c(1:end-1))./(S_g(2:end)-S_g(1:end-1));
                dP_c = [dP_c(1);dP_c];

                obj.scanning_curve_Pc(i,:) = P_c';
                obj.scanning_curve_dPc(i,:) = dP_c';
            end


            obj.scanning_curve_Sgi = zeros(size(obj.scanning_curve_Sg, 1), size(obj.scanning_curve_Sg, 2)) + S_gi';

            obj.scanning_curve_kr_table{faciesID} = scatteredInterpolant(obj.scanning_curve_Sg(:), obj.scanning_curve_Sgi(:), obj.scanning_curve_kr(:));
            obj.scanning_curve_dkr_table{faciesID} = scatteredInterpolant(obj.scanning_curve_Sg(:), obj.scanning_curve_Sgi(:), obj.scanning_curve_dkr(:));

            obj.scanning_curve_Pc_table{faciesID} = scatteredInterpolant(obj.scanning_curve_Sg(:), obj.scanning_curve_Sgi(:), obj.scanning_curve_Pc(:));
            obj.scanning_curve_dPc_table{faciesID} = scatteredInterpolant(obj.scanning_curve_Sg(:), obj.scanning_curve_Sgi(:), obj.scanning_curve_dPc(:));

        end

        %% Compute kr -- No Hysteresis
        function kr = ComputeRelPerm(obj, Phases, Reservoir)
            status_S_1 = Reservoir.State.Properties('S_1');
            S = status_S_1.Value;

            FaciesIndex = Reservoir.FaciesIndex;
            N_Facies = max(FaciesIndex);

            for faciesID = 1:N_Facies
                % safe-guard % make sure Srg < Sg < 1-Sro
                % S = max(S, Phases(1).sr);
                % S = min(S, 1-Phases(2).sr);
                Index = faciesID == FaciesIndex;
                % S(Index) = max(S(Index), obj.Sgr_Facies{faciesID});
                S(Index) = min(S(Index), 1-obj.Swr_Facies{faciesID});

                % phase 1 (gas phase)
                obj.kr(Index,1) = obj.krg_Sg_table{faciesID}(S(Index));
                % phase 2 (oil phase)
                obj.kr(Index,2) = obj.krw_Sg_table{faciesID}(S(Index));
            end
            kr = obj.kr;

            % % safe-guard
            % S = max(S, Phases(1).sr);
            % S = min(S, 1-Phases(2).sr);
            %
            % % phase 1 (gas phase)
            % kr(:,1) = obj.krg_Sg_table(S);
            % % phase 2 (oil phase)
            % kr(:,2) = obj.krw_Sg_table(S);

            status_S_1.Value = S;
            status_S_2 = Reservoir.State.Properties('S_2');
            status_S_2.Value = 1.0-S;

        end

        %% Compute dkr -- No Hysteresis
        function dkr = ComputeDerivative(obj, Phases, Reservoir)
            status_S_1 = Reservoir.State.Properties('S_1');
            S = status_S_1.Value;

            FaciesIndex = Reservoir.FaciesIndex;
            N_Facies = max(FaciesIndex);

            for faciesID = 1:N_Facies
                % safe-guard % make sure Srg < Sg < 1-Sro
                % S = max(S, Phases(1).sr);
                % S = min(S, 1-Phases(2).sr);
                Index = faciesID == FaciesIndex;
                % S(Index) = max(S(Index), obj.Sgr_Facies{faciesID});
                S(Index) = min(S(Index), 1-obj.Swr_Facies{faciesID});

                % phase 1 (gas phase)
                obj.dkr(Index,1) = obj.dkrg_Sg_table{faciesID}(S(Index));
                % phase 2 (oil phase)
                obj.dkr(Index,2) = obj.dkrw_Sg_table{faciesID}(S(Index));
            end
            dkr = obj.dkr;

            % % safe-guard
            % S = max(S, Phases(1).sr);
            % S = min(S, 1-Phases(2).sr);
            %
            % % phase 1 (gas phase)
            % dkr(:,1) = obj.dkrg_Sg_table(S);
            % % phase 2 (oil phase)
            % dkr(:,2) = obj.dkrw_Sg_table(S);

        end

        %% Compute kr -- With Hysteresis
        function kr = ComputeRelPermHysteresis(obj, Phases, Reservoir) % dkr is also computted in this function
            status_S_1 = Reservoir.State.Properties('S_1');
            S = Reservoir.State.Properties('S_1').Value;
            % S_old = Reservoir.State_old.Properties('S_1').Value;

            FaciesIndex = Reservoir.FaciesIndex;
            N_Facies = max(FaciesIndex);

            for faciesID = 1:N_Facies
                % safe-guard
                % S = max(S, Phases(1).sr);
                % S = min(S, 1-Phases(2).sr); % make sure Srg < Sg < 1-Sro
                Index = faciesID == FaciesIndex;
                S(Index) = min(S(Index), 1-obj.Swr_Facies{faciesID}); % make sure Srg < Sg < 1-Sro

                % phase 2 (water phase, no hysteresis)
                obj.kr(Index,2) = obj.krw_Sg_table{faciesID}(S(Index));
                obj.dkr(Index,2) = obj.dkrw_Sg_table{faciesID}(S(Index)); % dkr

                % phase 1 (gas phase)
                % no hysteresis
                % obj.kr(:,1) = obj.krg_Sg_table(S);
                % obj.dkr(:,1) = obj.dkrg_Sg_table(S); % dkr


                % with hysteresis
                % on primary drainage
                Index = and(faciesID == FaciesIndex, obj.process_label == 0);
                obj.kr(Index, 1) = obj.krg_Sg_table{faciesID}(S(Index));
                obj.dkr(Index, 1) = obj.dkrg_Sg_table{faciesID}(S(Index)); % dkr

                % on Scanning curves
                % S(obj.process_label == 1) = max(S(obj.process_label == 1), obj.Sg_r(obj.process_label == 1)); % safe-guard
                Index = and(faciesID == FaciesIndex, obj.process_label == 1);
                S(Index) = min(S(Index), obj.reverse_point(Index)); % safe-guard
                obj.kr(Index, 1) = obj.scanning_curve_kr_table{faciesID}(S(Index), obj.reverse_point(Index));
                obj.dkr(Index, 1) = obj.scanning_curve_dkr_table{faciesID}(S(Index), obj.reverse_point(Index)); % dkr
            end
            kr = obj.kr;

            % update saturation and pressure

            % obj.safeguard_S_label(:) = 0; % reset label
            % obj.safeguard_S_label(status_S_1.Value ~= S) = 1;
            % status_safeguard_label = Reservoir.State.Properties('safeguard_label');
            % status_safeguard_label.Value = obj.safeguard_S_label;

            status_S_1.Value = S;
            status_S_2 = Reservoir.State.Properties('S_2');
            status_S_2.Value = 1.0-S;

            obj.immobile_label = 0;
            obj.immobile_label = (obj.process_label == 1).*(S <= obj.Sg_r);
            status_immobile_label = Reservoir.State.Properties('immobile_label');
            status_immobile_label.Value = obj.immobile_label;

        end

        %% Compute dkr -- With Hysteresis
        function dkr = ComputeDerivativeHysteresis(obj)
            dkr = obj.dkr;
        end

        %% Compute TransportLabel
        function ComputeTransportLabel(obj, Phases, Reservoir)

            S_old = Reservoir.State_old.Properties('S_1').Value;
            S_old_old = Reservoir.State_old_old.Properties('S_1').Value;

            FaciesIndex = Reservoir.FaciesIndex;
            N_Facies = max(FaciesIndex);

            % Transition from drainage bounding curve to scanning curve
            drainage_to_scanning_label = (obj.process_label == 0) .* (S_old < S_old_old) .* (abs(S_old - S_old_old) > obj.ReversalTolerance) .* (S_old_old > obj.MinReversalSaturation );
            obj.process_label(drainage_to_scanning_label == 1) = 1;
            obj.reverse_point(drainage_to_scanning_label == 1) = S_old_old(drainage_to_scanning_label == 1); % update reverse point
            obj.Sg_r(drainage_to_scanning_label == 1) = obj.fr * S_old_old(drainage_to_scanning_label == 1); % update residual saturation

            % Transition from scanning curve back to drainage bounding curve
            scanning_to_drainage_label = (obj.process_label == 1) .* (S_old > S_old_old) .* (S_old >= obj.reverse_point);
            obj.process_label(scanning_to_drainage_label == 1) = 0;
            % obj.reverse_point(scanning_to_drainage_label == 1) = 1 - Phases(2).sr; % update reverse point
            % obj.Sg_r(scanning_to_drainage_label == 1) = Phases(1).sr; % update residual saturation
            for faciesID = 1:N_Facies
                Index = and(faciesID == FaciesIndex, scanning_to_drainage_label == 1);
                obj.reverse_point(Index) = 1 - obj.Swr_Facies{faciesID};
                obj.Sg_r(Index) = obj.Sgr_Facies{faciesID};
            end

            % update Sg_r when S is smaller than Sg_r
            obj.Sg_r(obj.process_label == 1) = min(obj.Sg_r(obj.process_label == 1), Reservoir.State.Properties('S_1').Value(obj.process_label == 1));

            status_process_label = Reservoir.State.Properties('process_label'); % update process_label
            status_process_label.Value = obj.process_label;
            status_reverse_point = Reservoir.State.Properties('reverse_point'); % update reverse_point
            status_reverse_point.Value = obj.reverse_point;
            status_residual_saturation = Reservoir.State.Properties('residual_saturation'); % update residual saturation
            status_residual_saturation.Value = obj.Sg_r;
        end

        %% reset relative permeability model
        function ResetRelPermModel(obj, State)
            obj.process_label = State.Properties('process_label').Value;
            obj.reverse_point = State.Properties('reverse_point').Value;
            obj.Sg_r = State.Properties('residual_saturation').Value;
        end
    end
end



