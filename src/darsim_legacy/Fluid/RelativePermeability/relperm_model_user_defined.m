 % Table relative permeability model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Yuhang Wang
%TU Delft
%Created: 6 Jan 2021
%Last modified: 14 July 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef relperm_model_user_defined < relperm_model
    properties
        krg_Sg_table
        nonwetting_imbibation_table
        kro_Sg_table
        dkrg_Sg_table
        dkro_Sg_table
        
        Pc_drainage_Sg_table
        dPc_drainage_Sg_table
        Pc_imbibation_Sg_table
        dPc_imbibation_Sg_table
        
        scanning_curve_kr_table
        scanning_curve_dkr_table
        scanning_curve_Pc_table
        scanning_curve_dPc_table
        
        kr
        dkr

        scanning_curve_Sg
        scanning_curve_kr
        scanning_curve_dkr
        scanning_curve_Pc
        scanning_curve_dPc
        scanning_curve_Sgi
        
        reverse_point
        process_label % 0: primary drainage; 1: scanning curve
        % immobile_label; % 0; mobile; 1: immobile
        Sg_r % residual saturation of gas phase
        Sw_r0
        Sg_r0
        fr0
        
        C % Land trapping coefficient
        S_gcr % critical gas phase saturation
        S_gimax % maximum gas phase saturation
        S_gtmax % maximum residual saturation of gas phase
        
        
        immobile_label
        SandType_index
        
        n = [2, 2]
   
        swr = [0.10, 0.32, 0.14, 0.12, 0.12, 0.12, 0.10, 0]
        sgr = [0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10]
        
        pce = [0, 1500, 300, 100, 25, 10, 1, 0]

          
    end
    methods
        
        function Initialise(obj, N)

            obj.process_label = zeros(N,1);
            obj.immobile_label = zeros(N,1);
            obj.reverse_point = zeros(N,1);
            obj.Sg_r = zeros(N,1);
            obj.Sg_r0 = zeros(N,1);
            obj.Sw_r0 = zeros(N,1);
            obj.fr0 = zeros(N,1);
            
            for S = 1 : 8   
                obj.reverse_point(obj.SandType_index == S) = (1 - obj.swr(S));
                obj.Sg_r(obj.SandType_index == S) = obj.sgr(S);
                % obj.Sg_r0(obj.SandType_index == S) = obj.sgr(S);
                obj.Sw_r0(obj.SandType_index == S) = obj.swr(S);
                % obj.fr0(obj.SandType_index == S) = obj.fr(S);
            end
        end
        
        function obj = relperm_model_table_old(SimulationInput, Sr)
            % allocate memory
            N = SimulationInput.ReservoirProperties.Grid.N_ActiveCells;
            
            obj.SandType_index = SimulationInput.ReservoirProperties.FaciesIndex;
            N_surfaces = 8;
            N_scanning_curves = 10;
            N_scanning_points = 25;
            
            obj.kr = zeros(N,2);
            obj.dkr = zeros(N,2);

            obj.scanning_curve_Sg  = zeros(N_scanning_curves, N_scanning_points);
            obj.scanning_curve_kr  = zeros(N_scanning_curves, N_scanning_points);
            obj.scanning_curve_dkr = zeros(N_scanning_curves, N_scanning_points);
            obj.scanning_curve_Pc  = zeros(N_scanning_curves, N_scanning_points);
            obj.scanning_curve_dPc = zeros(N_scanning_curves, N_scanning_points);
 
            obj.reverse_point = zeros(N,1);
            obj.process_label = zeros(N,1);
            obj.immobile_label = zeros(N,1);
            obj.Sg_r = zeros(N,1);
            
            obj.S_gimax = 1-Sr(2);
            obj.S_gtmax = obj.S_gimax * 0.4; % linear model
            
            
            for S = 1 : N_surfaces
                % relations are functions of resident wetting phase;
                % read from LookupTable
                wetting_kr=[];
                nonwetting_drainage_kr=[];
                nonwetting_imbibation_kr=[];
                Pc_drainage=[];
                Pc_imbibation=[];

                currentfolder=pwd;
                cd(SimulationInput.Directory);
                
                cd(fullfile(SimulationInput.Directory, 'LookupTable'));
                files = dir(fullfile(pwd,'*.txt*'));
                fid = fopen(fullfile(files(S).folder,files(S).name),'r');
                
                while true
                    thisline = fgetl(fid);
                    if (strcmp(thisline, 'end')) break; end

                    if thisline(1) == '-'
                        keyword = thisline(2:end);
                        if (strcmp(keyword, 'wetting_kr'))
                            thisline = fgetl(fid);
                            while ~isempty(thisline)
                                wetting_kr =[wetting_kr; sscanf(thisline,'%f',[1 2])];
                                thisline = fgetl(fid);
                            end
                        end
                        if (strcmp(keyword, 'nonwetting_drainage_kr'))
                            thisline = fgetl(fid);
                            while ~isempty(thisline)
                                nonwetting_drainage_kr =[nonwetting_drainage_kr; sscanf(thisline,'%f',[1 2])];
                                thisline = fgetl(fid);
                            end
                        end
                        if (strcmp(keyword, 'nonwetting_imbibation_kr'))
                            thisline = fgetl(fid);
                            while ~isempty(thisline)
                                nonwetting_imbibation_kr =[nonwetting_imbibation_kr; sscanf(thisline,'%f',[1 2])];
                                thisline = fgetl(fid);
                            end
                        end
                        if (strcmp(keyword, 'Pc_drainage'))
                            thisline = fgetl(fid);
                            while ~isempty(thisline)
                                Pc_drainage =[Pc_drainage; sscanf(thisline,'%f',[1 2])];
                                thisline = fgetl(fid);
                            end
                        end
                        if (strcmp(keyword, 'Pc_imbibation'))
                            thisline = fgetl(fid);
                            while ~isempty(thisline)
                                Pc_imbibation =[Pc_imbibation; sscanf(thisline,'%f',[1 2])];
                                thisline = fgetl(fid);
                            end
                        end
                    end

                end
                fclose(fid);
                cd(currentfolder);

                %% relative permeability
                % convert relations to functions of S_g;
                wetting_kr(:,1) = 1.0 - wetting_kr(:,1);
                nonwetting_drainage_kr(:,1) = 1.0 - nonwetting_drainage_kr(:,1);
                nonwetting_imbibation_kr(:,1) = 1.0 - nonwetting_imbibation_kr(:,1);

                wetting_kr = flip(wetting_kr); % stored in ascending order
                nonwetting_imbibation_kr = flip(nonwetting_imbibation_kr); % stored in ascending order
                nonwetting_drainage_kr = flip(nonwetting_drainage_kr); % stored in ascending order

                obj.nonwetting_imbibation_table{S} = griddedInterpolant(nonwetting_imbibation_kr(:,1),nonwetting_imbibation_kr(:,2));                  
                obj.krg_Sg_table{S} = griddedInterpolant(nonwetting_drainage_kr(:,1),nonwetting_drainage_kr(:,2));
                obj.kro_Sg_table{S}  = griddedInterpolant(wetting_kr(:,1),wetting_kr(:,2));

                dkrg_Sg_value = (nonwetting_drainage_kr(2:end,2)-nonwetting_drainage_kr(1:end-1,2))./(nonwetting_drainage_kr(2:end,1)-nonwetting_drainage_kr(1:end-1,1));
                dkrg_Sg_value = [dkrg_Sg_value(1);dkrg_Sg_value];
                obj.dkrg_Sg_table{S} = griddedInterpolant(nonwetting_drainage_kr(:,1),dkrg_Sg_value);
                dkro_Sg_value = (wetting_kr(2:end,2)-wetting_kr(1:end-1,2))./(wetting_kr(2:end,1)-wetting_kr(1:end-1,1));
                dkro_Sg_value = [dkro_Sg_value;dkro_Sg_value(end)];
                obj.dkro_Sg_table{S} = griddedInterpolant(wetting_kr(:,1),dkro_Sg_value);        

                %% capillary pressure
                % convert relations to functions of S_g;
                Pc_drainage(:,1) = 1.0 - Pc_drainage(:,1);
                Pc_imbibation(:,1) = 1.0 - Pc_imbibation(:,1);
                Pc_imbibation = flip(Pc_imbibation); % stored in ascending order
                Pc_drainage = flip(Pc_drainage); % stored in ascending order

                drainage_dPc_value = (Pc_drainage(2:end,2)-Pc_drainage(1:end-1,2))./(Pc_drainage(2:end,1)-Pc_drainage(1:end-1,1));
                drainage_dPc_value = [drainage_dPc_value(1);drainage_dPc_value];

                obj.Pc_drainage_Sg_table{S}  = griddedInterpolant(Pc_drainage(:,1),Pc_drainage(:,2));
                obj.dPc_drainage_Sg_table{S}  = griddedInterpolant(Pc_drainage(:,1),drainage_dPc_value);
                obj.Pc_imbibation_Sg_table{S}  = griddedInterpolant(Pc_imbibation(:,1),Pc_imbibation(:,2));

                %% Compute scanning curves mapping kr(Sg, Sgi), dkr(Sg, Sgi)
                % relative permeability
                S_gi = linspace(2e-3, 1-Sr(2), N_scanning_curves);

                for i = 1:length(S_gi)
                    % compute residual saturation corresponding to S_gi
    %                 S_gt = obj.S_gcr + (S_gi-obj.S_gcr)/(1+obj.C*(S_gi-obj.S_gcr)); % Land's model
                    S_gt = S_gi(i) * 0.4; % linear model

                    % Saturations of scanning curves range from Sgt to Sgi
    %                 S_g = linspace(S_gt, S_gi(i), N_scanning_points)';

                    % Extending the saturations of scanning curves (from 0 to Sgi)
                    S_g = linspace(S_gt, S_gi(i), N_scanning_points - 5);
                    S_g0 = linspace(0, S_gt, 6);
                    S_g = [S_g0(1:end-1), S_g]';

                    S_gstar = obj.S_gtmax + (S_g-S_gt).*(obj.S_gimax-obj.S_gtmax)/(S_gi(i)-S_gt);
                    k_rgi = obj.nonwetting_imbibation_table{S}(S_gstar)*obj.krg_Sg_table{S}(S_gi(i))/(obj.krg_Sg_table{S}(obj.S_gimax)); % Killough's model
                    k_rgi = max(0, k_rgi);

                    dk_rgi = (k_rgi(2:end)-k_rgi(1:end-1))./(S_g(2:end)-S_g(1:end-1)); 
                    dk_rgi = [dk_rgi;dk_rgi(end)];

                    obj.scanning_curve_Sg(i, :) = S_g';
                    obj.scanning_curve_kr(i, :) = k_rgi';
                    obj.scanning_curve_dkr(i, :) = dk_rgi';


                    %% capillary pressure
                    epsilon = 0.1;
                    F = (1./(1.0-S_g-(1.0-S_gi(i))+epsilon)-1.0/epsilon)./(1./((1.0-S_gt)-(1.0-S_gi(i))+epsilon)-1.0/epsilon); % the fraction of difference between Pc_im(Sg) and Pc_dr(Sg), 0<F<1
                    P_c = obj.Pc_drainage_Sg_table{S}(S_g) + F.*(obj.Pc_imbibation_Sg_table{S}(S_g)-obj.Pc_drainage_Sg_table{S}(S_g)); % Killough's model
                    P_c = max(0, P_c);

                    dP_c = (P_c(2:end)-P_c(1:end-1))./(S_g(2:end)-S_g(1:end-1));
                    dP_c = [dP_c(1);dP_c];

                    obj.scanning_curve_Pc(i, :) = P_c';
                    obj.scanning_curve_dPc(i, :) = dP_c';
                end

                obj.scanning_curve_Sgi = zeros(size(obj.scanning_curve_Sg, 1), size(obj.scanning_curve_Sg, 2)) + S_gi';

                obj.scanning_curve_kr_table{S} = scatteredInterpolant(obj.scanning_curve_Sg(:), obj.scanning_curve_Sgi(:), obj.scanning_curve_kr(:));
                obj.scanning_curve_dkr_table{S} = scatteredInterpolant(obj.scanning_curve_Sg(:), obj.scanning_curve_Sgi(:), obj.scanning_curve_dkr(:));

                obj.scanning_curve_Pc_table{S} = scatteredInterpolant(obj.scanning_curve_Sg(:), obj.scanning_curve_Sgi(:), obj.scanning_curve_Pc(:));
                obj.scanning_curve_dPc_table{S} = scatteredInterpolant(obj.scanning_curve_Sg(:), obj.scanning_curve_Sgi(:), obj.scanning_curve_dPc(:));
            end
        end
        
        function obj = relperm_model_user_defined(SimulationInput, Sr)
            % allocate memory
            N = SimulationInput.ReservoirProperties.Grid.N_ActiveCells;
            
            N_surfaces = 8;
            N_scanning_curves = 1;
            N_scanning_points = 35;
            
            obj.kr = zeros(N,2);
            obj.dkr = zeros(N,2);

            obj.scanning_curve_Sg  = zeros(N_scanning_curves, N_scanning_points);
            obj.scanning_curve_kr  = zeros(N_scanning_curves, N_scanning_points);
            obj.scanning_curve_dkr = zeros(N_scanning_curves, N_scanning_points);
            obj.scanning_curve_Pc  = zeros(N_scanning_curves, N_scanning_points);
            obj.scanning_curve_dPc = zeros(N_scanning_curves, N_scanning_points);

            obj.Initialise(N);  
            
            obj.S_gimax = 1- [0.10, 0.32, 0.14, 0.12, 0.12, 0.12, 0.10, 0];
            
            penmax = 2500; % Maximum capillary pressure, adjust as needed
            
            for S = 1 : N_surfaces
                
                % obj.Sg_r(obj.SandType_index == S) = obj.Sg_r(obj.SandType_index == S) + obj.sgr(S);

                %% relative permeability
                % wetting and non-wetting drainage
                %sw = (exp(flip(linspace(0, 5.0, 1000))) - 1) /(exp(5.0) - 1);
                %sg = linspace(0, 1, 1000);%1-sw;
                %sw = 1-sg;
                %sg = 1 - sw;
                sg = linspace(0, 1, 1000);
                sw = 1 - sg;
                wetting_kr = [(sg); (max(0, (sw - obj.swr(S)) / (1 - obj.swr(S)))).^2]'; % Updated wetting rel perm
                nonwetting_drainage_kr = [(sg); (max(0, (sg - obj.sgr(S)) / (1 - obj.sgr(S)))).^2]'; % Updated non-wetting rel perm

                obj.krg_Sg_table{S} = griddedInterpolant(nonwetting_drainage_kr(:,1),nonwetting_drainage_kr(:,2));
                obj.kro_Sg_table{S}  = griddedInterpolant(wetting_kr(:,1),wetting_kr(:,2));
            
                dkrg_Sg_value = (nonwetting_drainage_kr(2:end,2)-nonwetting_drainage_kr(1:end-1,2))./(nonwetting_drainage_kr(2:end,1)-nonwetting_drainage_kr(1:end-1,1));
                dkrg_Sg_value = [dkrg_Sg_value(1);dkrg_Sg_value];
                obj.dkrg_Sg_table{S} = griddedInterpolant(nonwetting_drainage_kr(:,1),dkrg_Sg_value);
                dkro_Sg_value = (wetting_kr(2:end,2)-wetting_kr(1:end-1,2))./(wetting_kr(2:end,1)-wetting_kr(1:end-1,1));
                dkro_Sg_value = [dkro_Sg_value;dkro_Sg_value(end)];
                obj.dkro_Sg_table{S} = griddedInterpolant(wetting_kr(:,1),dkro_Sg_value);        

                %% capillary pressure

                sg = linspace(0,1,1000);
                Pc_drainage(:, 1) = sg;
 
                if any(S == [2, 3, 4, 5, 6, 7])
                    swe = (1 - sg - obj.swr(S))/(1 - obj.swr(S));
                    swe = max(swe,0.01);
                    swe = min(swe,1);
                    Pc_drainage(:,2) = penmax*erf(obj.pce(S) * swe.^(-1/2) * sqrt(pi) / (penmax*2));
                else
                    Pc_drainage(:,2) = 0;
                end
            
                drainage_dPc_value = (Pc_drainage(2:end,2)-Pc_drainage(1:end-1,2))./(Pc_drainage(2:end,1)-Pc_drainage(1:end-1,1));
                drainage_dPc_value = [drainage_dPc_value(1);drainage_dPc_value];

                obj.Pc_drainage_Sg_table{S}  = griddedInterpolant(Pc_drainage(:,1),Pc_drainage(:,2));
                obj.dPc_drainage_Sg_table{S}  = griddedInterpolant(Pc_drainage(:,1),drainage_dPc_value);

            end
        end      

              
        function kr = ComputeRelPerm(obj, Phases, Medium)
            status_S_1 = Medium.State.Properties('S_1');
            S = Medium.State.Properties('S_1').Value;
            % safe-guard
            S = max(S, Phases(1).sr);
            S = min(S, 1-Phases(2).sr);
%             % phase 1 (gas phase) 
%             kr(:,1) = obj.krg_Sg_table(S);
%             % phase 2 (oil phase) 
%             kr(:,2) = obj.kro_Sg_table(S);
 
            % phase 1 (gas phase) 
            kr(:,1) = S/(1-Phases(2).sr);
            % phase 2 (oil phase) 
            kr(:,2) = 1.0 - S/(1-Phases(2).sr);           
  
            status_S_1.Value = S; 
            status_S_2 = Medium.State.Properties('S_2');
            status_S_2.Value = 1.0-S;
            
        end
        
        function dkr = ComputeDerivative(obj, Phases, S)
            % safe-guard
            S = max(S, Phases(1).sr);
            S = min(S, 1-Phases(2).sr);
%             % phase 1 (gas phase) 
%             dkr(:,1) = obj.dkrg_Sg_table(S);
%             % phase 2 (oil phase) 
%             dkr(:,2) = obj.dkro_Sg_table(S);
            
            % phase 1 (gas phase)
            dkr(:,1) = zeros(length(S),1) + 1/(1-Phases(2).sr);
            % phase 2 (oil phase) 
            dkr(:,2) = zeros(length(S),1) - 1/(1-Phases(2).sr);
        end 
        
        function kr = ComputeRelPermHysteresis_old(obj, Phases, Reservoir) % dkr is also computted in this function
            status_S_1 = Reservoir.State.Properties('S_1');
            S = Reservoir.State.Properties('S_1').Value;
            % S_old = Reservoir.State_old.Properties('S_1').Value;
            
            for i = 1:8
                % safe-guard
                % S = max(S, Phases(1).sr); 
                % S = min(S, 1-Phases(2).sr); % make sure Srg < Sg < 1-Sro

                % phase 2 (oil phase, no hysteresis)
                sand_label = obj.SandType_index == i;
                obj.kr(obj.SandType_index == i,2) = obj.kro_Sg_table{i}(S(sand_label));
                obj.dkr(obj.SandType_index == i,2) = obj.dkro_Sg_table{i}(S(sand_label)); % dkr
        
                % phase 1 (gas phase)       
                % no hysteresis
                % obj.kr(:,1) = obj.krg_Sg_table(S);
                % obj.dkr(:,1) = obj.dkrg_Sg_table(S); % dkr


                % with hysteresis
                % on primary drainage
                sand_dra_label = logical((obj.SandType_index == i) .* (obj.process_label == 0));
                obj.kr(sand_dra_label, 1) = obj.krg_Sg_table{i}(S(sand_dra_label));
                obj.dkr(sand_dra_label, 1) = obj.dkrg_Sg_table{i}(S(sand_dra_label)); % dkr

                % on Scanning curves
                % S(obj.process_label == 1) = max(S(obj.process_label == 1), obj.Sg_r(obj.process_label == 1)); % safe-guard
                sand_imb_label = logical((obj.SandType_index == i) .* (obj.process_label == 1));
                % S(sand_imb_label) = min(S(sand_imb_label), obj.reverse_point(sand_imb_label)); % safe-guard
                obj.kr(sand_imb_label, 1) = obj.scanning_curve_kr_table{i}(S(sand_imb_label), obj.reverse_point(sand_imb_label));
                obj.dkr(sand_imb_label, 1) = obj.scanning_curve_dkr_table{i}(S(sand_imb_label), obj.reverse_point(sand_imb_label)); % dkr

            end
            
            kr = obj.kr;
            
            kr = min(kr, 0.71);
            kr = max(kr, 0);
            
            
            % update saturation and pressure         
            status_S_1.Value = S; 
            status_S_2 = Reservoir.State.Properties('S_2');
            status_S_2.Value = 1.0-S;
            
            obj.immobile_label = 0;
            obj.immobile_label = (obj.process_label == 1).*(S <= obj.Sg_r);
            status_immobile_label = Reservoir.State.Properties('immobile_label');
            status_immobile_label.Value = obj.immobile_label;
            

        end
        function kr = ComputeRelPermHysteresis(obj, Phases, Reservoir) % dkr is also computted in this function
            status_S_1 = Reservoir.State.Properties('S_1');
            S = Reservoir.State.Properties('S_1').Value;

            for i = 1:8

                % Phase 2 relative permeability
                sand_label = obj.SandType_index == i;
                % safe-guard
                % S = max(S, Phases(1).sr);
                % S(sand_label == 1) = min(S(sand_label == 1), 1 - obj.Sw_r0(sand_label == 1)); % make sure Srg < Sg < 1-Sr
                obj.kr(sand_label == 1,2) = obj.kro_Sg_table{i}(S(sand_label == 1));
                obj.dkr(sand_label == 1,2) = obj.dkro_Sg_table{i}(S(sand_label == 1)); % dkr
                obj.kr(sand_label == 1, 1) = obj.krg_Sg_table{i}(S(sand_label == 1));
                obj.dkr(sand_label == 1, 1) = obj.dkrg_Sg_table{i}(S(sand_label == 1));
                %                 obj.kr(logical((obj.SandType_index == i) .* (S < obj.sgr(i))), 2) = 1;
                %                 obj.kr(logical((obj.SandType_index == i) .* (S > 1 - obj.swr(i))), 2) = 0;
                %                 obj.dkr(logical((obj.SandType_index == i) .* (S < obj.sgr(i))), 2) = 0;
                %                 obj.dkr(logical((obj.SandType_index == i) .* (S > 1 - obj.swr(i))), 2) = 0;

                % Phase 1 relative permeability
                % with hysteresis
                % on primary drainage
%                 sand_dra_label = logical((obj.SandType_index == i) .* (obj.process_label == 0));
%                 obj.kr(sand_dra_label == 1, 1) = obj.krg_Sg_table{i}(S(sand_dra_label == 1));
%                 obj.dkr(sand_dra_label == 1, 1) = obj.dkrg_Sg_table{i}(S(sand_dra_label == 1)); % dkr
% 
%                 % on Scanning curves
%                 sand_imb_label = logical((obj.SandType_index == i) .* (obj.process_label == 1));
%                 % S(obj.process_label == 1) = max(S(obj.process_label == 1), obj.Sg_r(obj.process_label == 1)); % safe-guard
%                 S(sand_imb_label == 1) = min(S(sand_imb_label == 1), obj.reverse_point(sand_imb_label == 1)); % safe-guard
% 
%                 obj.kr(sand_imb_label == 1, 1) = obj.scanning_curve_kr_table{i}(S(sand_imb_label == 1));
%                 obj.dkr(sand_imb_label == 1, 1) = obj.scanning_curve_dkr_table{i}(S(sand_imb_label == 1)); % dkr

                %                 obj.kr(logical((obj.SandType_index == i) .* (S < obj.sgr(i))), 1) = 0;
                %                 obj.kr(logical((obj.SandType_index == i) .* (S > 1 - obj.swr(i))), 1) = 1;
                %                 obj.dkr(logical((obj.SandType_index == i) .* (S < obj.sgr(i))), 1) = 0;
                %                 obj.dkr(logical((obj.SandType_index == i) .* (S > 1 - obj.swr(i))), 1) = 0;

            end

            kr = obj.kr;
            
            % update saturation and pressure         
            status_S_1.Value = S; 
            status_S_2 = Reservoir.State.Properties('S_2');
            status_S_2.Value = 1.0-S;
            
            obj.immobile_label = 0;
            obj.immobile_label = (obj.process_label == 1).*(S <= obj.Sg_r);
            status_immobile_label = Reservoir.State.Properties('immobile_label');
            status_immobile_label.Value = obj.immobile_label;
            

        end  
        function kr = ComputeRelPermHysteresis_old_old(obj, Phases, Reservoir) % dkr is also computted in this function
            status_S_1 = Reservoir.State.Properties('S_1');
            S = Reservoir.State.Properties('S_1').Value;
            S_old = Reservoir.State_old.Properties('S_1').Value;
            for i = 1:8
                
                % Phase 2 relative permeability
                sand_label = obj.SandType_index == i;       
                % safe-guard
                % S = max(S, Phases(1).sr); 
                S(sand_label == 1) = min(S(sand_label == 1), 1 - obj.Sw_r0(sand_label == 1)); % make sure Srg < Sg < 1-Sro
                
                obj.kr(sand_label == 1,2) = obj.kro_Sg_table{i}(S(sand_label == 1));
                obj.dkr(sand_label == 1,2) = obj.dkro_Sg_table{i}(S(sand_label == 1)); % dkr
                obj.kr(sand_label == 1,1) = obj.krg_Sg_table{i}(S(sand_label == 1));
                obj.dkr(sand_label == 1,1) = obj.dkrg_Sg_table{i}(S(sand_label == 1));
                 
%                 obj.kr(logical((obj.SandType_index == i) .* (S < obj.sgr(i))), 2) = 1;
%                 obj.kr(logical((obj.SandType_index == i) .* (S > 1 - obj.swr(i))), 2) = 0;
%                 obj.dkr(logical((obj.SandType_index == i) .* (S < obj.sgr(i))), 2) = 0;
%                 obj.dkr(logical((obj.SandType_index == i) .* (S > 1 - obj.swr(i))), 2) = 0;

                % Phase 1 relative permeability
                % with hysteresis
                % on primary drainage
%                 sand_dra_label = logical((obj.SandType_index == i) .* (obj.process_label == 0));
%                 obj.kr(sand_dra_label == 1, 1) = obj.krg_Sg_table{i}(S(sand_dra_label == 1));
%                 obj.dkr(sand_dra_label == 1, 1) = obj.dkrg_Sg_table{i}(S(sand_dra_label == 1)); % dkr
%                 
%                 % on Scanning curves                
%                 sand_imb_label = logical((obj.SandType_index == i) .* (obj.process_label == 1));                               
%                 % S(obj.process_label == 1) = max(S(obj.process_label == 1), obj.Sg_r(obj.process_label == 1)); % safe-guard
%                % S(sand_imb_label == 1) = min(S(sand_imb_label == 1), obj.reverse_point(sand_imb_label == 1)); % safe-guard
%                 
% %                 transition_below_Sgr = logical((sand_imb_label==1) .* (S_old>=obj.Sg_r) .* (S<=obj.Sg_r));
% %                 consistent_below_Sgr = logical((sand_imb_label==1) .* (S_old<=obj.Sg_r) .* (S<=obj.Sg_r));
% %                 S(transition_below_Sgr == 1) = obj.Sg_r(transition_below_Sgr==1);
% %                 S(consistent_below_Sgr == 1) = S_old(consistent_below_Sgr == 1);  
% 
%                 obj.kr(sand_imb_label == 1, 1) = obj.scanning_curve_kr_table{i}(S(sand_imb_label == 1));
%                 obj.dkr(sand_imb_label == 1, 1) = obj.scanning_curve_dkr_table{i}(S(sand_imb_label == 1)); % dkr
%                 
%                 obj.kr(logical((obj.SandType_index == i) .* (S < obj.sgr(i))), 1) = 0;
%                 obj.kr(logical((obj.SandType_index == i) .* (S > 1 - obj.swr(i))), 1) = 1;
%                 obj.dkr(logical((obj.SandType_index == i) .* (S < obj.sgr(i))), 1) = 0;
%                 obj.dkr(logical((obj.SandType_index == i) .* (S > 1 - obj.swr(i))), 1) = 0; 
            end            
            kr = obj.kr;          
            % update saturation and pressure         
            status_S_1.Value = S; 
            status_S_2 = Reservoir.State.Properties('S_2');
            status_S_2.Value = 1.0-S;
            
            obj.immobile_label = 0;
            obj.immobile_label = (obj.process_label == 1).*(S <= obj.Sg_r);
            status_immobile_label = Reservoir.State.Properties('immobile_label');
            status_immobile_label.Value = obj.immobile_label;
           
        end  
 
        function dkr = ComputeDerivativeHysteresis(obj)
            dkr = obj.dkr;
        end
  
        %% Vectorizing for loops in function ComputeTransportLabel
        function ComputeTransportLabel(obj, Phases, Reservoir)
            
            S_old = Reservoir.State_old.Properties('S_1').Value;
            S_old_old = Reservoir.State_old_old.Properties('S_1').Value;
                        
            % Transition from drainage bounding curve to scanning curve
            drainage_to_scanning_label = (obj.process_label == 0) .* (S_old < S_old_old) .* (abs(S_old - S_old_old) > 1.0e-4) .* (S_old_old > 2e-3);
            obj.process_label(drainage_to_scanning_label == 1) = 1;
            obj.reverse_point(drainage_to_scanning_label == 1) = S_old_old(drainage_to_scanning_label == 1); % update reverse point
            obj.Sg_r(drainage_to_scanning_label == 1) = 0.1;%obj.fr0(drainage_to_scanning_label == 1) .* S_old_old(drainage_to_scanning_label == 1); % update residual saturation

            % Transition from scanning curve back to drainage bounding curve
            scanning_to_drainage_label = (obj.process_label == 1) .* (S_old > S_old_old) .* (S_old == obj.reverse_point);
            obj.process_label(scanning_to_drainage_label == 1) = 0;
            obj.reverse_point(scanning_to_drainage_label == 1) = 1 - obj.Sw_r0(scanning_to_drainage_label == 1); % update reverse point
            obj.Sg_r(scanning_to_drainage_label == 1) = 0.1; % update residual saturation
          
            % update Sg_r when S is smaller than Sg_r
            obj.Sg_r(obj.process_label == 1) = min(obj.Sg_r(obj.process_label == 1), Reservoir.State.Properties('S_1').Value(obj.process_label == 1));
             
            status_process_label = Reservoir.State.Properties('process_label'); % update process_label
            status_process_label.Value = obj.process_label;
            status_reverse_point = Reservoir.State.Properties('reverse_point'); % update reverse_point
            status_reverse_point.Value = obj.reverse_point;
            status_residual_saturation = Reservoir.State.Properties('residual_saturation'); % update residual saturation
            status_residual_saturation.Value = obj.Sg_r;
        end

        % reset relative permeability model
        function ResetRelPermModel(obj, State)
            obj.process_label = State.Properties('process_label').Value;
            obj.reverse_point = State.Properties('reverse_point').Value;
            obj.Sg_r = State.Properties('residual_saturation').Value;
        end
              
    end
end
 
 

