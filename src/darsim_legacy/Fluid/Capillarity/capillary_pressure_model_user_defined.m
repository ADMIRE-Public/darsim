% Table capillary pressure model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Yuhang Wang
%TU Delft
%Created: 3 Feb 2021
%Last modified: 14 July 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef capillary_pressure_model_user_defined < capillary_pressure_model
    properties
        %%%%%%%%%%%
        Pc
        dPc
        %%%%%%%%%%%
    end
    methods
        function Initialise(obj, ProductionSystem)
            N = length(ProductionSystem.Reservoir.Por);
            % allocate memory
            obj.Pc = zeros(N,1);
            obj.dPc = zeros(N,1); 
            
        end

        function ComputePcHysteresis_old(obj, RelPermModel, DiscretizationModel, State)
            S = State.Properties('S_1').Value;
            status_P_1 = State.Properties('P_1');
            status_P_2 = State.Properties('P_2');
            status_Pc = State.Properties('Pc');
            % Compute Pc as a function of S_1
            
%             % no hysteresis
%             obj.Pc =  RelPermModel.Pc_drainage_Sg_table(S);
%             obj.dPc = -RelPermModel.dPc_drainage_Sg_table(S);

            % with hysteresis
            for i = 1:length(S)
                if (RelPermModel.process_label(i)==0)
                    obj.Pc(i) =  RelPermModel.Pc_drainage_Sg_table(S(i));
                    obj.dPc(i) = -RelPermModel.dPc_drainage_Sg_table(S(i)); 
                else
                   scanning_curve_Pc_table = griddedInterpolant(RelPermModel.scanning_curve_Sg(i,:),RelPermModel.scanning_curve_Pc(i,:));
                   scanning_curve_dPc_table = griddedInterpolant(RelPermModel.scanning_curve_Sg(i,:),RelPermModel.scanning_curve_dPc(i,:)); 
                   obj.Pc(i) = scanning_curve_Pc_table(S(i));
                   obj.dPc(i) = -scanning_curve_dPc_table(S(i));
                end 
            end
            
            obj.Pc(obj.Pc<0)=0;
            % convergence consideration
            obj.Pc(DiscretizationModel.ReservoirGrid.ListOfFracturedReservoirCells) = obj.Pc(DiscretizationModel.ReservoirGrid.ListOfFracturedReservoirCells)*0.01;
            obj.dPc(DiscretizationModel.ReservoirGrid.ListOfFracturedReservoirCells) = obj.dPc(DiscretizationModel.ReservoirGrid.ListOfFracturedReservoirCells)*0.01;
            
            status_Pc.Value = obj.Pc;
            % update pressure;
            status_P_1.Value = status_P_2.Value + obj.Pc;
        end
        
        %% Vectorizing for loops in ComputePcHysteresis
        function ComputePcHysteresis(obj, RelPermModel, Reservoir, DiscretizationModel)
            S = Reservoir.State.Properties('S_1').Value;
            status_P_1 = Reservoir.State.Properties('P_1');
            status_P_2 = Reservoir.State.Properties('P_2');
            status_Pc = Reservoir.State.Properties('Pc');
            % Compute Pc as a function of S_1
            
%             % no hysteresis
%             obj.Pc =  RelPermModel.Pc_drainage_Sg_table(S);
%             obj.dPc = -RelPermModel.dPc_drainage_Sg_table(S);
            for i = 2:7
 
                % % with hysteresis
                % % on primary drainage
                % sand_dra_label = logical((obj.Sand_index == i) .* (RelPermModel.process_label == 0));
                % obj.Pc(sand_dra_label == 1) =  RelPermModel.Pc_drainage_Sg_table{i}(S(sand_dra_label == 1));
                % obj.dPc(sand_dra_label == 1) = -RelPermModel.dPc_drainage_Sg_table{i}(S(sand_dra_label == 1));
                % 
                % % on scanning curves
                % sand_imb_label = logical((obj.Sand_index == i) .* (RelPermModel.process_label == 1));
                % obj.Pc(sand_imb_label == 1) = RelPermModel.scanning_curve_Pc_table{i}(S(sand_imb_label == 1), RelPermModel.reverse_point(sand_imb_label == 1));
                % obj.dPc(sand_imb_label == 1) = -RelPermModel.scanning_curve_dPc_table{i}(S(sand_imb_label == 1), RelPermModel.reverse_point(sand_imb_label == 1));
          
                sand_label = RelPermModel.SandType_index == i;
                obj.Pc(sand_label) =  RelPermModel.Pc_drainage_Sg_table{i}(S(sand_label));
                obj.dPc(sand_label) = -RelPermModel.dPc_drainage_Sg_table{i}(S(sand_label));

%                 obj.Pc(obj.Sand_index == 2) = 1471.5;
%                 obj.Pc(obj.Sand_index == 3) = 294.3; 
%                 obj.Pc(obj.Sand_index == 4) = 98.1; 
%                 
%                 obj.dPc(obj.Sand_index == i) = 0;
                
                obj.Pc(obj.Pc<0)=0;
                % convergence consideration
                obj.Pc(DiscretizationModel.ReservoirGrid.ListOfFracturedReservoirCells) = obj.Pc(DiscretizationModel.ReservoirGrid.ListOfFracturedReservoirCells)*0.01;
                obj.dPc(DiscretizationModel.ReservoirGrid.ListOfFracturedReservoirCells) = obj.dPc(DiscretizationModel.ReservoirGrid.ListOfFracturedReservoirCells)*0.01;
            end
            
            status_Pc.Value = obj.Pc;
            % update pressure;
            status_P_1.Value = status_P_2.Value + obj.Pc;
        end

        function ComputePcnoHysteresis(obj, RelPermModel, Reservoir, DiscretizationModel)
            S = Reservoir.State.Properties('S_1').Value;
            status_P_1 = Reservoir.State.Properties('P_1');
            status_P_2 = Reservoir.State.Properties('P_2');
            status_Pc = Reservoir.State.Properties('Pc');
            % Compute Pc as a function of S_1
            
            % no hysteresis
            Pc_value =  RelPermModel.Pc_drainage_Sg_table(S);
            % dPc = -RelPermModel.dPc_drainage_Sg_table(S);

            Pc_value(Pc_value<0)=0;
            status_Pc.Value = Pc_value;
            % update pressure;
            status_P_1.Value = status_P_2.Value + Pc_value;
        end        
        function dPc = dPcdSHysteresis(obj)
            dPc = obj.dPc;
        end
        function dPc = dPcdSnoHysteresis(obj, RelPermModel, State)
            S = State.Properties('S_1').Value; 
            % no hysteresis
            dPc = -RelPermModel.dPc_drainage_Sg_table(S);
        end 
        function d2Pc = dPcdSdS(obj, S)
        end
        
    end
    
end