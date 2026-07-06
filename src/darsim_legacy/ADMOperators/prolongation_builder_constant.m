%  Prolongation builder constant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef prolongation_builder_constant < prolongation_builder
    properties
        Pdelta
        Pdeltac
    end
    methods
        function obj = prolongation_builder_constant(n)
            obj@prolongation_builder(n)
            obj.R = cell(n, 1);
            obj.P = cell(n, 1);
            obj.C = cell(n, 1);
        end
        function BuildStaticOperators(obj, ProductionSystem, FluidModel, FineGrid, CrossConnections, maxLevel, CoarseGrid)
            % Build Restriction and Prolongation operators for static grids
            % Initialise
            obj.R = cell(maxLevel(1), 1);
            obj.P = cell(maxLevel(1), 1);
            
            %% Level 1
            % Build static restriction operator (FV)
            disp('Building Pressure Restriction 1');
            start1 = tic;
            obj.R{1} = obj.MsRestriction(FineGrid, CoarseGrid(:,1));
            % Build Prolongation operator1
            disp('Building Pressure Prolongation 1');
            obj.P{1} = obj.R{1}';
            obj.C{1} = [];
            
            %% Level 2+
            for x = 2:maxLevel(1)
                % Build static restriction operator (FV)
                disp(['Building Pressure Restriction ', num2str(x)]);
                obj.R{x} = obj.MsRestriction(CoarseGrid(:, x-1), CoarseGrid(:, x));
                % Build Prolongation operator
                disp(['Building Pressure Prolongation ', num2str(x)]);
                obj.P{x} = obj.R{x}';
                obj.C{x} = [];
            end
            
            StaticOperators = toc(start1);
            disp(['Pressure static operators built in: ', num2str(StaticOperators), ' s']);
        end
        function UpdateProlongationOperator(obj, FineGrid, CoarseGrid, ProductionSystem)
            % For constant bf no update is necessary
        end
        function ADMProl = ADMProlongation(obj, ADMGrid, GlobalGrids, ADMRest)
            % Since it s constant interpolation it is just transpose(R)
            ADMProl = ADMRest';
        end
        function AverageMassOnCoarseBlocks(obj, Formulation, DiscretizationModel, ProductionSystem, FineGrid, FluidModel, ADMRest)
            % Average mass on a coarse block
            Formulation.AverageMassOnCoarseBlocks(DiscretizationModel, ProductionSystem, FineGrid, FluidModel, ADMRest);
        end
        function StaticMultilevelPressureGuess(obj, ProductionSystem, FluidModel, DiscretizationModel)
            % virtual call
        end
    end
end