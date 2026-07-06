% VTK Plotter
%--------------------------------------------------------------------------
% DARSim 2.5 Reservoir Simulator
% Author: Artur Castiel Reis de Souza
% TU Delft
%--------------------------------------------------------------------------
classdef VTK_Plotter < Plotter
    properties
        FileName
        isBinary
        ReferenceDepth
        PlotInterfaces
    end
    methods
        function obj = VTK_Plotter(Directory, Problem, NumOfPreviousReports)
            if ~exist(strcat(Directory,'/Output/VTK/'), 'dir')
                mkdir(Directory,'/Output/VTK');
            else
                delete(strcat(Directory,'/Output/VTK/*.vtk'));
            end
            obj.FileName = strcat(Directory, '/Output/VTK/', Problem);
            obj.VTKindex = NumOfPreviousReports;
        end
        function PlotWells(obj, Inj, Prod, Grid)
            %Injectors
            for i=1:length(Inj)
                name = strcat(obj.FileName,num2str(i),'Inj','.vtk');
                obj.WriteAWell(Inj{i}, name, Grid);
            end 
            %Producers
            for i=1:length(Prod)
                name = strcat(obj.FileName, num2str(i),'Prod','.vtk');
                obj.WriteAWell(Prod{i}, name, Grid);
            end
            
        end
        function WriteAWell(obj, Well, name, Grid)
            fileID = fopen(name, 'w');
            fprintf(fileID, '# vtk DataFile Version 2.0\n');
            fprintf(fileID, 'Matteo Simulator: Well s file\n');
            fprintf(fileID, '\n');
            fprintf(fileID, 'ASCII\n');
            fprintf(fileID, '\n');
            fprintf(fileID, 'DATASET POLYDATA\n');
            fprintf(fileID, ['POINTS ', num2str(length(Well.Cells) + 1), ' float\n']);
            %loop over perforations
            ztop = (Grid.Nz+1)*Grid.dz;
            for c=1:length(Well.Cells)
                [i, j, k] = ind2sub([Grid.Nx, Grid.Ny, Grid.Nz], Well.Cells(c));
                x = (i - 1)* Grid.dx + Grid.dx/2;
                y = (j - 1)* Grid.dy + Grid.dy/2;
                z = (k -1) * Grid.dz + Grid.dz/2;
                fprintf(fileID, '%10.5f %10.5f %10.5f\n', x, y, z);
            end
            fprintf(fileID, '%10.5f %10.5f %10.5f\n', x, y, ztop);
            fprintf(fileID, ['LINES ', num2str(length(Well.Cells)), ' ',num2str(3*(length(Well.Cells))), '\n']);
            for c = 1:length(Well.Cells)
                fprintf(fileID, '%d %d %d\n', [2, c-1, c]);
            end
            fclose(fileID);
        end
        function PlotSolution(obj, ProductionSystem, DiscretizationModel)
           obj.PlotReservoirSolution(ProductionSystem.Reservoir, DiscretizationModel.ReservoirGrid);
            for f = 1 : length(ProductionSystem.FracturesNetwork.Fractures)
                obj.PlotFractureSolution(ProductionSystem.FracturesNetwork.Fractures(f), DiscretizationModel.FracturesGrid.Grids(f), f);
            end
            if isempty(ProductionSystem.FracturesNetwork.Fractures)==0 && obj.VTKindex==0
                obj.PlotAllFracturesMesh(ProductionSystem.FracturesNetwork.Fractures, DiscretizationModel.FracturesGrid.Grids);
            end
        end
        function PlotReservoirSolution(obj, Reservoir, Grid)
            %Write a VTK file for Reservoir
            fileID = fopen(strcat(obj.FileName, '_Reservoir', num2str(obj.VTKindex,'%04d'),'.vtk'), 'w');
            fprintf(fileID, '# vtk DataFile Version 2.0\n');
            fprintf(fileID, 'DARSim2 Reservoir Simulator\n');
            if obj.isBinary
                fprintf(fileID, 'BINARY\n');
            else
            	fprintf(fileID, 'ASCII\n');
            end
            fprintf(fileID, '\n');
            %
            fprintf(fileID, 'DATASET RECTILINEAR_GRID\n');
            fprintf(fileID, 'DIMENSIONS    %d   %d   %d\n', Grid.Nx+1, Grid.Ny+1, Grid.Nz+1);
            fprintf(fileID, '\n');
            fprintf(fileID, ['X_COORDINATES ' num2str(Grid.Nx+1) ' float\n']);
            if obj.isBinary
                fwrite(fileID,[0:Grid.dx:Grid.dx * Grid.Nx],'float', 'b');
            else
            	fprintf(fileID, '%f ', 0:Grid.dx:Grid.dx * Grid.Nx);
            end
            fprintf(fileID, '\n');
            fprintf(fileID, ['Y_COORDINATES ' num2str(Grid.Ny+1) ' float\n']);
            if obj.isBinary
                fwrite(fileID,[0:Grid.dy:Grid.dy * Grid.Ny],'float', 'b');
            else
            	fprintf(fileID, '%f ', 0:Grid.dy:Grid.dy * Grid.Ny);
            end
            fprintf(fileID, '\n');
            fprintf(fileID, ['Z_COORDINATES ' num2str(Grid.Nz+1) ' float\n']);
            if obj.isBinary
                fwrite(fileID,[0:Grid.dz:Grid.dz * Grid.Nz],'float','b');
            else
            	fprintf(fileID, '%d ', 0:Grid.dz:Grid.dz * Grid.Nz);
            end
            
            fprintf(fileID, '\n');
            fprintf(fileID, '\n');
            %
            fprintf(fileID, 'CELL_DATA %d\n', Grid.N);
            fprintf(fileID, '\n');
            
            % Print all existing variables
            N_var = double(Reservoir.State.Properties.Count);
            Names = Reservoir.State.Properties.keys;
            for i=1:N_var
                if strcmp(Reservoir.State.Properties(Names{i}).Type, 'scalar')
                    obj.PrintScalar2VTK(fileID, Reservoir.State.Properties(Names{i}).Value, [' ',Names{i}]);
                else
                    obj.PrintVector2VTK(fileID, Reservoir.State.Properties(Names{i}).Value, [' ',Names{i}]);
                end
                fprintf(fileID, '\n');
            end
            
            % Add the "isFractured" flag for reservoir grid cells that are overlapped by a fracture (if any)
            if ~isempty(Grid.ListOfFracturedReservoirCells)
                FracturedFlag = zeros(Grid.N,1);
                FracturedFlag(Grid.ListOfFracturedReservoirCells) = 1;
                obj.PrintScalar2VTK(fileID, FracturedFlag, ' isFractured');
                fprintf(fileID, '\n');
            end
            
            % Add the "isPerforated" flag for reservoir grid cells that are overlapped by a fracture (if any)
            if ~isempty(Grid.ListOfPerforatedCells)
                isPerforated = zeros(Grid.N,1);
                isPerforated(Grid.ListOfPerforatedCells) = 1;
                obj.PrintScalar2VTK(fileID, isPerforated, ' isPerforated');
                fprintf(fileID, '\n');
            end
            
            % Add the cell indices
            obj.PrintScalar2VTK(fileID, [1:Grid.N]', ' Index');
            fprintf(fileID, '\n');
            
            % Add the cell volumes
            obj.PrintScalar2VTK(fileID, Grid.Volume, ' Volume');
            fprintf(fileID, '\n');

            % Add ADM ACTIVEFine (coarse grids)
            obj.PrintScalar2VTK(fileID, Grid.Active, ' ACTIVEFine');
            fprintf(fileID, '\n');
            
            % Add ADM ACTIVETime
            obj.PrintScalar2VTK(fileID, Grid.ActiveTime, ' ACTIVETime');
            fprintf(fileID, '\n');

            
            % Adding the fine cell type based on dual coarse grid construction
            obj.PrintScalar2VTK(fileID, Grid.DualCoarseType, ' DualCoarseType');
            fprintf(fileID, '\n');

            fclose(fileID);

            if obj.PlotInterfaces
               obj.PlotInternalFaces(Grid);
            end

        end

        function PlotInternalFaces(obj, Grid)

            fileID = fopen(strcat(obj.FileName, '_Interfaces', num2str(obj.VTKindex,'%04d'),'.vtk'), 'w');
            fprintf(fileID, '# vtk DataFile Version 2.0\n');
            fprintf(fileID, 'DARSim 2 Reservoir Simulator\n');
            if obj.isBinary
                fprintf(fileID, 'BINARY\n');
            else
                fprintf(fileID, 'ASCII\n');
            end
            fprintf(fileID, '\n');
            fprintf(fileID, 'DATASET UNSTRUCTURED_GRID\n');

            VTK_Nodes = Grid.GridCoords;
            if Grid.N ~= sum(Grid.Active)
                obj.ReferenceDepth = max(VTK_Nodes(:,3));
                VTK_Nodes(:,3) = obj.ReferenceDepth - VTK_Nodes(:,3);
            end

            fprintf(fileID, 'POINTS %d double\n', size(VTK_Nodes,1));
            if obj.isBinary
                fwrite(fileID, VTK_Nodes', 'double', 0, 'b');
            else
                fprintf(fileID, '%0.15g %0.15g %0.15g\n', VTK_Nodes');
            end
            fprintf(fileID, '\n');

            Nx = Grid.Nx;
            Ny = Grid.Ny;
            Nz = Grid.Nz;

            nx = Nx + 1;
            ny = Ny + 1;
            nz = Nz + 1;

            vid = @(i,j,k) i + j*nx + k*nx*ny + 1;

            nFx = (Nx-1)*Ny*Nz;
            nFy = Nx*(Ny-1)*Nz;
            nFz = Nx*Ny*(Nz-1);
            
            N_InternalFaces = nFx + nFy + nFz;
            faceVertsCell  = zeros(N_InternalFaces, 4);

            f = 1;
            
            for k = 0:Nz-1
                for j = 0:Ny-1
                    for i = 1:Nx-1
                        
                        v1 = vid(i,j,k);
                        v2 = vid(i,j+1,k);
                        v3 = vid(i,j+1,k+1);
                        v4 = vid(i,j,k+1);
                        
                        faceVertsCell(f,:) = [v1 v2 v3 v4];
                        f = f + 1;
                    end
                end
            end

            for k = 0:Nz-1
                for j = 1:Ny-1
                    for i = 0:Nx-1
                        
                        v1 = vid(i,j,k);
                        v2 = vid(i,j,k+1);
                        v3 = vid(i+1,j,k+1);
                        v4 = vid(i+1,j,k);
                        
                        faceVertsCell(f,:) = [v1 v2 v3 v4];
                        f = f + 1;
                    end
                end
            end

            for k = 1:Nz-1
                for j = 0:Ny-1
                    for i = 0:Nx-1
                        
                        v1 = vid(i,j,k);
                        v2 = vid(i+1,j,k);
                        v3 = vid(i+1,j+1,k);
                        v4 = vid(i,j+1,k);
                        
                        faceVertsCell(f,:) = [v1 v2 v3 v4];
                        f = f + 1;
                    end
                end
            end
        
            Nv = zeros(N_InternalFaces,1,'int32');
            cleanFace = cell(N_InternalFaces,1);
            for f = 1:N_InternalFaces
                v = faceVertsCell(f,:);
                v = v(:)';
                v = v(~isnan(v) & v>0);
                if ~isempty(v)
                  
                    v = v([true, diff(v)~=0]);
                end
                cleanFace{f} = v;
                Nv(f) = int32(numel(v));
            end

            totalInts = sum(double(Nv) + 1);
            fprintf(fileID, 'CELLS %d %d\n', N_InternalFaces, totalInts);

            if obj.isBinary
                for f = 1:N_InternalFaces
                    v = int32(cleanFace{f} - 1);
                    fwrite(fileID, int32([Nv(f), v]), 'int32', 0, 'b');
                end
            else
                for f = 1:N_InternalFaces
                    fprintf(fileID, '%d ', Nv(f));
                    fprintf(fileID, '%d ', cleanFace{f} - 1);
                    fprintf(fileID, '\n');
                end
            end
            fprintf(fileID, '\n');
            
            fprintf(fileID, 'CELL_TYPES %d\n', N_InternalFaces);
            if obj.isBinary
                fwrite(fileID, int32(7*ones(N_InternalFaces,1))', 'int32', 0, 'b'); % VTK_POLYGON
            else
                fprintf(fileID, '%d\n', 7*ones(N_InternalFaces,1));
            end
            fprintf(fileID, '\n\n');

            fprintf(fileID, 'CELL_DATA %d\n', N_InternalFaces);
            T_inner = [reshape(Grid.Tx(2:end-1,:,:),1,[])';...
                 reshape(Grid.Ty(:,2:end-1,:),1,[])';...
                 reshape(Grid.Tz(:,:,2:end-1),1,[])'];
            alpha_inner = [reshape(Grid.pEDFM_alpha_Tx(2:end-1,:,:),1,[])'; ...
                           reshape(Grid.pEDFM_alpha_Ty(:,2:end-1,:),1,[])'; ...
                           reshape(Grid.pEDFM_alpha_Tz(:,:,2:end-1),1,[])'];
            obj.PrintScalar2VTK(fileID, T_inner, ' Trans');
            obj.PrintScalar2VTK(fileID, alpha_inner, ' alpha');
            fprintf(fileID, '\n');

            fclose(fileID);
        end

        function PlotFractureSolution(obj, Fracture, Grid, f)
            %Write a VTK file for each fracture
            fileID = fopen(strcat(obj.FileName, '_Fracture', num2str(f,'%03d'), '_', num2str(obj.VTKindex,'%04d'),'.vtk'), 'w');
            fprintf(fileID, '# vtk DataFile Version 2.0\n');
            fprintf(fileID, 'DARSim 2 Reservoir Simulator\n');
            
            if obj.isBinary
                fprintf(fileID, 'BINARY\n');
            else
            	fprintf(fileID, 'ASCII\n');
            end
            fprintf(fileID, '\n');
            
            fprintf(fileID, 'DATASET STRUCTURED_GRID\n');
            fprintf(fileID, 'DIMENSIONS    %d   %d   %d\n', Grid.Nx+1, Grid.Ny+1, 1);
            fprintf(fileID, '\n');
            
            fprintf(fileID, 'POINTS    %d   double\n', size(Grid.GridCoords,1) );
            if obj.isBinary
                fwrite(fileID, Grid.GridCoords', 'double', 'b');
            else
                fprintf(fileID, '%f %f %f\n' , Grid.GridCoords'); 
            end
            fprintf(fileID, '\n');
            
            fprintf(fileID, '\n');
            fprintf(fileID, 'CELL_DATA %d\n', Grid.N);
            fprintf(fileID, '\n');
            
            %Add ADM coarse grids
            obj.PrintScalar2VTK(fileID, Grid.Active, ' ACTIVEFine');
            fprintf(fileID, '\n');
            N_var = double(Fracture.State.Properties.Count);
            Names = Fracture.State.Properties.keys;
            
            for i=1:N_var
                if strcmp(Fracture.State.Properties(Names{i}).Type, 'scalar')
                    obj.PrintScalar2VTK(fileID, Fracture.State.Properties(Names{i}).Value, [' ',Names{i}]);
                else
                    obj.PrintVector2VTK(fileID, Fracture.State.Properties(Names{i}).Value, [' ',Names{i}]);
                end
                fprintf(fileID, '\n');
            end
            
            % Adding the fine cell type based on dual coarse grid construction
            obj.PrintScalar2VTK(fileID, Grid.DualCoarseType, ' DualCoarseType');
            fprintf(fileID, '\n');
            
            fclose(fileID);
        end
        function PlotPermeability(obj, ProductionSystem, DiscretizationModel)
            obj.PlotReservoirPermeability(DiscretizationModel.ReservoirGrid, ProductionSystem.Reservoir.K);
            for f = 1 : length(ProductionSystem.FracturesNetwork.Fractures)
                obj.PlotFracturePermeability(DiscretizationModel.FracturesGrid.Grids(f), ProductionSystem.FracturesNetwork.Fractures(f).K, f);
            end
        end
        function PlotReservoirPermeability(obj, ReservoirGrid, K)
            % Permeability of the reservoir
            fileID = fopen(strcat(obj.FileName, '_Reservoir', num2str(obj.VTKindex,'%04d'),'.vtk'), 'a');
            obj.PrintScalar2VTK(fileID, reshape(K(:,1), ReservoirGrid.N, 1), ' PERMX');
            obj.PrintScalar2VTK(fileID, reshape(K(:,2), ReservoirGrid.N, 1), ' PERMY');
            obj.PrintScalar2VTK(fileID, reshape(K(:,3), ReservoirGrid.N, 1), ' PERMZ');
            fprintf(fileID, '\n');
            fclose(fileID);
        end
        function PlotFracturePermeability(obj, FractureGrid, K, f)
            % Permeability of each fracture
            fileID = fopen(strcat(obj.FileName, '_Fracture', num2str(f,'%03d'), '_', num2str(obj.VTKindex,'%04d'),'.vtk'), 'a');
            obj.PrintScalar2VTK(fileID, reshape(K(:,1), FractureGrid.N, 1), ' PERMX');
            obj.PrintScalar2VTK(fileID, reshape(K(:,2), FractureGrid.N, 1), ' PERMY');
            obj.PrintScalar2VTK(fileID, reshape(K(:,3), FractureGrid.N, 1), ' PERMZ');
            fprintf(fileID, '\n');
            fclose(fileID);
        end
        function PlotPorosity(obj, ProductionSystem, DiscretizationModel)
            obj.PlotReservoirPorosity(DiscretizationModel.ReservoirGrid, ProductionSystem.Reservoir.Por);
            for f = 1 : length(ProductionSystem.FracturesNetwork.Fractures)
                obj.PlotFracturePorosity(DiscretizationModel.FracturesGrid.Grids(f), ProductionSystem.FracturesNetwork.Fractures(f).Por, f);
            end
        end
        function PlotReservoirPorosity(obj, ReservoirGrid, phi)
            % Porosity of the reservoir
            fileID = fopen(strcat(obj.FileName, '_Reservoir', num2str(obj.VTKindex,'%04d'),'.vtk'), 'a');
            obj.PrintScalar2VTK(fileID, reshape(phi, ReservoirGrid.N, 1), ' Porosity');
            fprintf(fileID, '\n');
            fclose(fileID);
        end
        function PlotFracturePorosity(obj, FractureGrid, phi, f)
            % Porosity of each fracture
            fileID = fopen(strcat(obj.FileName, '_Fracture', num2str(f,'%03d'), '_', num2str(obj.VTKindex,'%04d'),'.vtk'), 'a');
            obj.PrintScalar2VTK(fileID, reshape(phi.*ones(FractureGrid.N,1), FractureGrid.N, 1), ' Porosity');
            fprintf(fileID, '\n');
            fclose(fileID);
        end
        function PlotBasisFunctions(obj,FineGrid, CoarseGrid, Prolp, Nf, Nc)
            for m = 1 : size(Nc,1)
                for L = 1 : size(Nc,2)
                    Nc(m,L) = Nc(m,L) * CoarseGrid(m,L).hasCoarseNodes;
                end
            end
            obj.PlotReservoirBF(FineGrid(1), CoarseGrid(1,:), Prolp);
            for i=2:length(FineGrid)
                obj.PlotFractureBF(FineGrid(i), CoarseGrid(i,:), Prolp, sum(Nf(1:i-1)), Nc, i);
            end
        end
        function PlotReservoirBF(obj, Grid, CoarseGrid, Prolp)
            %% 1. Level 1
            fileID = fopen(strcat(obj.FileName, '_Reservoir', '_BF_Level1.vtk'), 'w');
            fprintf(fileID, '# vtk DataFile Version 2.0\n');
            fprintf(fileID, 'DARSim 2 Reservoir Simulator\n');
            fprintf(fileID, 'BINARY\n');
            fprintf(fileID, '\n');
            fprintf(fileID, 'DATASET RECTILINEAR_GRID\n');
            fprintf(fileID, 'DIMENSIONS    %d   %d   %d\n', Grid.Nx +1, Grid.Ny+1, Grid.Nz+1);
            fprintf(fileID, '\n');
            fprintf(fileID, ['X_COORDINATES ' num2str(Grid.Nx+1) ' float\n']);
            %fprintf(fileID, '%f ', 0:Grid.dx:Grid.dx * Grid.Nx);
            fwrite(fileID, 0:Grid.dx:Grid.dx * Grid.Nx, 'float', 'b');
            fprintf(fileID, '\n');
            fprintf(fileID, ['Y_COORDINATES ' num2str(Grid.Ny+1) ' float\n']);
            %fprintf(fileID, '%f ', 0:Grid.dy:Grid.dy * Grid.Ny);
            fwrite(fileID, 0:Grid.dy:Grid.dy * Grid.Ny, 'float', 'b');
            fprintf(fileID, '\n');
            fprintf(fileID, ['Z_COORDINATES ' num2str(Grid.Nz+1) ' float\n']);
            %fprintf(fileID, '%d ', 0:Grid.dz:Grid.dz * Grid.Nz);
            fwrite(fileID, 0:Grid.dz:Grid.dz * Grid.Nz, 'float', 'b');
            fprintf(fileID, '\n');
            fprintf(fileID, '\n');
            fprintf(fileID, 'CELL_DATA   %d\n', Grid.N);
            fprintf(fileID, '\n');
            [~, n_columns] = size(Prolp{1});
            % Matrix basis functions in the matrix
            for j = 1:CoarseGrid(1).N
                obj.PrintScalar2VTK(fileID, full(Prolp{1}(1:Grid.N,j)),strcat(' BF',num2str(j,'%03d')));
                fprintf(fileID, '\n');
            end
            % Fracture basis functions in the matrix
            for j = CoarseGrid(1).N+1:n_columns
                obj.PrintScalar2VTK(fileID, full(Prolp{1}(1:Grid.N, j)), strcat(' Frac_BF',num2str(j,'%03d')));
                fprintf(fileID, '\n');
            end
            
            fclose(fileID);
            %% 2. Levels > 1
            for i=2:length(CoarseGrid)
                fileID = fopen(strcat(obj.FileName, '_Reservoir', '_BF_Level',num2str(i),'.vtk'), 'w');
                fprintf(fileID, '# vtk DataFile Version 2.0\n');
                fprintf(fileID, 'DARSim 2 Reservoir Simulator\n');
                fprintf(fileID, 'BINARY\n');
                fprintf(fileID, '\n');
                fprintf(fileID, 'DATASET RECTILINEAR_GRID\n');
                fprintf(fileID, 'DIMENSIONS    %d   %d   %d\n', CoarseGrid(i-1).Nx +1, CoarseGrid(i-1).Ny+1, CoarseGrid(i-1).Nz+1);
                fprintf(fileID, '\n');
                fprintf(fileID, ['X_COORDINATES ' num2str(CoarseGrid(i-1).Nx+1) ' float\n']);
                % fprintf(fileID, '%f ', 0 : Grid.dx * CoarseGrid(i-1).CoarseFactor(1) : Grid.dx * Grid.Nx);
                fwrite(fileID, 0 : Grid.dx * CoarseGrid(i-1).CoarseFactor(1) : Grid.dx * Grid.Nx, 'float', 'b');
                fprintf(fileID, '\n');
                fprintf(fileID, ['Y_COORDINATES ' num2str(CoarseGrid(i-1).Ny+1) ' float\n']);
                %fprintf(fileID, '%f ', 0 : Grid.dy * CoarseGrid(i-1).CoarseFactor(2) : Grid.dy * Grid.Ny);
                fwrite(fileID, 0 : Grid.dy * CoarseGrid(i-1).CoarseFactor(2) : Grid.dy * Grid.Ny, 'float', 'b');
                fprintf(fileID, '\n');
                fprintf(fileID, ['Z_COORDINATES ' num2str(CoarseGrid(i-1).Nz+1) ' float\n']);
                %fprintf(fileID, '%f ', 0 : Grid.dz * CoarseGrid(i-1).CoarseFactor(3) : Grid.dz * Grid.Nz);
                fwrite(fileID, 0 : Grid.dz * CoarseGrid(i-1).CoarseFactor(3) : Grid.dz * Grid.Nz, 'float', 'b');
                fprintf(fileID, '\n');
                fprintf(fileID, '\n');
                fprintf(fileID, 'CELL_DATA   %d\n', CoarseGrid(i-1).N);
                [~, n_columns] = size(Prolp{i});
                % Matrix basis functions in the matrix
                for j = 1:CoarseGrid(i).N
                    obj.PrintScalar2VTK(fileID, full(Prolp{i}(1:CoarseGrid(i-1).N, j)), strcat(' BF',num2str(j,'%03d')));
                    fprintf(fileID, '\n');
                end
                % Fracture basis functions in the matrix
                for j = CoarseGrid(i).N+1:n_columns
                    obj.PrintScalar2VTK(fileID, full(Prolp{i}(1:CoarseGrid(i-1).N, j)), strcat(' Frac_BF',num2str(j,'%03d')));
                    fprintf(fileID, '\n');
                end
                fclose(fileID);
            end
        end
        function PlotFractureBF(obj, Grid, CoarseGrid, Prolp, Nf, Nc, f)
            %% 1. Level 1
            fileID = fopen(strcat(obj.FileName, '_Fracture', num2str(f-1),'_BF_Level1.vtk'), 'w');
            fprintf(fileID, '# vtk DataFile Version 2.0\n');
            fprintf(fileID, 'DARSim 2 Reservoir Simulator\n');
            fprintf(fileID, 'BINARY\n');
            fprintf(fileID, '\n');
            fprintf(fileID, 'DATASET STRUCTURED_GRID\n');
            fprintf(fileID, 'DIMENSIONS    %d   %d   %d\n', Grid.Nx+1, Grid.Ny+1, 1);
            fprintf(fileID, '\n');
            fprintf(fileID, 'POINTS    %d   double\n', size(Grid.GridCoords,1) );
            %fprintf(fileID, '%f %f %f\n' , Fracture.GridCoords'); 
            fwrite(fileID, Grid.GridCoords', 'double', 'b');
            fprintf(fileID, '\n');
            fprintf(fileID, '\n');
            fprintf(fileID, 'CELL_DATA %d\n', Grid.N);
            fprintf(fileID, '\n');
            [~, n_columns] = size(Prolp{1});
            % Matrix basis functions in the fracture f
            for j = 1:Nc(1,1)
                obj.PrintScalar2VTK(fileID, full(Prolp{1}(Nf+1:Nf+Grid.N,j)),strcat(' BF',num2str(j,'%03d')));
                fprintf(fileID, '\n');
            end
            % Fractures basis functions in the fracture
            if CoarseGrid(1).hasCoarseNodes
                for j = Nc(1,1)+1:n_columns
                    obj.PrintScalar2VTK(fileID, full(Prolp{1}(Nf+1:Nf+Grid.N, j)), strcat(' Frac_BF',num2str(j,'%03d')));
                    fprintf(fileID, '\n');
                end
            end
            fclose(fileID);
            %% 2. Levels > 1
            for i=2:length(CoarseGrid)
                fileID = fopen(strcat(obj.FileName, '_Fracture', num2str(f-1),'_BF_Level',num2str(i),'.vtk'), 'w');
                fprintf(fileID, '# vtk DataFile Version 2.0\n');
                fprintf(fileID, 'DARSim 2 Reservoir Simulator\n');
                fprintf(fileID, 'BINARY\n');
                fprintf(fileID, '\n');
                fprintf(fileID, 'DATASET STRUCTURED_GRID\n');
                fprintf(fileID, 'DIMENSIONS    %d   %d   %d\n', CoarseGrid(i-1).Nx+1, CoarseGrid(i-1).Ny+1, 1);
                fprintf(fileID, '\n');
                fprintf(fileID, 'POINTS    %d   double\n', size(CoarseGrid(i-1).GridCoords, 1) );
                fwrite(fileID, CoarseGrid(i-1).GridCoords', 'double', 'b');
                fprintf(fileID, '\n');
                fprintf(fileID, 'CELL_DATA %d\n', CoarseGrid(i-1).N);
                fprintf(fileID, '\n');
                [~, n_columns] = size(Prolp{i});
                % Matrix basis functions in the fracture f
                for j = 1:Nc(1,i)
                    obj.PrintScalar2VTK(fileID, full(Prolp{i}( sum(Nc(1:f-1, i-1))+1:sum(Nc(1:f-1, i-1)) + CoarseGrid(i-1).N*CoarseGrid(i-1).hasCoarseNodes, j) ), strcat(' BF',num2str(j,'%03d')));
                    fprintf(fileID, '\n');
                end
                % Fracture basis functions in the fracture f
                if CoarseGrid(i).hasCoarseNodes
                    for j = Nc(1,i)+1:n_columns
                        obj.PrintScalar2VTK(fileID, full(Prolp{i}(sum(Nc(1:f-1, i-1))+1:sum(Nc(1:f-1, i-1)) + CoarseGrid(i-1).N, j)), strcat(' Frac_BF',num2str(j,'%03d')));
                        fprintf(fileID, '\n');
                    end
                end
                fclose(fileID);
            end
        end
        function PlotDynamicBasisFunctions(obj, Grid, Prolp)
            fileID = fopen(strcat(obj.FileName, '_Reservoir', 'Dynamic_BF_', num2str(obj.VTKindex,'%04d'),'.vtk'), 'w');
            fprintf(fileID, '# vtk DataFile Version 2.0\n');
            fprintf(fileID, 'DARSim 2 Reservoir Simulator\n');
            fprintf(fileID, 'BINARY\n');
            fprintf(fileID, '\n');
            fprintf(fileID, 'DATASET RECTILINEAR_GRID\n');
            fprintf(fileID, 'DIMENSIONS    %d   %d   %d\n', Grid.Nx +1, Grid.Ny+1, Grid.Nz+1);
            fprintf(fileID, '\n');
            fprintf(fileID, ['X_COORDINATES ' num2str(Grid.Nx+1) ' float\n']);
            %fprintf(fileID, '%f ', 0:Grid.dx:Grid.dx * Grid.Nx);
            fwrite(fileID, 0:Grid.dx:Grid.dx * Grid.Nx, 'float', 'b');
            fprintf(fileID, '\n');
            fprintf(fileID, ['Y_COORDINATES ' num2str(Grid.Ny+1) ' float\n']);
            %fprintf(fileID, '%f ', 0:Grid.dy:Grid.dy * Grid.Ny);
            fwrite(fileID, 0:Grid.dy:Grid.dy * Grid.Ny, 'float', 'b');
            fprintf(fileID, '\n');
            fprintf(fileID, ['Z_COORDINATES ' num2str(Grid.Nz+1) ' float\n']);
            %fprintf(fileID, '%d ', 0:Grid.dz:Grid.dz * Grid.Nz);
            fwrite(fileID, 0:Grid.dz:Grid.dz * Grid.Nz, 'float', 'b');
            fprintf(fileID, '\n');
            fprintf(fileID, '\n');
            fprintf(fileID, 'CELL_DATA   %d\n', Grid.N);
            fprintf(fileID, '\n');
            if obj.VTKindex > 2
                [~, col] = size(Prolp);
                for j = 1:col
                    obj.PrintScalar2VTK(fileID, full(Prolp(:,j)),strcat(' DynamicBF',num2str(j)));
                    fprintf(fileID, '\n');
                end
                fclose(fileID);
            end
        end
        function PlotSaturationInterpolator(obj, Grid, ProlS, Pdelta, Pdeltac)
            %Write a VTK file for Reservoir
            fileID = fopen(strcat(obj.FileName, '_Reservoir' ,'SatInterp', num2str(obj.VTKindex,'%04d'),'.vtk'), 'w');
            fprintf(fileID, '# vtk DataFile Version 2.0\n');
            fprintf(fileID, 'DARSim 2 Reservoir Simulator\n');
            %fprintf(fileID, 'ASCII\n');
            fprintf(fileID, 'BINARY\n');
            fprintf(fileID, '\n');
            fprintf(fileID, 'DATASET RECTILINEAR_GRID\n');
            fprintf(fileID, 'DIMENSIONS    %d   %d   %d\n', Grid.Nx+1, Grid.Ny+1, Grid.Nz+1);
            fprintf(fileID, '\n');
            fprintf(fileID, ['X_COORDINATES ' num2str(Grid.Nx+1) ' float\n']);
            %fprintf(fileID, '%f ', 0:Grid.dx:Grid.dx * Grid.Nx);
            fwrite(fileID,[0:Grid.dx:Grid.dx * Grid.Nx],'float', 'b');
            fprintf(fileID, '\n');
            fprintf(fileID, ['Y_COORDINATES ' num2str(Grid.Ny+1) ' float\n']);
            fwrite(fileID,[0:Grid.dy:Grid.dy * Grid.Ny],'float', 'b');
            %fprintf(fileID, '%f ', 0:Grid.dy:Grid.dy * Grid.Ny);
            fprintf(fileID, '\n');
            fprintf(fileID, ['Z_COORDINATES ' num2str(Grid.Nz+1) ' float\n']);
            %fprintf(fileID, '%d ', 0:Grid.dz:Grid.dz * Grid.Nz);
            fwrite(fileID,[0:Grid.dz:Grid.dz * Grid.Nz],'float','b');
            fprintf(fileID, '\n');
            fprintf(fileID, '\n');
            fprintf(fileID, 'CELL_DATA   %d\n', Grid.N);
            fprintf(fileID, '\n');
            [~, col] = size(ProlS);
            if obj.VTKindex >= 2
                %obj.PrintScalar2VTK(fileID, Pdelta, ' Delta_fine');
                obj.PrintScalar2VTK(fileID, full(sum(ProlS, 2)), ' BF_S');
            %for i=1:length(ProlS)
             %   [~, col] = size(ProlS);
              %  obj.PrintScalar2VTK(fileID, Pdeltac{i}, strcat(' Level',num2str(i),'Deltac'));
                %for j = 1:col
                 %   obj.PrintScalar2VTK(fileID, full(ProlS{i}(:,j)), strcat(' Level',num2str(i),'Node',num2str(j)));
                  %  fprintf(fileID, '\n');
                %end 
            %end
            fclose(fileID);
            end
        end
        function PlotADMGrid(obj, ProductionSystem, DiscretizationModel)
            %% Plot ADM Grid
            % 1. Reservoir
            obj.PlotReservoirADMGrid(DiscretizationModel.CoarseGrid(1,:));
            if ~isempty(ProductionSystem.Reservoir.K_coarse)
                obj.PlotCoarsePermeability(DiscretizationModel.CoarseGrid(1,:), ProductionSystem.Reservoir.K_coarse);
            end
            % 2. Fractures
            for f = 1 : length(ProductionSystem.FracturesNetwork.Fractures)
                obj.PlotFractureADMGrid(DiscretizationModel.CoarseGrid(1+f,:), f);
            end
        end
        function PlotReservoirADMGrid(obj, CoarseGrid)
            for i=1:length(CoarseGrid)
                fileID = fopen(strcat(obj.FileName, '_Reservoir', num2str(i),'Level', num2str(obj.VTKindex,'%04d'),'.vtk'), 'w');
                fprintf(fileID, '# vtk DataFile Version 2.0\n');
                fprintf(fileID, 'DARSim 2 Reservoir Simulator\n');
                fprintf(fileID, 'BINARY\n');
                fprintf(fileID, '\n');
                fprintf(fileID, 'DATASET STRUCTURED_GRID\n');
                fprintf(fileID, 'DIMENSIONS    %d   %d   %d\n', CoarseGrid(i).Nx+1, CoarseGrid(i).Ny+1, CoarseGrid(i).Nz+1);
                fprintf(fileID, '\n');				  
                fprintf(fileID, 'POINTS    %d   double\n', size(CoarseGrid(i).GridCoords, 1) );
                fwrite(fileID, CoarseGrid(i).GridCoords', 'double', 'b');														 
                fprintf(fileID, '\n');
                fprintf(fileID, '\n');
                fprintf(fileID, 'CELL_DATA   %d\n', CoarseGrid(i).N);
                obj.PrintScalar2VTK(fileID, CoarseGrid(i).Active.*CoarseGrid(i).hasCoarseNodes, ' ActiveCoarse');
                fprintf(fileID, '\n');
                
                % Adding the fine cell type based on dual coarse grid construction
                obj.PrintScalar2VTK(fileID, CoarseGrid(i).DualCoarseType, ' DualCoarseType');
                fprintf(fileID, '\n');
                
                % Delta_S
                obj.PrintScalar2VTK(fileID, CoarseGrid(i).DeltaS, ' Delta_S');
                fprintf(fileID, '\n');
                fclose(fileID);
            end
        end
        function PlotFractureADMGrid(obj, CoarseGrid, f)
            %Write a VTK file for each
            for i=1:length(CoarseGrid)
                fileID = fopen(strcat(obj.FileName, '_Fracture', num2str(f,'%03d'),...
                    '_',num2str(i),'Level_',num2str(obj.VTKindex,'%04d'),'.vtk'), 'w');
                fprintf(fileID, '# vtk DataFile Version 2.0\n');
                fprintf(fileID, 'DARSim 2 Reservoir Simulator\n');
                fprintf(fileID, 'BINARY\n');
                fprintf(fileID, '\n');
                fprintf(fileID, 'DATASET STRUCTURED_GRID\n');
                fprintf(fileID, 'DIMENSIONS    %d   %d   %d\n', CoarseGrid(i).Nx+1, CoarseGrid(i).Ny+1, 1);
                fprintf(fileID, '\n');
                fprintf(fileID, 'POINTS    %d   double\n', size(CoarseGrid(i).GridCoords, 1) );
                fwrite(fileID, CoarseGrid(i).GridCoords', 'double', 'b');
                fprintf(fileID, '\n');
                fprintf(fileID, '\n');
                fprintf(fileID, 'CELL_DATA %d\n', CoarseGrid(i).N);
                fprintf(fileID, '\n');
                obj.PrintScalar2VTK(fileID, CoarseGrid(i).Active.*CoarseGrid(i).hasCoarseNodes, ' ActiveCoarse');
                fprintf(fileID, '\n');
                
                obj.PrintScalar2VTK(fileID, CoarseGrid(i).DualCoarseType, ' DualCoarseType');
                fprintf(fileID, '\n');
                
                fclose(fileID);
            end
        end
        function PlotCoarsePermeability(obj, CoarseGrid, K_coarse)
            for i=1:length(CoarseGrid)
                fileID = fopen(strcat(obj.FileName, '_Reservoir', num2str(i),'Level', num2str(obj.VTKindex,'%04d'),'.vtk'), 'a');
                obj.PrintScalar2VTK(fileID, reshape(K_coarse{i+1}(:,1), CoarseGrid(i).N, 1), ' PERMX');
                obj.PrintScalar2VTK(fileID, reshape(K_coarse{i+1}(:,2), CoarseGrid(i).N, 1), ' PERMY');
                obj.PrintScalar2VTK(fileID, reshape(K_coarse{i+1}(:,3), CoarseGrid(i).N, 1), ' PERMZ');
                fprintf(fileID, '\n');
                fclose(fileID);
            end
        end
        function PrintScalar2VTK(obj, fileID, scalar, name)
            %Print a scalar in VTK format
            fprintf(fileID, '\n');
            fprintf(fileID, strcat('SCALARS  ', name,' double 1\n'));
            fprintf(fileID, 'LOOKUP_TABLE default\n');
            if obj.isBinary
                fwrite(fileID, scalar', 'double', 'b');
            else
            	fprintf(fileID,'%1.5e ', scalar);
            end
        end
        function PrintVector2VTK(obj, fileID, vector, name)
            %Print a vector in VTK format
            fprintf(fileID, '\n');
            fprintf(fileID, strcat('VECTORS  ', name,' double \n'));
            if obj.isBinary
                fwrite(fileID, vector','double', 'b');
            else
            	fprintf(fileID,'%1.5e ', vector);
            end
        end
        function PlotAllFracturesMesh(obj, Fractures, Grids)
            %Write a VTK file for each fracture


            %%
            %warning("plot all fractures needs to be rewritten")
            fileID = fopen(strcat(obj.FileName, '_AllFracturesMesh','.vtk'), 'w');
            
            VTK_Nodes = []; VTK_Cells = [];
            
            for i=1:size(Grids,1)
                
                Coord = Grids(i).GridCoords;
                
                posmin = find( Coord( :, 1 ) == min( Coord( :, 1 ) ) );
                
                posmax = find( Coord( :, 1 ) == max( Coord( :, 1 ) ) );
                
                Cells = [ 1 2 3 4 ] + 4*(i-1);
                
                Nodes = [ min(Coord(:,1)) Coord(posmin(1),2) min(Coord(:,3)); ...
                          max(Coord(:,1)) Coord(posmax(1),2) min(Coord(:,3)); ...
                          max(Coord(:,1)) Coord(posmax(1),2) max(Coord(:,3)); ...
                          min(Coord(:,1)) Coord(posmin(1),2) max(Coord(:,3))];
            
                VTK_Nodes = [ VTK_Nodes; Nodes ];
                
                VTK_Cells = [ VTK_Cells; Cells ];
                
            end
                        
            if ~isempty(obj.ReferenceDepth)
                VTK_Nodes(:,3) = obj.ReferenceDepth - VTK_Nodes(:,3);
            end
            
            fprintf(fileID, '# vtk DataFile Version 2.0\n');
            fprintf(fileID, 'DARSim 2 Reservoir Simulator\n');
            
            if obj.isBinary
                fprintf(fileID, 'BINARY\n');
            else
            	fprintf(fileID, 'ASCII\n');
            end
            fprintf(fileID, '\n');
            
            fprintf(fileID, 'DATASET UNSTRUCTURED_GRID\n');

            fprintf(fileID, ['POINTS ' num2str(size(VTK_Nodes,1)) ' double\n']);
            fprintf(fileID, '%f %f %f\n',VTK_Nodes');
            fprintf(fileID, '\n');
            
            N_Ret = size(VTK_Cells,1);
            maxNumVertices = 4;
            fprintf(fileID, ['CELLS ' num2str(N_Ret) ' ' num2str(N_Ret*(1+maxNumVertices)) '\n']);
            
            indMatrix = horzcat(maxNumVertices*ones(N_Ret,1) , VTK_Cells-1 );
            
            FormatSpec =[];
            for n = 1 : maxNumVertices + 1
                FormatSpec = strcat(FormatSpec,"%d ");
            end
            FormatSpec = char(FormatSpec);
            FormatSpec(end)=[];
            FormatSpec = strcat(FormatSpec,'\n');
            fprintf(fileID, FormatSpec,indMatrix');
            fprintf(fileID, '\n');
            
            fprintf(fileID, ['CELL_TYPES ' num2str(N_Ret) '\n']);
            fprintf(fileID, '%d ', 9*ones(1,N_Ret));
            fprintf(fileID, '\n\n');

            fclose(fileID);
        end
    end
end