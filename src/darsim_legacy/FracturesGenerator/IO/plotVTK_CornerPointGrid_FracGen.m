% VTK Plotter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 11 July 2016
%Last modified: 17 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef plotVTK_CornerPointGrid_FracGen < plotVTK_FracGen
    properties
    end
    methods
        function PlotReservoirSolution(obj, Grid)
            %Write a VTK file for Reservoir
            fileID = fopen(strcat(obj.FileName, '_Reservoir', num2str(obj.VTKindex,'%04d'),'.vtk'), 'w');
            fprintf(fileID, '# vtk DataFile Version 2.0\n');
            fprintf(fileID, 'DARSim 2 Reservoir Simulator\n');
            if obj.isBinary
                fprintf(fileID, 'BINARY\n');
            else
            	fprintf(fileID, 'ASCII\n');
            end
            fprintf(fileID, '\n');
            %
            fprintf(fileID, 'DATASET UNSTRUCTURED_GRID\n');
            VTK_Nodes = Grid.CornerPointGridData.Nodes;
            obj.ReferenceDepth = max(VTK_Nodes(:,3));
            VTK_Nodes(:,3) = obj.ReferenceDepth - VTK_Nodes(:,3);
            fprintf(fileID, ['POINTS ' num2str(Grid.CornerPointGridData.N_Nodes) ' double\n']);
            if obj.isBinary
                fwrite(fileID, VTK_Nodes','float', 'b');
            else
                fprintf(fileID, '%f %f %f\n',VTK_Nodes');
            end
            fprintf(fileID, '\n');
            %
            fprintf(fileID, ['CELLS ' num2str(Grid.N) ' ' num2str(Grid.N*(1+8)) '\n']);
            indMatrix = horzcat(8*ones(Grid.N,1) , Grid.CornerPointGridData.Cells.Vertices-1 );
            if obj.isBinary
                fwrite(fileID, indMatrix','float', 'b');
            else
                fprintf(fileID, '%d %d %d %d %d %d %d %d %d\n', indMatrix');
            end
            
            fprintf(fileID, '\n');
            %
            fprintf(fileID, ['CELL_TYPES ' num2str(Grid.N) '\n']);
            if obj.isBinary
                fwrite(fileID, 11*ones(1,Grid.N),'float', 'b');
            else
                fprintf(fileID, '%d ', 11*ones(1,Grid.N));
            end
            
            fprintf(fileID, '\n\n');
            
            %%% Adding properties
            fprintf(fileID, 'CELL_DATA %d\n', Grid.N);
            
            % Add the "isFractured" flag for reservoir grid cells that are overlapped by a fracture (if any)
            if ~isempty(Grid.OvelapList)
                isFractured = zeros(Grid.N,1);
                isFractured(Grid.OvelapList) = 1;
                obj.PrintScalar2VTK(fileID, isFractured, ' isFractured');
                fprintf(fileID, '\n');
            end
            
            % Add the cell indices
            obj.PrintScalar2VTK(fileID, [1:Grid.N]', ' Index');
            fprintf(fileID, '\n');

            % Add the cell volumes
            obj.PrintScalar2VTK(fileID, Grid.CellVolume, ' Volume');
            fprintf(fileID, '\n');

            fclose(fileID);
        end 
        function PlotFracturesSolution(obj, Fracture, f)
            %Write a VTK file for each
            fileID = fopen(strcat(obj.FileName, '_Fracture', num2str(f,'%03d'), '_', num2str(obj.VTKindex),'.vtk'), 'w');
            fprintf(fileID, '# vtk DataFile Version 2.0\n');
            fprintf(fileID, 'DARSim 2 Reservoir Simulator\n');
            
            if obj.isBinary
                fprintf(fileID, 'BINARY\n');
            else
            	fprintf(fileID, 'ASCII\n');
            end
            fprintf(fileID, '\n');
            
            fprintf(fileID, 'DATASET STRUCTURED_GRID\n');
            fprintf(fileID, 'DIMENSIONS    %d   %d   %d\n', Fracture.N_Length_AB+1, Fracture.N_Width_AD+1, 1);
            fprintf(fileID, '\n');
            
            fprintf(fileID, 'POINTS    %d   double\n', size(Fracture.GridCoords,1)*size(Fracture.GridCoords,2) );
            GridCoords_X = Fracture.GridCoords(:,:,1);
            GridCoords_Y = Fracture.GridCoords(:,:,2);
            GridCoords_Z = obj.ReferenceDepth - Fracture.GridCoords(:,:,3);
            if obj.isBinary
                fwrite(fileID, [GridCoords_X(:),GridCoords_Y(:),GridCoords_Z(:)]', 'double', 'b');
            else
                fprintf(fileID, '%f %f %f\n' , [GridCoords_X(:),GridCoords_Y(:),GridCoords_Z(:)]' );
            end
            fprintf(fileID, '\n');

            fprintf(fileID, '\n');
            fprintf(fileID, 'CELL_DATA %d\n', Fracture.N_Length_AB*Fracture.N_Width_AD);
            fprintf(fileID, '\n');
            fclose(fileID);
        end
    end
end