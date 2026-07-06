% VTK Plotter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 11 July 2016
%Last modified: 17 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef plotVTK_CarstesianGrid_FracGen < plotVTK_FracGen
    properties
    end
    methods
        function PlotReservoirSolution(obj, Grid)
            %Write a VTK file for Reservoir
            fileID = fopen(strcat(obj.FileName, '_Reservoir_', num2str(obj.VTKindex),'.vtk'), 'w');
            fprintf(fileID, '# vtk DataFile Version 2.0\n');
            fprintf(fileID, 'DARSim 2 Reservoir Simulator\n');
            
            if obj.isBinary
                fprintf(fileID, 'BINARY\n');
            else
            	fprintf(fileID, 'ASCII\n');
            end
            fprintf(fileID, '\n');
            
            fprintf(fileID, 'DATASET RECTILINEAR_GRID\n');
            fprintf(fileID, 'DIMENSIONS    %d   %d   %d\n', Grid.Nx+1, Grid.Ny+1, Grid.Nz+1);
            fprintf(fileID, '\n');
            
            fprintf(fileID, ['X_COORDINATES ' num2str(Grid.Nx+1) ' float\n']);
            if obj.isBinary
                fwrite(fileID,(0:Grid.dx:Grid.dx * Grid.Nx),'float', 'b');
            else
            	fprintf(fileID, '%f ', 0:Grid.dx:Grid.dx * Grid.Nx);
            end
            fprintf(fileID, '\n');
            
            fprintf(fileID, ['Y_COORDINATES ' num2str(Grid.Ny+1) ' float\n']);
            if obj.isBinary
                fwrite(fileID,(0:Grid.dy:Grid.dy * Grid.Ny),'float', 'b');
            else
            	fprintf(fileID, '%f ', 0:Grid.dy:Grid.dy * Grid.Ny);
            end
            fprintf(fileID, '\n');
            
            fprintf(fileID, ['Z_COORDINATES ' num2str(Grid.Nz+1) ' float\n']);
            if obj.isBinary
                fwrite(fileID,(0:Grid.dz:Grid.dz * Grid.Nz),'float', 'b');
            else
            	fprintf(fileID, '%f ', 0:Grid.dz:Grid.dz * Grid.Nz);
            end
            fprintf(fileID, '\n');
            
            fprintf(fileID, '\n');
            
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
            GridCoords_Z = Fracture.GridCoords(:,:,3);
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