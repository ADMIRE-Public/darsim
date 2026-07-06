% VTK Plotter
%--------------------------------------------------------------------------
% DARSim 2 Reservoir Simulator
% Author: Túlio de M. Cavalcante
% UFPE
%--------------------------------------------------------------------------
classdef TetrahedralGrid_VTK_Plotter < VTK_Plotter
    methods
        function WriteAWell(obj, Well, name, Grid)
        end
        function PlotReservoirSolution(obj, Reservoir, Grid)
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
            VTK_Nodes = Grid.TetrahedralGridData.Vertex.Coordinates;
            if Grid.N ~= size(Grid.Active,1)
                obj.ReferenceDepth = max(VTK_Nodes(:,3));
                % If the domain is a Cartesian box, there is no need to invert the Z coordinates.
                % But if the domain is not a Cartesian boxTherefore, the Z coordinates are inverted only for the sake of the visualization.
                VTK_Nodes(:,3) = obj.ReferenceDepth - VTK_Nodes(:,3);
            end
            % POINTS ------------------------------------------------------
            fprintf(fileID, ['POINTS ' num2str(size(VTK_Nodes,1)) ' double\n']);
            if obj.isBinary
                fwrite(fileID, VTK_Nodes','float', 'b');
            else
                fprintf(fileID, '%f %f %f\n',VTK_Nodes');
            end
            fprintf(fileID, '\n');
            %--------------------------------------------------------------
            % CELLS -------------------------------------------------------
            fprintf(fileID, ['CELLS ' num2str(Grid.N) ' ' num2str(Grid.N*(1+4)) '\n']);
            indMatrix = horzcat(4*ones(Grid.N,1) , Grid.TetrahedralGridData.Cell.Vertices-1 );
            if obj.isBinary
                fwrite(fileID, indMatrix','float', 'b');
            else
                fprintf(fileID, '%d %d %d %d %d\n', indMatrix');
            end            
            fprintf(fileID, '\n');
            %--------------------------------------------------------------
            % CELL TYPES --------------------------------------------------
            fprintf(fileID, ['CELL_TYPES ' num2str(Grid.N) '\n']);
            if obj.isBinary
                fwrite(fileID, 10*ones(1,Grid.N),'float', 'b');
            else
                fprintf(fileID, '%d ', 10*ones(1,Grid.N));
            end            
            fprintf(fileID, '\n\n');
            %--------------------------------------------------------------
            % PRINT ALL VARIABLES -----------------------------------------
            fprintf(fileID, 'CELL_DATA %d\n', Grid.N);
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
            %--------------------------------------------------------------
            % IS FRACTURED? -----------------------------------------------
            if ~isempty(Grid.ListOfFracturedReservoirCells)
                isFractured = zeros(Grid.N,1);
                isFractured(Grid.ListOfFracturedReservoirCells) = 1;
                obj.PrintScalar2VTK(fileID, isFractured, ' isFractured');
                fprintf(fileID, '\n');
            end
            %--------------------------------------------------------------
            % IS PERFORATED? ----------------------------------------------
            if ~isempty(Grid.ListOfPerforatedCells)
                isPerforated = zeros(Grid.N,1);
                isPerforated(Grid.ListOfPerforatedCells) = 1;
                obj.PrintScalar2VTK(fileID, isPerforated, ' isPerforated');
                fprintf(fileID, '\n');
            end
            %--------------------------------------------------------------
            % CELL INDEX --------------------------------------------------
            obj.PrintScalar2VTK(fileID, [1:Grid.N]', ' Index');
            fprintf(fileID, '\n');
            %--------------------------------------------------------------
            % CELL VOLUME -------------------------------------------------
            obj.PrintScalar2VTK(fileID, Grid.Volume, ' Volume');
            fprintf(fileID, '\n');
            %--------------------------------------------------------------
            % ADM ACTIVE TIME ---------------------------------------------
            obj.PrintScalar2VTK(fileID, Grid.ActiveTime, ' ACTIVETime');
            fprintf(fileID, '\n');
            %--------------------------------------------------------------
            % ADM ACTIVE TIME FINE ----------------------------------------
            obj.PrintScalar2VTK(fileID, Grid.Active, ' ACTIVEFine');
            fprintf(fileID, '\n');
            %--------------------------------------------------------------
                        
            fclose(fileID);
            
            if obj.PlotInterfaces
                obj.PlotInternalFaces(Reservoir, Grid);
            end
        end
        function PlotInternalFaces(obj, Reservoir, Grid)
            %Write a VTK file for Reservoir
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
            
            VTK_Nodes = Grid.TetrahedralGridData.Nodes;
            
            if Grid.TetrahedralGridData.N_TotalCells ~= Grid.TetrahedralGridData.N_ActiveCells
                obj.ReferenceDepth = max(VTK_Nodes(:,3));
                % If the domain is a Cartesian box, there is no need to invert the Z coordinates.
                % But if the domain is not a Cartesian boxTherefore, the Z coordinates are inverted only for the sake of the visualization.
                VTK_Nodes(:,3) = obj.ReferenceDepth - VTK_Nodes(:,3);
            end

            fprintf(fileID, ['POINTS ' num2str(Grid.TetrahedralGridData.N_Nodes) ' double\n']);
            fprintf(fileID, '%f %f %f\n',VTK_Nodes');
            fprintf(fileID, '\n');
            
            N_InternalFaces = size(Grid.TetrahedralGridData.Face.Inner.Vertices,1);
            maxNumVertices = 3;
            FaceVertices = cell2mat( cellfun(@(x) [x ones(1,maxNumVertices-numel(x))*x(end)] , Grid.TetrahedralGridData.Internal_Faces.Vertices , 'uni' , 0 ) );
            fprintf(fileID, ['CELLS ' num2str(N_InternalFaces) ' ' num2str(N_InternalFaces*(1+maxNumVertices)) '\n']);
            
            indMatrix = horzcat(maxNumVertices*ones(N_InternalFaces,1) , FaceVertices-1 );
            
            FormatSpec =[];
            for n = 1 : maxNumVertices
                FormatSpec = strcat(FormatSpec,"%d ");
            end
            FormatSpec = char(FormatSpec);
            FormatSpec(end)=[];
            FormatSpec = strcat(FormatSpec,'\n');
            fprintf(fileID, FormatSpec,indMatrix');
            fprintf(fileID, '\n');
            
            fprintf(fileID, ['CELL_TYPES ' num2str(N_InternalFaces) '\n']);
            fprintf(fileID, '%d ', 5*ones(1,N_InternalFaces));
            fprintf(fileID, '\n\n');
            
            % Print all existing variables
            fprintf(fileID, 'CELL_DATA %d\n', N_InternalFaces);
            N_var = 1;
            Names{1} = 'Trans';
            for i=1:N_var
                obj.PrintScalar2VTK(fileID, Grid.Trans, [' ',Names{i}]);
                fprintf(fileID, '\n');
            end

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
            
            VTK_Nodes = Grid.GridCoords;
            if ~isempty(obj.ReferenceDepth)
                % If the domain is a Cartesian box, there is no need to invert the Z coordinates.
                % But if the domain is not a Cartesian boxTherefore, the Z coordinates are inverted only for the sake of the visualization.
                % The objects in this function do not know the corner-point grid information.
                % Therefore, the "if" condition here checks whether the "obj.ReferenceDepth" is not empty.
                % Because a similar "if condition" in the "PlotReservoirSolution" function above adds a value
                % to the "obj.ReferenceDepth" property. This way, we know that the condition already applies.
                VTK_Nodes(:,3) = obj.ReferenceDepth - VTK_Nodes(:,3);
            end
            fprintf(fileID, 'POINTS    %d   double\n', size(Grid.GridCoords,1) );
            if obj.isBinary
                fwrite(fileID, VTK_Nodes', 'double', 'b');
            else
                fprintf(fileID, '%f %f %f\n' , VTK_Nodes'); 
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
            fclose(fileID);
        end
        function PlotAllFracturesMesh(obj, Fractures, Grids)
            
            for i=1:size(Grids,1)
               
                fracK(i,1) = max( max( Fractures(i).K ) ); 
                
            end
            
            unifracK = unique(fracK); t = 0;
            
            for ik=unifracK'
               
                lista = find( fracK == ik ); t = t + 1;
                
                fileID = fopen(strcat(obj.FileName, '_AllFracturesMesh_',sprintf('%u',t),'.vtk'), 'w');
            
                VTK_Nodes = []; VTK_Cells = [];

                for i=1:size(lista,1)

                    Coord = Grids(lista(i)).GridCoords;

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
end

%--------------------------------------------------------------------------
% CELL TYPE
% vertex = 1;
% vertices set = 2;
% P1 segment = 3;
% polyline = 4;
% P1 triangle = 5;
% triangle strip = 6;
% polygon = 7;
% Q1 rectangle = 8;
% Q1 quadrangle = 9;
% P1 tetrahedron = 10; 
% Q1 cube = 11;
% Q1 hexahedron = 12;
% P1 prism = 13;
% P1 pyramid = 14;
% P2 segment = 21;
% P2 triangle = 22;
% Q2 serendip quadrangle = 23;
% P2 tetrahedron = 24;
% Q2 serendip hexahedron = 25;
%--------------------------------------------------------------------------
% vtk format is designed for first and second order element, unfortunately 
% for quadrangle and hexahedron it is not the standard quadratic element! 
% To avoid intricate conversions from Q2 to Q2 serendip element, only 
% first order elements mesh may be exported with vtk format. When Mesh 
% is of order greater than one, a one order mesh is created by spliting 
% elements into first order elements; the underlying mesh being attached 
% to its mesh parent using Mesh::firstOrderMesh_p pointer. Note that if 
% the mesh is order one, firstOrderMesh_p points to it.
%--------------------------------------------------------------------------
% vtk format does not manage multiple domains. By default only the main 
% (whole) domain is exported If boolean argument multipleDomains is true, 
% all the domains are exported to different "files"
%--------------------------------------------------------------------------