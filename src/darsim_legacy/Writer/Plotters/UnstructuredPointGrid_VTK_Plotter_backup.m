% VTK Plotter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Artur Castiel
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef UnstructuredPointGrid_VTK_Plotter < VTK_Plotter
    properties
        PlotInterfaces
    end
    methods
        function WriteAWell(obj, Well, name, Grid)
        end
        function PlotReservoirSolution(obj, Reservoir, Grid)
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
            %VTK_Nodes = Grid.CornerPointGridData.Nodes;
            VTK_Nodes = Grid.mesh.vertex.coord;
            % if Grid.CornerPointGridData.N_TotalCells ~= Grid.CornerPointGridData.N_ActiveCells
            %     obj.ReferenceDepth = max(VTK_Nodes(:,3));
            %     % If the domain is a Cartesian box, there is no need to invert the Z coordinates.
            %     % But if the domain is not a Cartesian boxTherefore, the Z coordinates are inverted only for the sake of the visualization.
            %     VTK_Nodes(:,3) = obj.ReferenceDepth - VTK_Nodes(:,3);
            % end
            %fprintf(fileID, ['POINTS ' num2str(Grid.mesh.vertex.num) ' float\n']);
            %if obj.isBinary
            % %   fwrite(fileID, VTK_Nodes','float', 'b');
            %else
             %   fprintf(fileID, '%f %f %f\n',VTK_Nodes');
            %end
            %fprintf(fileID, '\n');
            %
            fprintf(fileID, 'POINTS    %d   double\n', size(Grid.GridCoords,1) );
            if obj.isBinary
                fwrite(fileID, Grid.GridCoords', 'double', 'b');
            else
                fprintf(fileID, '%f %f %f\n' , Grid.GridCoords'); 
            end
            fprintf(fileID, '\n');
            
            lfaces = [Grid.mesh.face.inner.vertices; Grid.mesh.face.bound.vertices];
            num_nodes = sum(lfaces ~=0,2);
            node_group = unique(num_nodes)';
            [~,ref_order] = sort(num_nodes);
            num_nodes = num_nodes(ref_order);
            lfaces = lfaces(ref_order,:);
            %num_lfaces = size(lfaces, 1);
            %fprintf(fileID,'CELLS %d %d\n', num_lfaces, sum(num_nodes)+num_lfaces);
            %fprintf(fileID, ['CELLS ' num2str(num_lfaces) ' ' num2str(sum(num_nodes) + num_lfaces) '\n']);
            %fprintf(fileID, ['CELLS ' num2str(Grid.N) ' ' num2str(Grid.N * 4 + Grid.N) '\n']);
            %indMatrix = horzcat(num_nodes, lfaces-1);
            %indMatrix = horzcat(8*ones(Grid.N,1) , Grid.CornerPointGridData.Cells.Vertices-1 );
            verts = Grid.mesh.element.vertices;      
            if any(isnan(verts(:)))
                nv = sum(~isnan(verts), 2);
                verts(isnan(verts)) = 0;
            else
                nv = sum(verts ~= 0, 2);
            end

            totalInts = sum(nv + 1);                 
            fprintf(fileID, 'CELLS %d %d\n', Grid.N, totalInts);

            if obj.isBinary
                for c = 1:Grid.N
                    v = verts(c,1:nv(c)) - 1;         % VTK 0-based
                    fwrite(fileID, int32([nv(c), v]), 'int32', 0, 'b');
                end
                fprintf(fileID, '\n');
            else
                for c = 1:Grid.N
                    fprintf(fileID, '%d ', nv(c));
                    fprintf(fileID, '%d ', verts(c,1:nv(c)) - 1);
                    fprintf(fileID, '\n');
                end
                fprintf(fileID, '\n');
            end

            cellTypes = zeros(Grid.N,1,'int32');
            if isfield(Grid.mesh, 'element') && isfield(Grid.mesh.element, 'type')
                t = Grid.mesh.element.type(:);
                cellTypes(t==1) = 10;   % VTK_TETRA
                cellTypes(t==2) = 12;   % VTK_HEXAHEDRON
                cellTypes(t==3) = 13;   % VTK_WEDGE
                cellTypes(t==4) = 14;
            else
                for c = 1:Grid.N
                    switch nv(c)
                        case 4
                            cellTypes(c) = 10;  % TETRA
                        case 8
                            cellTypes(c) = 12;  % HEXAHEDRON
                        case 6
                            cellTypes(c) = 13;  % WEDGE/PRISM
                        case 5
                            cellTypes(c) = 14;  % PYRAMID
                        otherwise
                            error('Unsupported cell: cell %d has %d vertices. Need polyhedron support.', c, nv(c));
                    end
                end
            end

            fprintf(fileID, 'CELL_TYPES %d\n', Grid.N);
            if obj.isBinary
                fwrite(fileID, cellTypes.', 'int32', 0, 'b');
            else
                fprintf(fileID, '%d ', cellTypes);
            end
            fprintf(fileID, '\n\n');

            % indMatrix = [repmat(4, Grid.N, 1),Grid.mesh.element.vertices-1];  
            % if obj.isBinary
            %     fwrite(fileID, int32(indMatrix.'), 'int32', 0, 'b');
            %     fprintf(fileID, '\n');
            % else
            %     fprintf(fileID, '%d %d %d %d %d\n', indMatrix.');
            %     fprintf(fileID, '\n');
            % end

            % for index = node_group
            %     lref = indMatrix(:,1) == index;
            %     lindMatrix = indMatrix(lref,1:index+1);
            %     string_pattern = [repmat('%d ', 1, index +1) '\n'];
            %     if obj.isBinary
            %         fwrite(fileID, lindMatrix','int32', 'b');
            %     else
            %         fprintf(fileID, string_pattern, lindMatrix');
            %     end
            % end
            %fprintf(fileID, '\n');
            %
            %fprintf(fileID, ['CELL_TYPES ' num2str(num_lfaces) '\n']);
            % fprintf(fileID, 'CELL_TYPES %d\n', Grid.N);
            % if obj.isBinary
            %     fwrite(fileID, 10*ones(1,Grid.N),'int32', 'b');
            % else
            %     fprintf(fileID, '%d ', 10*ones(1,Grid.N));
            % end
            
            %fprintf(fileID, '\n\n');
            
            % Print all existing variables
            fprintf(fileID, 'CELL_DATA %d\n', Grid.N);
            N_var = double(Reservoir.State.Properties.Count);
            Names = Reservoir.State.Properties.keys;
            all_neigh = [Grid.mesh.face.inner.upstream Grid.mesh.face.inner.downstream; Grid.mesh.face.bound.upstream Grid.mesh.face.bound.upstream*0];
            all_neigh = all_neigh(ref_order,:);
            for i=1:N_var
                type_var = Reservoir.State.Properties(Names{i}).Type;                
                lvalue = Reservoir.State.Properties(Names{i}).Value;
                %[face_value] = face_average(all_neigh,lvalue);                          
                if strcmp(type_var, 'scalar')
                    obj.PrintScalar2VTK(fileID, lvalue, [' ',Names{i}]);
                else
                    obj.PrintVector2VTK(fileID, lvalue, [' ',Names{i}]);
                end
                fprintf(fileID, '\n');
            end

            % Add the "isFractured" flag for reservoir grid cells that are overlapped by a fracture (if any)
            if ~isempty(Grid.ListOfFracturedReservoirCells)
                %isFractured = any(ismember(all_neigh,Grid.ListOfFracturedReservoirCells),2);
                isFractured = zeros(Grid.N,1);
                isFractured(Grid.ListOfFracturedReservoirCells) = 1;
                obj.PrintScalar2VTK(fileID, isFractured, ' isFractured');
                fprintf(fileID, '\n');
            end
            
            % Add the "isPerforated" flag for reservoir grid cells that are overlapped by a fracture (if any)
            if ~isempty(Grid.ListOfPerforatedCells)
                isPerforated = zeros(Grid.N,1);
                isPerforated(Grid.ListOfPerforatedCells) = 1;
                %isPerforated = any(ismember(all_neigh,Grid.ListOfPerforatedCells),2);
                obj.PrintScalar2VTK(fileID, isPerforated, ' isPerforated');
                fprintf(fileID, '\n');
            end
            
            % Add the cell indices
            
           % obj.PrintScalar2VTK(fileID, [1:num_lfaces]', ' Face_Index');
           % fprintf(fileID, '\n');

            % Add the cell volumes
           % [face_value] = face_average(all_neigh,Grid.Volume);   
            obj.PrintScalar2VTK(fileID, Grid.Volume, ' Volume');
            fprintf(fileID, '\n');
            
            % Add ADM ACTIVETime
            
            %act_time = face_average(all_neigh,Grid.ActiveTime);
            obj.PrintScalar2VTK(fileID, Grid.ActiveTime, ' ACTIVETime');
            fprintf(fileID, '\n');
            
            % Add ADM ACTIVEFine (coarse grids)
            %active = face_average(all_neigh,Grid.Active);
            obj.PrintScalar2VTK(fileID, Grid.Active, ' ACTIVEFine');
            fprintf(fileID, '\n');
            
            % Add manual flag for R5 model
%             if Grid.N == 1276945
%                 Flagged = zeros(Grid.N,1);
%                 Cell_Indices = [959588;1029937;1033734;1193170;1266538];
%                 Flagged(Cell_Indices)=1;
%                 for n = 1 : length(Cell_Indices)
%                     Flagged( Grid.CornerPointGridData.Cells.Neighbors{n} ) = 2;
%                 end
%                 obj.PrintScalar2VTK(fileID, Flagged, ' isFlagged');
%                 fprintf(fileID, '\n');
%             end
            
            fclose(fileID);
            
            if obj.PlotInterfaces
                obj.PlotInternalFaces(Reservoir, Grid);
            end
        end
        function PlotReservoirSolution2(obj, Reservoir, Grid)
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
            %VTK_Nodes = Grid.CornerPointGridData.Nodes;
            VTK_Nodes = Grid.mesh.vertex.coord;
            % if Grid.CornerPointGridData.N_TotalCells ~= Grid.CornerPointGridData.N_ActiveCells
            %     obj.ReferenceDepth = max(VTK_Nodes(:,3));
            %     % If the domain is a Cartesian box, there is no need to invert the Z coordinates.
            %     % But if the domain is not a Cartesian boxTherefore, the Z coordinates are inverted only for the sake of the visualization.
            %     VTK_Nodes(:,3) = obj.ReferenceDepth - VTK_Nodes(:,3);
            % end
            fprintf(fileID, ['POINTS ' num2str(Grid.mesh.vertex.num) ' double\n']);
            if obj.isBinary
                fwrite(fileID, VTK_Nodes','float', 'b');
            else
                fprintf(fileID, '%f %f %f\n',VTK_Nodes');
            end
            fprintf(fileID, '\n');
            %
            fprintf(fileID, ['CELLS ' num2str(Grid.N) ' ' num2str(Grid.N*(1+8)) '\n']);
            indMatrix = horzcat(8*ones(Grid.N,1) , Grid.mesh.element.vertices(:,1:8)-1);
            %indMatrix = horzcat(8*ones(Grid.N,1) , Grid.CornerPointGridData.Cells.Vertices-1 );
            
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
            
            % Print all existing variables
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

            % Add the "isFractured" flag for reservoir grid cells that are overlapped by a fracture (if any)
            if ~isempty(Grid.ListOfFracturedReservoirCells)
                isFractured = zeros(Grid.N,1);
                isFractured(Grid.ListOfFracturedReservoirCells) = 1;
                obj.PrintScalar2VTK(fileID, isFractured, ' isFractured');
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
            
            % Add ADM ACTIVETime
            obj.PrintScalar2VTK(fileID, Grid.ActiveTime, ' ACTIVETime');
            fprintf(fileID, '\n');
            
            % Add ADM ACTIVEFine (coarse grids)
            obj.PrintScalar2VTK(fileID, Grid.Active, ' ACTIVEFine');
            fprintf(fileID, '\n');
            
            % Add manual flag for R5 model
%             if Grid.N == 1276945
%                 Flagged = zeros(Grid.N,1);
%                 Cell_Indices = [959588;1029937;1033734;1193170;1266538];
%                 Flagged(Cell_Indices)=1;
%                 for n = 1 : length(Cell_Indices)
%                     Flagged( Grid.CornerPointGridData.Cells.Neighbors{n} ) = 2;
%                 end
%                 obj.PrintScalar2VTK(fileID, Flagged, ' isFlagged');
%                 fprintf(fileID, '\n');
%             end
            
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
            
            VTK_Nodes = Grid.CornerPointGridData.Nodes;
            
            if Grid.CornerPointGridData.N_TotalCells ~= Grid.CornerPointGridData.N_ActiveCells
                obj.ReferenceDepth = max(VTK_Nodes(:,3));
                % If the domain is a Cartesian box, there is no need to invert the Z coordinates.
                % But if the domain is not a Cartesian boxTherefore, the Z coordinates are inverted only for the sake of the visualization.
                VTK_Nodes(:,3) = obj.ReferenceDepth - VTK_Nodes(:,3);
            end

            fprintf(fileID, ['POINTS ' num2str(Grid.CornerPointGridData.N_Nodes) ' double\n']);
            fprintf(fileID, '%f %f %f\n',VTK_Nodes');
            fprintf(fileID, '\n');
            
            N_InternalFaces = Grid.CornerPointGridData.N_InternalFaces;
            % maxNumVertices = max(Grid.CornerPointGridData.Internal_Faces.N_Vertices);
            % FaceVertices = cell2mat( cellfun(@(x) [x ones(1,maxNumVertices-numel(x))*x(end)] , Grid.CornerPointGridData.Internal_Faces.Vertices , 'uni' , 0 ) );
            % fprintf(fileID, ['CELLS ' num2str(N_InternalFaces) ' ' num2str(N_InternalFaces*(1+maxNumVertices)) '\n']);
            % 
            % indMatrix = horzcat(maxNumVertices*ones(N_InternalFaces,1) , FaceVertices-1 );
            Nv = Grid.CornerPointGridData.Internal_Faces.N_Vertices(:);
            totalInts = sum(Nv + 1);
            fprintf(fileID, 'CELLS %d %d\n', N_InternalFaces, totalInts);

            for f = 1:N_InternalFaces
                v = Grid.CornerPointGridData.Internal_Faces.Vertices{f} - 1; % 0-based
                fprintf(fileID, '%d ', numel(v));
                fprintf(fileID, '%d ', v);
                fprintf(fileID, '\n');
            end

            fprintf(fileID, 'CELL_TYPES %d\n', N_InternalFaces);
            fprintf(fileID, '%d\n', 7*ones(N_InternalFaces,1)); % 7 = VTK_POLYGON

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
            fprintf(fileID, '%d ', 9*ones(1,N_InternalFaces));
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
            
            len_max_frac = size(Grid.fracmesh.elements.vertices,2);
            fprintf(fileID, '\n');
            
            fprintf(fileID, 'DATASET UNSTRUCTURED_GRID\n');
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
            %%
            fprintf(fileID, ['CELLS ' num2str(Grid.N) ' ' num2str(Grid.N*(1+len_max_frac)) '\n']);
            indMatrix = horzcat((len_max_frac)*ones(Grid.N,1) , Grid.fracmesh.elements.vertices-1);
            %indMatrix = horzcat(8*ones(Grid.N,1) , Grid.CornerPointGridData.Cells.Vertices-1 );
            if obj.isBinary
                fwrite(fileID, int32(indMatrix'), 'int32', 0, 'b');
            else
                if len_max_frac == 3
                    fprintf(fileID, '%d %d %d %d\n', indMatrix');
                else
                    fprintf(fileID, '%d %d %d %d %d\n', indMatrix');
                end
            end
            fprintf(fileID, '\n');
            %
            if len_max_frac == 3
                celltype  = 5;
            else
                celltype = 9;
            end

            fprintf(fileID, ['CELL_TYPES ' num2str(Grid.N) '\n']);
            if obj.isBinary
                fwrite(fileID, 5*ones(1,Grid.N),'int32', 'b');
            else
                fprintf(fileID, '%d ', celltype*ones(1,Grid.N));
            end

            fprintf(fileID, '\n\n');
           
            %% CELL DATA
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
    end
end