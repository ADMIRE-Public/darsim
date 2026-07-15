% VTK Plotter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Mengjie
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef UnstructuredPointGrid_VTK_Plotter < VTK_Plotter

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
            fprintf(fileID, 'DATASET UNSTRUCTURED_GRID\n');

            verts_raw = Grid.mesh.element.vertices;
            maxVid = obj.local_max_vertex_id(verts_raw);

            [VTK_Nodes, whichPts] = obj.local_pick_points(Grid, maxVid);

            fprintf(fileID, 'POINTS    %d   double\n', size(VTK_Nodes,1));
            if obj.isBinary
                fwrite(fileID, VTK_Nodes', 'double', 0, 'b');
            else
                fprintf(fileID, '%0.15g %0.15g %0.15g\n', VTK_Nodes');
            end
            fprintf(fileID, '\n');

     
            [cellVerts, nv] = obj.local_clean_connectivity(verts_raw);

            cellTypes = zeros(Grid.N,1,'int32');
            if isprop(Grid.mesh, 'element') && isprop(Grid.mesh.element, 'type')
                t = Grid.mesh.element.type(:);
                cellTypes(t==1) = int32(10);   % VTK_TETRA
                cellTypes(t==2) = int32(12);   % VTK_HEXAHEDRON
                cellTypes(t==3) = int32(13);   % VTK_WEDGE
                cellTypes(t==4) = int32(14);   % VTK_PYRAMID
                miss = (cellTypes==0);
                if any(miss)
                    cellTypes(miss) = obj.local_infer_cell_types(nv(miss)); % fallback
                end
            else
                cellTypes = obj.local_infer_cell_types(nv);
            end

            cellVerts = obj.wedge_order_correction(VTK_Nodes, cellTypes, cellVerts);

            totalInts = sum(double(nv) + 1);
            fprintf(fileID, 'CELLS %d %d\n', Grid.N, totalInts);

            if obj.isBinary
                for c = 1:Grid.N
                    v = int32(cellVerts{c} - 1);          % could be row or column
                    fwrite(fileID, int32([nv(c); v(:)]), 'int32', 0, 'b');
                end
            else
                for c = 1:Grid.N
                    fprintf(fileID, '%d ', nv(c));
                    fprintf(fileID, '%d ', cellVerts{c} - 1);
                    fprintf(fileID, '\n');
                end
            end
            fprintf(fileID, '\n');

            fprintf(fileID, 'CELL_TYPES %d\n', Grid.N);
            if obj.isBinary
                fwrite(fileID, cellTypes.', 'int32', 0, 'b');
            else
                fprintf(fileID, '%d\n', cellTypes);
            end
            fprintf(fileID, '\n\n');

    
            fprintf(fileID, 'CELL_DATA %d\n', Grid.N);
            N_var = double(Reservoir.State.Properties.Count);
            Names = Reservoir.State.Properties.keys;

            lfaces = [Grid.mesh.face.inner.vertices; Grid.mesh.face.bound.vertices];
            num_nodes = sum(lfaces ~=0,2);
            node_group = unique(num_nodes)';
            [~,ref_order] = sort(num_nodes);
            num_nodes = num_nodes(ref_order);
            lfaces = lfaces(ref_order,:); 

            all_neigh = [Grid.mesh.face.inner.upstream Grid.mesh.face.inner.downstream; ...
                         Grid.mesh.face.bound.upstream Grid.mesh.face.bound.upstream*0];
            all_neigh = all_neigh(ref_order,:);

            for i=1:N_var
                type_var = Reservoir.State.Properties(Names{i}).Type;
                lvalue = Reservoir.State.Properties(Names{i}).Value;
                if strcmp(type_var, 'scalar')
                    obj.PrintScalar2VTK(fileID, lvalue, [' ',Names{i}]);
                else
                    obj.PrintVector2VTK(fileID, lvalue, [' ',Names{i}]);
                end
                fprintf(fileID, '\n');
            end

            if ~isempty(Grid.ListOfFracturedReservoirCells)
                isFractured = zeros(Grid.N,1);
                isFractured(Grid.ListOfFracturedReservoirCells) = 1;
                obj.PrintScalar2VTK(fileID, isFractured, ' isFractured');
                fprintf(fileID, '\n');
            end

            if ~isempty(Grid.ListOfPerforatedCells)
                isPerforated = zeros(Grid.N,1);
                isPerforated(Grid.ListOfPerforatedCells) = 1;
                obj.PrintScalar2VTK(fileID, isPerforated, ' isPerforated');
                fprintf(fileID, '\n');
            end

            obj.PrintScalar2VTK(fileID, Grid.Volume, ' Volume');        fprintf(fileID, '\n');
            obj.PrintScalar2VTK(fileID, Grid.ActiveTime, ' ACTIVETime');fprintf(fileID, '\n');
            obj.PrintScalar2VTK(fileID, Grid.Active, ' ACTIVEFine');    fprintf(fileID, '\n');

            fclose(fileID);

            if obj.PlotInterfaces
               obj.PlotInternalFaces(Grid);
            end
        end

        function PlotReservoirSolution2(obj, Reservoir, Grid)

            fileID = fopen(strcat(obj.FileName, '_Reservoir', num2str(obj.VTKindex,'%04d'),'.vtk'), 'w');

            fprintf(fileID, '# vtk DataFile Version 2.0\n');
            fprintf(fileID, 'DARSim 2 Reservoir Simulator\n');
            if obj.isBinary
                fprintf(fileID, 'BINARY\n');
            else
                fprintf(fileID, 'ASCII\n');
            end
            fprintf(fileID, '\n');
            fprintf(fileID, 'DATASET UNSTRUCTURED_GRID\n');

            verts_raw = Grid.mesh.element.vertices;
            maxVid = obj.local_max_vertex_id(verts_raw);
            [VTK_Nodes, whichPts] = obj.local_pick_points(Grid, maxVid);

            fprintf(fileID, 'POINTS %d double\n', size(VTK_Nodes,1));
            if obj.isBinary
                fwrite(fileID, VTK_Nodes', 'double', 0, 'b');
            else
                fprintf(fileID, '%0.15g %0.15g %0.15g\n', VTK_Nodes');
            end
            fprintf(fileID, '\n');

            [cellVerts, nv] = obj.local_clean_connectivity(verts_raw);
            totalInts = sum(double(nv) + 1);
            fprintf(fileID, 'CELLS %d %d\n', Grid.N, totalInts);

            if obj.isBinary
                for c = 1:Grid.N
                    v = int32(cellVerts{c} - 1);
                    fwrite(fileID, int32([nv(c), v']), 'int32', 0, 'b');
                end
            else
                for c = 1:Grid.N
                    fprintf(fileID, '%d ', nv(c));
                    fprintf(fileID, '%d ', cellVerts{c} - 1);
                    fprintf(fileID, '\n');
                end
            end
            fprintf(fileID, '\n');

            cellTypes = zeros(Grid.N,1,'int32');
            if isprop(Grid.mesh, 'element') && isprop(Grid.mesh.element, 'type')
                t = Grid.mesh.element.type(:);
                cellTypes(t==1) = int32(10);
                cellTypes(t==2) = int32(12); %
                cellTypes(t==3) = int32(13);
                cellTypes(t==4) = int32(14);
                miss = (cellTypes==0);
                if any(miss)
                    cellTypes(miss) = obj.local_infer_cell_types(nv(miss));
                end
            else
                cellTypes = obj.local_infer_cell_types(nv);
            end

            fprintf(fileID, 'CELL_TYPES %d\n', Grid.N);
            if obj.isBinary
                fwrite(fileID, cellTypes.', 'int32', 0, 'b');
            else
                fprintf(fileID, '%d\n', cellTypes);
            end
            fprintf(fileID, '\n\n');

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

            if ~isempty(Grid.ListOfFracturedReservoirCells)
                isFractured = zeros(Grid.N,1);
                isFractured(Grid.ListOfFracturedReservoirCells) = 1;
                obj.PrintScalar2VTK(fileID, isFractured, ' isFractured');
                fprintf(fileID, '\n');
            end

            if ~isempty(Grid.ListOfPerforatedCells)
                isPerforated = zeros(Grid.N,1);
                isPerforated(Grid.ListOfPerforatedCells) = 1;
                obj.PrintScalar2VTK(fileID, isPerforated, ' isPerforated');
                fprintf(fileID, '\n');
            end

            obj.PrintScalar2VTK(fileID, [1:Grid.N]', ' Index');         fprintf(fileID, '\n');
            obj.PrintScalar2VTK(fileID, Grid.Volume, ' Volume');        fprintf(fileID, '\n');
            obj.PrintScalar2VTK(fileID, Grid.ActiveTime, ' ACTIVETime');fprintf(fileID, '\n');
            obj.PrintScalar2VTK(fileID, Grid.Active, ' ACTIVEFine');    fprintf(fileID, '\n');

            fclose(fileID);

            % if obj.PlotInterfaces
            %     obj.PlotInternalFaces(Reservoir, Grid);
            % end
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

            VTK_Nodes = Grid.mesh.vertex.coord;
            if Grid.mesh.element.num ~= sum(Grid.Active)
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

            N_InternalFaces = Grid.mesh.face.inner.num;
            faceVertsCell  = Grid.mesh.face.inner.vertices;

        
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
            obj.PrintScalar2VTK(fileID, Grid.T, ' Trans');
            obj.PrintScalar2VTK(fileID, Grid.pEDFM_alpha_T, ' alpha');
            fprintf(fileID, '\n');

            fclose(fileID);
        end

        function PlotFractureSolution(obj, Fracture, Grid, f)
            fileID = fopen(strcat(obj.FileName, '_Fracture', num2str(f,'%03d'), '_', num2str(obj.VTKindex,'%04d'),'.vtk'), 'w');

            fprintf(fileID, '# vtk DataFile Version 2.0\n');
            fprintf(fileID, 'DARSim 2 Reservoir Simulator\n');
            if obj.isBinary
                fprintf(fileID, 'BINARY\n');
            else
                fprintf(fileID, 'ASCII\n');
            end
            fprintf(fileID, '\n');
            fprintf(fileID, 'DATASET UNSTRUCTURED_GRID\n');
            fprintf(fileID, '\n');

            verts_raw = Grid.fracmesh.elements.vertices;
            maxVid = obj.local_max_vertex_id(verts_raw);

            [VTK_Nodes, whichPts] = obj.local_pick_fracture_points(Grid, maxVid);

            if ~isempty(obj.ReferenceDepth)
                VTK_Nodes(:,3) = obj.ReferenceDepth - VTK_Nodes(:,3);
            end

            fprintf(fileID, 'POINTS    %d   double\n', size(VTK_Nodes,1));
            if obj.isBinary
                fwrite(fileID, VTK_Nodes', 'double', 0, 'b');
            else
                fprintf(fileID, '%0.15g %0.15g %0.15g\n', VTK_Nodes');
            end
            fprintf(fileID, '\n');

            [cellVerts, nv] = obj.local_clean_connectivity(verts_raw);
            totalInts = sum(double(nv) + 1);
            fprintf(fileID, 'CELLS %d %d\n', Grid.N, totalInts);

            if obj.isBinary
                for c = 1:Grid.N
                    v = int32(cellVerts{c} - 1);          % could be row or column
                    fwrite(fileID, int32([nv(c); v(:)]), 'int32', 0, 'b');
                end
            else
                for c = 1:Grid.N
                    fprintf(fileID, '%d ', nv(c));
                    fprintf(fileID, '%d ', cellVerts{c} - 1);
                    fprintf(fileID, '\n');
                end
            end
            fprintf(fileID, '\n');

            % CELL_TYPES：tri/quad/polygon
            celltype = zeros(Grid.N,1,'int32');
            for c = 1:Grid.N
                switch double(nv(c))
                    case 3
                        celltype(c) = int32(5); % VTK_TRIANGLE
                    case 4
                        celltype(c) = int32(9); % VTK_QUAD
                    otherwise
                        celltype(c) = int32(7); % VTK_POLYGON
                end
            end

            fprintf(fileID, 'CELL_TYPES %d\n', Grid.N);
            if obj.isBinary
                fwrite(fileID, celltype.', 'int32', 0, 'b');
            else
                fprintf(fileID, '%d\n', celltype);
            end
            fprintf(fileID, '\n\n');

            fprintf(fileID, 'CELL_DATA %d\n', Grid.N);
            fprintf(fileID, '\n');

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

    methods (Access=private)
        function maxVid = local_max_vertex_id(~, verts_raw)
            v = verts_raw(:);
            v = v(~isnan(v) & v>0);
            if isempty(v), maxVid = 0; else, maxVid = max(v); end
        end

        function [VTK_Nodes, whichPts] = local_pick_points(~, Grid, maxVid)
            VTK_Nodes = [];
            whichPts = '';

            if isprop(Grid,'GridCoords') && ~isempty(Grid.GridCoords) && size(Grid.GridCoords,1) >= maxVid
                VTK_Nodes = Grid.GridCoords;
                whichPts = 'Grid.GridCoords';
                return;
            end

            if isprop(Grid,'mesh') && isprop(Grid.mesh,'vertex') && isprop(Grid.mesh.vertex,'coord') && ...
               ~isempty(Grid.mesh.vertex.coord) && size(Grid.mesh.vertex.coord,1) >= maxVid
                VTK_Nodes = Grid.mesh.vertex.coord;
                whichPts = 'Grid.mesh.vertex.coord';
                return;
            end
            error('VTK POINTS do not match connectivity: max vertex id=%d, but available point arrays are too short.', maxVid);
        end

        function [VTK_Nodes, whichPts] = local_pick_fracture_points(~, Grid, maxVid)
            VTK_Nodes = [];
            whichPts = '';

            if isprop(Grid,'fracmesh') && isprop(Grid.fracmesh,'vertex') && isprop(Grid.fracmesh.vertex,'coord') && ...
               ~isempty(Grid.fracmesh.vertex.coord) && size(Grid.fracmesh.vertex.coord,1) >= maxVid
                VTK_Nodes = Grid.fracmesh.vertex.coord;
                whichPts = 'Grid.fracmesh.vertex.coord';
                return;
            end

            if isprop(Grid,'GridCoords') && ~isempty(Grid.GridCoords) && size(Grid.GridCoords,1) >= maxVid
                VTK_Nodes = Grid.GridCoords;
                whichPts = 'Grid.GridCoords';
                return;
            end
            error('VTK POINTS do not match connectivity: max vertex id=%d, but available point arrays are too short.', maxVid);
        end

        function [cellVerts, nv] = local_clean_connectivity(~, verts_raw)
            N = size(verts_raw,1);
            cellVerts = cell(N,1);
            nv = zeros(N,1,'int32');

            for c = 1:N
                v = verts_raw(c,:);
                v = v(:)';
                v = v(~isnan(v) & v>0);

                if ~isempty(v)
                   
                    v = v([true, diff(v)~=0]);
                end

                cellVerts{c} = v;
                nv(c) = int32(numel(v));
            end
        end

        function cellVertsNew = wedge_order_correction( ...
                ~, VTK_Nodes, cellTypes, cellVertsOld)
        
            cellVertsNew = cellVertsOld;
        
            for c = 1:numel(cellVertsOld)
        
                if cellTypes(c) ~= 13
                    continue
                end
        
                ids = cellVertsOld{c};
                ids = ids(isfinite(ids) & ids > 0);
        
                if numel(ids) ~= 6
                    error('Cell %d is marked as a wedge but has %d vertices.', ...
                          c, numel(ids));
                end
        
                xyz = VTK_Nodes(ids,:);
        
                xyzCentered = xyz - mean(xyz,1);
                [~,~,V] = svd(xyzCentered, 0);
        
                bestGap = -inf;
        
                for k = 1:3
                    prismAxisCandidate = V(:,k);
                    projection = xyzCentered * prismAxisCandidate;
        
                    [projectionSorted, order] = sort(projection);
                    gap = projectionSorted(4) - projectionSorted(3);
        
                    if gap > bestGap
                        bestGap = gap;
                        prismAxis = prismAxisCandidate;
                        sortedOrder = order;
                    end
                end
        
                bottom = sortedOrder(1:3);
                top    = sortedOrder(4:6);
        
                if mean(xyz(bottom,3)) > mean(xyz(top,3))
                    temp = bottom;
                    bottom = top;
                    top = temp;
                    prismAxis = -prismAxis;
                end
        
                bottomProjected = xyz(bottom,:) ...
                    - (xyz(bottom,:) * prismAxis) * prismAxis';
        
                topProjected = xyz(top,:) ...
                    - (xyz(top,:) * prismAxis) * prismAxis';
        
                permutations = perms(1:3);
                bestDistance = inf;
                bestTop = top;
        
                for p = 1:size(permutations,1)
                    candidatePermutation = permutations(p,:);
                    candidateTop = top(candidatePermutation);
        
                    d = topProjected(candidatePermutation,:) ...
                        - bottomProjected;
        
                    totalDistance = sum(vecnorm(d,2,2));
        
                    if totalDistance < bestDistance
                        bestDistance = totalDistance;
                        bestTop = candidateTop;
                    end
                end
        
                top = bestTop;
        
                p1 = xyz(bottom(1),:);
                p2 = xyz(bottom(2),:);
                p3 = xyz(bottom(3),:);
                p4 = xyz(top(1),:);
        
                J = [p2-p1; p3-p1; p4-p1]';
        
                if det(J) < 0
                    bottom([2 3]) = bottom([3 2]);
                    top([2 3]) = top([3 2]);
                end
        
                orderedLocalIDs = [bottom(:); top(:)];
                cellVertsNew{c} = ids(orderedLocalIDs)';
            end
        end

        function cellTypes = local_infer_cell_types(~, nv)
            n = numel(nv);
            cellTypes = zeros(n,1,'int32');
            for i = 1:n
                switch double(nv(i))
                    case 4
                        cellTypes(i) = int32(10); 
                    case 5
                        cellTypes(i) = int32(14); % pyramid
                    case 6
                        cellTypes(i) = int32(13); % wedge
                    case 8
                        cellTypes(i) = int32(12); % hex
                    otherwise
                        error('Unsupported cell with %d vertices. If this is a 3D polyhedron, need VTK_POLYHEDRON(42).', nv(i));
                end
            end
        end
    end
end
