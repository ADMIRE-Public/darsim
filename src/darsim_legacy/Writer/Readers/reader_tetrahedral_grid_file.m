%  Read the tetrahedral grid file
%--------------------------------------------------------------------------
% DARSim 2 Reservoir Simulator
% Author: Túlio de M. Cavalcante
% UFPE
%--------------------------------------------------------------------------
classdef reader_tetrahedral_grid_file < handle
    properties
        Directory
        File
    end
    methods
        function obj = reader_tetrahedral_grid_file(dir, file)
            if nargin == 1
                obj.Directory = [];
                obj.File = dir;
            else
                obj.Directory = dir;
                obj.File = strcat(dir,'/',file);
            end
        end
        function TetrahedralGridData = ReadMSH(obj)   
            %"coord" matrix - The coordinates (x,y,z) of each point in a carts: 
            [ coord, nnode ] = getcoord(obj);

            %"elem" matrix - The points that constitute each element:
            [ elem, bedgeref ] = getelem( obj, nnode );
            
            %"bedge" and "inedge" - It get inform. about boundary and internal edges:
            [ TetrahedralGridData ] = getinfo( obj, coord, elem, bedgeref );
        end
        function [ coord, nnode ] = getcoord( obj )
            % write filepath
            filepath = obj.File;
            %Open the *.msh file
            readmsh = fopen(filepath);
            %"nnode" is the number of nodes in the discrete domain
            getmshdata = textscan(readmsh,'%u',1,'HeaderLines',4);
            %Attribute the data to "nnode"
            nnode = getmshdata{1};

            %"coord" is a matrix which contain the coordinate of each point of domain
            getmshdata = textscan(readmsh,'%*u %f64 %f64 %f64',nnode,'HeaderLines',1);
            %Fill the matrix "coord" with the x, y and z coordinates.
            coord = cell2mat(getmshdata);

            %Close the *.msh file
            fclose(readmsh);
        end
        function [ elem, bedgeref ] = getelem( obj, nnode )
            % write filepath
            filepath = obj.File;
            %Open the *.msh file
            readmsh = fopen(filepath);        
            %"nent" is the number of all entities understud as element in *.msh file
            getmshdata = textscan(readmsh,'%u',1,'HeaderLines',7 + nnode);
            %Attribute the data to "nnode"
            nent = getmshdata{1};
            %Swept the "elements" of "*.msh" in order to fill the parameters above. 
            getmshdata = textscan(readmsh,'%*n %n %*n %n %*[^\n]',nent);
            %Attribute the data to "auxmat"
            auxmat = cell2mat(getmshdata);
            %Fill "entitype"
            entitype = auxmat(:,1);
            %"nbe" verifies how many entities are edges which constitute the 
            nbedge = sum(entitype == 2);
            nelem = sum(entitype == 4);
            %Close the *.msh file
            fclose(readmsh);

            %---------------------------
            %Define the matrix ("elem"):
            %Open again the *.msh file
            readmsh = fopen(filepath);
            %Create and initialize the parameter "elem"
            elem = zeros(nelem,4);

            if nbedge > 0
                H1 = 8 + nnode + sum(entitype == 15) + sum(entitype == 1);
                getmshdata = textscan(readmsh,'%*u %*u %*u %u %*u %u %u %u',nbedge,...
                    'HeaderLines',H1);
                bedgeaux = cell2mat(getmshdata);
                bedgeref(1:nbedge,5) = bedgeaux(:,1); 
                bedgeref(1:nbedge,1:3) = bedgeaux(:,2:4); 
            end  %End of second IF

            if nelem > 0
                getmshdata = textscan(readmsh,'%*u %*u %*u %u %*u %u %u %u %u',nelem,...
                    'HeaderLines',1);
                elemaux = cell2mat(getmshdata);
                %Attribute to "elem" (another column)
                elem(1:nelem,1:4) = elemaux(:,2:5); 
            end  %End of second IF

            %Close the *.msh file
            fclose(readmsh);
        end
        function [ TetrahedralGridData ] = getinfo( obj, coord, elem, bedgeref )
            %
            fprintf('--> Creating data structure for vertices... ')
            
            TetrahedralGridData.Vertex.Coordinates = coord;
            
            TetrahedralGridData.Cell.Vertices = elem;
      
            fprintf('Done!\n')
            %--------------------------------------------------------------
            fprintf('--> Creating data structure for cells... ')
            
            centelem = reshape(mean(reshape(coord(elem(:,1:4),:),[size(elem,1),4,3]),2),[size(elem,1),3,1]);
            
            TetrahedralGridData.Cell.Centroid = centelem;
            %--------------------------------------------------------------
            
            %--------------------------------------------------------------
            I = coord(elem(:,1),:); J = coord(elem(:,2),:);
            K = coord(elem(:,3),:); Q = coord(elem(:,4),:);

            IQ = Q - I; IJ = J - I; IK = K - I;

            elemvolume = (1/6)*abs(sum(cross(IJ,IK).*IQ,2));
            
            TetrahedralGridData.Cell.Volume = elemvolume;
            
            fprintf('Done!\n')
            %--------------------------------------------------------------
            
            %--------------------------------------------------------------
            fprintf('--> Creating data structure for faces... ')
            
            elemaux = sort(elem(:,1:4),2);
            
            inedgeaux = [elemaux(:,1) elemaux(:,2) elemaux(:,3) [1:size(elem,1)]';...
                         elemaux(:,1) elemaux(:,2) elemaux(:,4) [1:size(elem,1)]';...
                         elemaux(:,1) elemaux(:,3) elemaux(:,4) [1:size(elem,1)]';...
                         elemaux(:,2) elemaux(:,3) elemaux(:,4) [1:size(elem,1)]'];
                
            inedgeaux = sortrows(inedgeaux,[1 2 3 4]); % Ordinate the faces matrix in ascending way of vertices;
            
            [ ~, idx ] = unique(inedgeaux(:,1:3),'rows');
            
            repet = setdiff(1:size(inedgeaux,1),idx)'; % Get the repeated faces;
            
            edges = zeros(size(inedgeaux,1),5);
            
            edges(idx,1:4) = inedgeaux(idx,1:4);
            
            edges(repet-1,5) = inedgeaux(repet,4); % Merge the equal faces;
            
            inidx = find(edges(:,5)~=0); bdidx = find(edges(:,5)==0 & edges(:,1)~=0); % Index for inner and boundary faces;
            
            inedge = edges(inidx,:); bedge = edges(bdidx,:); bedge(:,5) = 201; % Null flow flag;
                       
            I = coord(inedge(:,1),:); J = coord(inedge(:,2),:);
            K = coord(inedge(:,3),:); L = centelem(inedge(:,4),:);
            
            JI = I - J; JK = K - J; LJ = J - L;
            
            innormals = 0.5*(cross(JI,JK)); trocaLR = -0.5*(sign(sum(LJ.*innormals,2)) - 1);
            
            inedge(trocaLR==1,4:5) = [ inedge(trocaLR==1,5) inedge(trocaLR==1,4) ];
            
            I = coord(bedge(:,1),:); J = coord(bedge(:,2),:);
            K = coord(bedge(:,3),:); L = centelem(bedge(:,4),:);
            
            JI = I - J; JK = K - J; LJ = J - L;
            
            bounnormals = 0.5*(cross(JI,JK)); trocaJK = -0.5*(sign(sum(LJ.*bounnormals,2)) - 1);
            
            bedge(trocaJK==1,2:3) = [ bedge(trocaJK==1,3) bedge(trocaJK==1,2) ];
            
            bounnormals(trocaJK==1,:) = -bounnormals(trocaJK==1,:); 
            
            TetrahedralGridData.Face.Inner.Normals = innormals;
            
            TetrahedralGridData.Face.Bound.Normals = bounnormals;
            
            TetrahedralGridData.Face.Inner.Vertices = inedge(:,1:3);
            
            TetrahedralGridData.Face.Inner.MontCell = inedge(:,4);
            
            TetrahedralGridData.Face.Inner.JusCell = inedge(:,5);
            
            TetrahedralGridData.Face.Bound.Vertices = bedge(:,1:3);
            
            TetrahedralGridData.Face.Bound.MontCell = bedge(:,4);
            
            TetrahedralGridData.Face.Bound.Flag = bedge(:,5);
            
            fprintf('Done!\n')
            %--------------------------------------------------------------
            
            fprintf('--> Creating data structure for connectivities... ')
            %--------------------------------------------------------------
            eleminfacematrix1 = sparse(size(elem,1),size(inedge,1)); 
            
            eleminfacematrix2 = sparse(size(elem,1),size(inedge,1));
    
            eleminfacematrix1(sub2ind(size(eleminfacematrix1),inedge(:,4),[1:size(inedge,1)]')) = 1; 
            
            eleminfacematrix2(sub2ind(size(eleminfacematrix2),inedge(:,5),[1:size(inedge,1)]')) = 1;
            
            eleminfacematrix = eleminfacematrix1 + eleminfacematrix2;
            
            TetrahedralGridData.Cell.CellInnerFaceMatrix = eleminfacematrix;
            %--------------------------------------------------------------
            
            %--------------------------------------------------------------
            eleminfacematrix3 = sparse(size(elem,1),size(inedge,1));
            
            eleminfacematrix3(sub2ind(size(eleminfacematrix3),inedge(:,5),[1:size(inedge,1)]')) = -1;
            
            infaceelemmatrix = ( eleminfacematrix1 + eleminfacematrix3 )';
            
            TetrahedralGridData.Face.Inner.InnerFaceCellMatrix = infaceelemmatrix;
            %--------------------------------------------------------------
            
            %--------------------------------------------------------------
            elembdfacematrix = sparse(size(elem,1),size(bedge,1));
            
            elembdfacematrix(sub2ind(size(elembdfacematrix),bedge(:,4),[1:size(bedge,1)]')) = 1;
            
            TetrahedralGridData.Cell.CellBoundFaceMatrix = elembdfacematrix;
            %--------------------------------------------------------------
            
            %--------------------------------------------------------------
            vertinfacematrix1 = sparse(size(coord,1),size(inedge,1)); 
            
            vertinfacematrix2 = sparse(size(coord,1),size(inedge,1)); 
            
            vertinfacematrix3 = sparse(size(coord,1),size(inedge,1));         
            
            vertinfacematrix1(sub2ind(size(vertinfacematrix1),inedge(:,1),[1:size(inedge,1)]')) = 1; 
            
            vertinfacematrix2(sub2ind(size(vertinfacematrix2),inedge(:,2),[1:size(inedge,1)]')) = 1;
            
            vertinfacematrix3(sub2ind(size(vertinfacematrix3),inedge(:,3),[1:size(inedge,1)]')) = 1;
            
            vertinfacematrix = vertinfacematrix1 + vertinfacematrix2 + vertinfacematrix3;
            
            TetrahedralGridData.Vertex.VertInnerFaceMatrix = vertinfacematrix;
            
            TetrahedralGridData.Vertex.IVertInnerFaceMatrix = vertinfacematrix1;
            
            TetrahedralGridData.Vertex.JVertInnerFaceMatrix = vertinfacematrix2;
            
            TetrahedralGridData.Vertex.KVertInnerFaceMatrix = vertinfacematrix3;
            %--------------------------------------------------------------
            
            %--------------------------------------------------------------
            vertbdfacematrix1 = sparse(size(coord,1),size(bedge,1));
            
            vertbdfacematrix2 = sparse(size(coord,1),size(bedge,1));
            
            vertbdfacematrix3 = sparse(size(coord,1),size(bedge,1));
            
            vertbdfacematrix1(sub2ind(size(vertbdfacematrix1),bedge(:,1),[1:size(bedge,1)]')) = 1;
            
            vertbdfacematrix2(sub2ind(size(vertbdfacematrix2),bedge(:,2),[1:size(bedge,1)]')) = 1;
            
            vertbdfacematrix3(sub2ind(size(vertbdfacematrix3),bedge(:,3),[1:size(bedge,1)]')) = 1;
            
            vertbdfacematrix = vertbdfacematrix1 + vertbdfacematrix2 + vertbdfacematrix3;
            
            TetrahedralGridData.Vertex.VertBoundFaceMatrix = vertbdfacematrix;
            %--------------------------------------------------------------
            
            %--------------------------------------------------------------
            vertelemmatrix1 = sparse(size(coord,1),size(elem,1));
            
            vertelemmatrix2 = sparse(size(coord,1),size(elem,1));
            
            vertelemmatrix3 = sparse(size(coord,1),size(elem,1));
            
            vertelemmatrix4 = sparse(size(coord,1),size(elem,1));
            
            vertelemmatrix1(sub2ind(size(vertelemmatrix1),elem(:,1),[1:size(elem,1)]')) = 1;
            
            vertelemmatrix2(sub2ind(size(vertelemmatrix2),elem(:,2),[1:size(elem,1)]')) = 1;
            
            vertelemmatrix3(sub2ind(size(vertelemmatrix3),elem(:,3),[1:size(elem,1)]')) = 1;
            
            vertelemmatrix4(sub2ind(size(vertelemmatrix4),elem(:,4),[1:size(elem,1)]')) = 1;
            
            vertelemmatrix = vertelemmatrix1 + vertelemmatrix2 + vertelemmatrix3 + vertelemmatrix4;
            
            TetrahedralGridData.Vertex.VertCellMatrix = vertelemmatrix;
            %--------------------------------------------------------------
            fprintf('Done!\n')
        end     
    end
end
