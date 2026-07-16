% Simulation class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Mousa HosseiniMehr
%TU Delft
%Created: 12 July 2016
%Last modified: 17 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef geometry_CartesianGrid_FracGen < geometry_FracGen
    properties
    end
    methods
        function Compute_EDFM_Fracture_Matrix_Connectivities(obj)
            Fracture = obj.FracturesGrid.Fracture;
            %% Assigning reservoir properties
            Epsilon = obj.ReservoirGrid.Epsilon;
            % Number of reservoir grid cells in each direction
            Nx = obj.ReservoirGrid.Nx;
            Ny = obj.ReservoirGrid.Ny;
            Nz = obj.ReservoirGrid.Nz;
            % Coordinates of reservoir cell centers in each direction
            Xcm = obj.ReservoirGrid.Xcm;
            Ycm = obj.ReservoirGrid.Ycm;
            Zcm = obj.ReservoirGrid.Zcm;
            % Coordinates of reservoir grid nodes in each direction
            Xim = obj.ReservoirGrid.Xim;
            Yim = obj.ReservoirGrid.Yim;
            Zim = obj.ReservoirGrid.Zim;
            %% Intersections of fracture cells with each matrix cell
            fprintf('Obtaining fracture - matrix connectivities:\n');
            fprintf('---> Fracture ');
            for f = 1 : length(Fracture)
                if (f>1),  fprintf(repmat('\b', 1, 5+29));  end
                fprintf('%03d/%03d',f,length(Fracture));
                fprintf(' ---> Grid cell ');
                for j_f = 1 : Fracture(f).N_Width_AD
                    for i_f = 1 : Fracture(f).N_Length_AB
                        If = Index_Local_To_Global_3D( Fracture(f).N_Length_AB , Fracture(f).N_Width_AD , 1 , i_f , j_f , 1 );
                        if (If>1),  fprintf(repmat('\b', 1, 11));  end
                        fprintf('%05d/%05d',If,Fracture(f).N_Length_AB*Fracture(f).N_Width_AD);
                        Index_matIntersect = 0;
                        
                        % The matrix cell closest to fracture cell for the first intersection check
                        [ ~ , i_1st ] = min( ( Fracture(f).CellCenterCoordsV1(i_f,j_f,1) - Xcm ).^2 );
                        [ ~ , j_1st ] = min( ( Fracture(f).CellCenterCoordsV1(i_f,j_f,2) - Ycm ).^2 );
                        [ ~ , k_1st ] = min( ( Fracture(f).CellCenterCoordsV1(i_f,j_f,3) - Zcm ).^2 );
                        ImList  = Index_Local_To_Global_3D( Nx,Ny,Nz , i_1st , j_1st , k_1st );
                        Im_count = 1;
                        
                        while Im_count <= length(ImList)
                            Im = ImList(Im_count);
                            Im_count = Im_count +1;
                            
                            % Retrieving i,j,k from the "Im" index
                            [i,j,k] = Index_Global_To_Local_3D(Nx, Ny, Nz, Im);

                            % Corner Points of the fracture cell
                            PointA = [ Fracture(f).GridCoords(i_f  ,j_f  ,1) , Fracture(f).GridCoords(i_f  ,j_f  ,2) , Fracture(f).GridCoords(i_f  ,j_f  ,3) ];
                            PointB = [ Fracture(f).GridCoords(i_f  ,j_f+1,1) , Fracture(f).GridCoords(i_f  ,j_f+1,2) , Fracture(f).GridCoords(i_f  ,j_f+1,3) ];
                            PointC = [ Fracture(f).GridCoords(i_f+1,j_f+1,1) , Fracture(f).GridCoords(i_f+1,j_f+1,2) , Fracture(f).GridCoords(i_f+1,j_f+1,3) ];
                            PointD = [ Fracture(f).GridCoords(i_f+1,j_f  ,1) , Fracture(f).GridCoords(i_f+1,j_f  ,2) , Fracture(f).GridCoords(i_f+1,j_f  ,3) ];
                            Plane_fracCell = tetragon_DARSim(PointA,PointB,PointC,PointD);
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % Obtaining the intersection points between the matrix cell and fracture cell
                            PointA = [ Xim(i)   , Yim(j)   , Zim(k)   ];
                            PointB = [ Xim(i+1) , Yim(j+1) , Zim(k+1) ];
                            CubeMatCell = cube_DARSim(PointA,PointB);
                            [intersectCoord, areCoplanar] = CubeMatCell.Obtain_Cube_Tetragon_Intersection(Plane_fracCell,Epsilon);
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            if isempty(intersectCoord)
                                if Im_count == 2
                                    % If it is the first try and no intersection occurs between the fracture cell
                                    % and this matrix cell, add the neighboring matrix cells to the ImList for intersection check.
                                    indNeighbors = sort(obj.ReservoirGrid.Neighbours{Im});
                                    ImList = union(ImList, indNeighbors, 'stable');
                                end
                                continue;
                            else
                                % There was already an intersection before the fracture cell and this matrix cell.
                                % Therefore, add the neighboring matrix cells to the ImList for intersection check
                                % as the neighbors may intersect as well.
                                indNeighbors = sort(obj.ReservoirGrid.Neighbours{Im});  
                                ImList = union(ImList, indNeighbors, 'stable');
                            end
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % Writing the final intersection Points into the fracture data structure
                            Index_matIntersect = Index_matIntersect + 1;
                            Fracture(f).intersectCoord_rock{If}{Index_matIntersect,1} = Im;
                            Fracture(f).       Af_rock_EDFM{If}(Index_matIntersect,1) = Im;
                            Fracture(f).  avgDist_rock_EDFM{If}(Index_matIntersect,1) = Im;
                            Fracture(f).       CI_rock_EDFM{If}(Index_matIntersect,1) = Im;
                            Fracture(f).      Coplanar_EDFM{If}(Index_matIntersect,1) = Im;
                            
                            Fracture(f).intersectCoord_rock{If}{Index_matIntersect,2} = intersectCoord;
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % Calculating the area fraction
                            Fracture(f).Af_rock_EDFM{If}(Index_matIntersect,2) = 0;
                            for n = 2 : size( intersectCoord , 1 )  - 1
                                Fracture(f).Af_rock_EDFM{If}(Index_matIntersect,2) = Fracture(f).Af_rock_EDFM{If}(Index_matIntersect,2) + ...
                                    Triangle_Area_3D( intersectCoord(1,:) , intersectCoord(n,:) , intersectCoord(n+1,:) );
                            end
                            
                            % Doubling the area fraction if the cell is not coplanar to any matrix cube edges
                            if areCoplanar == 0
                                Fracture(f).Af_rock_EDFM{If}(Index_matIntersect,2) = Fracture(f).Af_rock_EDFM{If}(Index_matIntersect,2) * 2;
                            end
                            Fracture(f).Coplanar_EDFM{If}(Index_matIntersect,2) = areCoplanar;
                            Fracture(f).Af_rock_EDFM_sum = Fracture(f).Af_rock_EDFM_sum + Fracture(f).Af_rock_EDFM{If}(Index_matIntersect,2);
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % Obtaining the average distance between the matrix cube cell and the fracture plate cell
                            Fracture(f).avgDist_rock_EDFM{If}(Index_matIntersect,2) = 0;
                            Refinement = 6;
                            Fracture(f).avgDist_rock_EDFM{If}(Index_matIntersect,2) = CubeMatCell.Obtain_Average_Distance_Plane_From_Cube(Fracture(f), Refinement);
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % Obtaining the geometrical transmissibility between the matrix cube cell and the fracture plate cell ( CI = Af/<d> )
                            Fracture(f).CI_rock_EDFM{If}(Index_matIntersect,2) = ...
                            Fracture(f).Af_rock_EDFM{If}(Index_matIntersect,2) / ...
                            Fracture(f).avgDist_rock_EDFM{If}(Index_matIntersect,2);
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % Adding the matrix cell index to the list of
                            % all overlapped matrix cells inside the
                            % "reservoir" class and "fracture" class
                            if ~ismember(Im,obj.ReservoirGrid.OvelapList)
                                obj.ReservoirGrid.OvelapList = [obj.ReservoirGrid.OvelapList;Im];
                            end
                            if ~ismember(Im,Fracture(f).MatrixOverlapList)
                                Fracture(f).MatrixOverlapList = [Fracture(f).MatrixOverlapList;Im];
                            end
                            
                        end % End of while-loop
                        
                        % Obtaining the overlapped cells due to the intersection line between fractures f and g
                        max_dist = sqrt( Fracture(f).D_Length_AB^2 + Fracture(f).D_Width_AD^2 )/2;
                        for g = 1 : length(Fracture)
                            if ~isempty(Fracture(f).intersectCoord_fracturePlate{g})
                                p1 = Fracture(f).intersectCoord_fracturePlate{g}(1,:);
                                p2 = Fracture(f).intersectCoord_fracturePlate{g}(2,:);
                                Fracture(f).Overlap_frac{g}((j_f-1)*Fracture(f).N_Length_AB+i_f) = ( norm( cross(Fracture(f).CellCenterCoordsV2((j_f-1)*Fracture(f).N_Length_AB+i_f,:)-p1 , p2-p1) ) / norm(p2-p1) ) <= max_dist;
                            end
                        end
                        Fracture(f).     Af_rock_EDFM{If} = sortrows(Fracture(f).     Af_rock_EDFM{If});
                        Fracture(f).avgDist_rock_EDFM{If} = sortrows(Fracture(f).avgDist_rock_EDFM{If});
                        Fracture(f).     CI_rock_EDFM{If} = sortrows(Fracture(f).     CI_rock_EDFM{If});
                        Fracture(f).     Coplanar_EDFM{If} = sortrows(Fracture(f).    Coplanar_EDFM{If});
                    end % End of fractre 1st for-loop
                end % End of fractre 2nd for-loop
                Fracture(f).MatrixOverlapList = sort(Fracture(f).MatrixOverlapList);
            end % End of main fracture for-loop
            obj.ReservoirGrid.OvelapList = sort(obj.ReservoirGrid.OvelapList);
            
            fprintf('\n---------------------------------------------------------\n');
        end
        function Add_pEDFM_Connectivities(obj)
            Fracture = obj.FracturesGrid.Fracture;
            Epsilon = obj.ReservoirGrid.Epsilon;
            %% Assigning reservoir properties
            % Number of reservoir grid cells in each direction
            Nx = obj.ReservoirGrid.Nx;
            Ny = obj.ReservoirGrid.Ny;
            Nz = obj.ReservoirGrid.Nz;
            % Coordinates of reservoir cell centers in each direction
            Xcm = obj.ReservoirGrid.Xcm;
            Ycm = obj.ReservoirGrid.Ycm;
            Zcm = obj.ReservoirGrid.Zcm;
            % Coordinates of reservoir grid nodes in each direction
            Xim = obj.ReservoirGrid.Xim;
            Yim = obj.ReservoirGrid.Yim;
            Zim = obj.ReservoirGrid.Zim;
            %% Adding pEDFM Connectivities
            fprintf('Adding pEDFM Connectivities:\n');
            fprintf('---> Fracture ');
            
            % Looping through fractures
            for f = 1 : length(Fracture)
                if (f>1),  fprintf(repmat('\b', 1, 5+27));  end
                fprintf('%02d/%02d',f,length(Fracture));
                fprintf(' ---> Grid cell ');
                for If = 1 : Fracture(f).N_Total
                    if (If>1),  fprintf(repmat('\b', 1, 11));  end
                    fprintf('%05d/%05d',If,Fracture(f).N_Total);
                    
                    %% 1- Obtaining pEDFM connectivities of this frature with matrix
                    Af = sortrows(Fracture(f).Af_rock_EDFM{If});
                    areCoplanar = sortrows(Fracture(f).Coplanar_EDFM{If});
                    if isempty(Af),  continue;  end
                    ImList = Af(:,1);
                    % Initializing the necessary variable to collect the conectivities data
                    counter = 0;
                    pEDFM_List = [];
                    pEDFM_Af = [];
                    pEDFM_AvgDist = [];
                    
                    % Looping over the overlapped matrix cells
                    for m = 1 : length(ImList)
                        Im_This = ImList(m);
                        [ i_This , j_This , k_This ] = Index_Global_To_Local_3D(Nx, Ny, Nz, Im_This);
                        CellCenterCoord_This = [ Xcm(i_This) , Ycm(j_This) , Zcm(k_This) ];
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % Finding if the fracture plate crosses between the matrix cell and any of its neighbors.
                        % If it crosses, the connectivity between these two matrix cells is affected and the
                        % transmissiblity assigned to that interface is modified by "alpha" factor of pEDFM.
                        indNeighbors = obj.ReservoirGrid.Neighbours{Im_This};
                        for n = 1 : length(indNeighbors)
                            Im_Other = indNeighbors(n);
                            [ i_Other , j_Other , k_Other ] = Index_Global_To_Local_3D(Nx, Ny, Nz, Im_Other);
                            CellCenterCoord_Other = [ Xcm(i_Other) , Ycm(j_Other) , Zcm(k_Other) ];
                            LineSegmentConnectingCellCenters = lineSegment_DARSim( CellCenterCoord_This , CellCenterCoord_Other );
                            LineSegmentConnectingCellCenters.AddEquation(CellCenterCoord_This);
                            
                            % Check if the fracture plane crosses
                            [Geostatus, IntersectPoint] = Fracture(f).Obtain_Polygon_LineSegment_Intersection(LineSegmentConnectingCellCenters, Epsilon);
                            
                            if Geostatus.haveIntersect == 1 && ~isempty(IntersectPoint)
                                % [~ , Dimension] = setdiff([i_This,j_This,k_This] , [i_Other,j_Other,k_Other] ); % setdiff is not working properily.
                                Dimension = find( [i_This,j_This,k_This] ~= [i_Other,j_Other,k_Other] );
                                if isempty(Dimension) ,  continue;  end
                                if length(Dimension)>1,  error('Only one of the i,j,k indices can be different.');  end
                                
                                % Check if the intersected cell is already on the list or not
                                if ~ismember(Im_Other,pEDFM_List)
                                    pEDFM_List = [pEDFM_List;Im_Other];
                                    pEDFM_Af = [pEDFM_Af;0];
                                    pEDFM_AvgDist = [pEDFM_AvgDist;0];
                                    counter = counter+1;
                                    Index = counter;
                                else
                                    Index = find(pEDFM_List==Im_Other);
                                end
                                
                                % Calculating the projected area fraction Af_proj:
                                % The Cosine of the projection angle can be calculated as dot product between the
                                % normal vectors of fracture plate and projection face.
                                % Note: In EDFM, Af is doubled for fracture cells
                                % that are not coplanar with matrix mesh
                                % face. This is not done for pEDFM. 
                                % If the fracture is coplanar with matrix mesh face, 
                                % set pEDFM area to zero to ignore pEDFM connection altogether. 
                                switch Dimension
                                    case(1) % x direction
                                        if areCoplanar(m,2) == 1
                                            Af_proj = (Af(m,2)) * abs(dot(Fracture(f).nVec,[1,0,0]));
                                        else 
                                            Af_proj = (Af(m,2) / 2) * abs(dot(Fracture(f).nVec,[1,0,0]));
                                        end
                                    case(2) % y direction
                                        if areCoplanar(m,2) == 1
                                            Af_proj = (Af(m,2)) * abs(dot(Fracture(f).nVec,[0,1,0]));
                                        else
                                            Af_proj = (Af(m,2) / 2) * abs(dot(Fracture(f).nVec,[0,1,0]));
                                        end
                                    case(3) % z direction
                                        if areCoplanar(m,2) == 1
                                            Af_proj = (Af(m,2)) * abs(dot(Fracture(f).nVec,[0,0,1]));
                                        else
                                            Af_proj = (Af(m,2) / 2) * abs(dot(Fracture(f).nVec,[0,0,1]));
                                        end
                                end
                                pEDFM_Af(Index) = pEDFM_Af(Index) + (1-areCoplanar(m,2))*Af_proj;
                                
                                % Calculating average distance AvgDist_proj:
                                PointA = [ Xim(i_Other)   , Yim(j_Other)   , Zim(k_Other)   ];
                                PointB = [ Xim(i_Other+1) , Yim(j_Other+1) , Zim(k_Other+1) ];
                                CubeMatCell = cube_DARSim(PointA,PointB);
                                Refinement = 6;
                                AvgDist_proj = CubeMatCell.Obtain_Average_Distance_Plane_From_Cube(Fracture(f), Refinement);
                                pEDFM_AvgDist(Index) = AvgDist_proj;
                                
                                % Measuring the alpha:
                                % alpha is a the fraction of blockage (from 0 to 1) between the
                                % neighboring matrix cells due to existence of fracture plates when
                                % considering pEDFM method. the alpha=0 means no blockage and
                                % alpha=1 means total blockage.
                                switch Dimension
                                    case(1)
                                        alpha = (Af_proj) / (obj.ReservoirGrid.dy * obj.ReservoirGrid.dz);
                                        obj.ReservoirGrid.pEDFM_alpha_Tx( max(i_This,i_Other) , j_This , k_This ) =      obj.ReservoirGrid.pEDFM_alpha_Tx( max(i_This,i_Other) , j_This , k_This ) + alpha;
                                        obj.ReservoirGrid.pEDFM_alpha_Tx( max(i_This,i_Other) , j_This , k_This ) = min( obj.ReservoirGrid.pEDFM_alpha_Tx( max(i_This,i_Other) , j_This , k_This ) , 1 );
                                    case(2)
                                        alpha = (Af_proj) / (obj.ReservoirGrid.dx * obj.ReservoirGrid.dz);
                                        obj.ReservoirGrid.pEDFM_alpha_Ty( i_This , max(j_This,j_Other) , k_This ) =      obj.ReservoirGrid.pEDFM_alpha_Ty( i_This , max(j_This,j_Other) , k_This ) + alpha;
                                        obj.ReservoirGrid.pEDFM_alpha_Ty( i_This , max(j_This,j_Other) , k_This ) = min( obj.ReservoirGrid.pEDFM_alpha_Ty( i_This , max(j_This,j_Other) , k_This ) , 1 );
                                    case(3)
                                        alpha = (Af_proj) / (obj.ReservoirGrid.dx * obj.ReservoirGrid.dy);
                                        obj.ReservoirGrid.pEDFM_alpha_Tz( i_This , j_This , max(k_This,k_Other) ) =      obj.ReservoirGrid.pEDFM_alpha_Tz( i_This , j_This , max(k_This,k_Other) ) + alpha;
                                        obj.ReservoirGrid.pEDFM_alpha_Tz( i_This , j_This , max(k_This,k_Other) ) = min( obj.ReservoirGrid.pEDFM_alpha_Tz( i_This , j_This , max(k_This,k_Other) ) , 1 );
                                end
                            end
                        end
                    end
                    
                    % Adding the pEDFM information to the Fracture objects
                    pEDFM_Temp = [pEDFM_List,pEDFM_Af,pEDFM_AvgDist,pEDFM_Af./pEDFM_AvgDist];
                    if ~isempty(pEDFM_Temp)
                        pEDFM_Temp = sortrows(pEDFM_Temp,1);
                        Fracture(f).     Af_rock_pEDFM{If} = pEDFM_Temp(:,[1,2]);
                        Fracture(f).avgDist_rock_pEDFM{If} = pEDFM_Temp(:,[1,3]);
                        Fracture(f).     CI_rock_pEDFM{If} = pEDFM_Temp(:,[1,4]);
                    end
                    
                    %% 2- Obtaining pEDFM connectivities of this frature with other fractures
                    if Fracture(f).NumOf_fracConn_EDFM(If) > 0
                        for g = 1 : length(Fracture)
                            if g==f,  continue;  end
                            if isempty(Fracture(f).Af_frac_EDFM{g}),  continue;  end
                            if ~isempty(Fracture(f).Af_frac_EDFM{g}{If})
                                if isempty(Fracture(f).Af_frac_pEDFM{g})
                                    Fracture(f).     Af_frac_pEDFM{g} = cell(Fracture(f).N_Total,1);
                                    Fracture(f).avgDist_frac_pEDFM{g} = cell(Fracture(f).N_Total,1);
                                    Fracture(f).     CI_frac_pEDFM{g} = cell(Fracture(f).N_Total,1);
                                end
                                Aff = sortrows(Fracture(f).Af_frac_EDFM{g}{If});
                                Ig_List = Aff(:,1);
                                counter = 0;
                                pEDFM_List = [];
                                pEDFM_Af = [];
                                pEDFM_AvgDist = [];
                                IntersectLine = lineSegment_DARSim( Fracture(f).intersectCoord_fracturePlate{g}(2,:),Fracture(f).intersectCoord_fracturePlate{g}(1,:) );
                                for m = 1 : length(Ig_List)
                                    Ig_This = Ig_List(m);
                                    [ i_This , j_This , ~ ] = Index_Global_To_Local_3D(Fracture(g).N_Length_AB, Fracture(g).N_Width_AD  , 1, Ig_This);
                                    CellCenterCoord_This = Fracture(g).CellCenterCoordsV2(Ig_This,:);
                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                    % Finding if the fracture plate crosses between the matrix cell and any of its neighbors.
                                    % If it crosses, the connectivity between these two matrix cells is affected and the
                                    % transmissiblity assigned to that interface is modified by "alpha" factor of pEDFM.
                                    indNeighbors = Fracture(g).Neighbours{Ig_This};
                                    for n = 1 : length(indNeighbors)
                                        Ig_Other = indNeighbors(n);
                                        [ i_Other , j_Other , ~ ] = Index_Global_To_Local_3D(Fracture(g).N_Length_AB, Fracture(g).N_Width_AD  , 1, Ig_Other);
                                        CellCenterCoord_Other = Fracture(g).CellCenterCoordsV2(Ig_Other,:);
                                        LineSegmentConnectingCellCenters = lineSegment_DARSim( CellCenterCoord_This , CellCenterCoord_Other );
                                        LineSegmentConnectingCellCenters.AddEquation(CellCenterCoord_This);

                                        % Check if the fracture plan crosses
                                        [Geostatus, IntersectPoint] = Fracture(f).Obtain_Polygon_LineSegment_Intersection(LineSegmentConnectingCellCenters, Epsilon);

                                        if Geostatus.haveIntersect == 1 && ~isempty(IntersectPoint)
                                            [~ , Dimension] = setdiff( [i_This,j_This] , [i_Other,j_Other] );
                                            if isempty(Dimension) ,  continue;  end
                                            if length(Dimension)>1,  error('Only one of the i,j,k indices can be different.');  end

                                            % Check if the intersected cell is already on the list or not
                                            if ~ismember(Ig_Other,pEDFM_List)
                                                pEDFM_List = [pEDFM_List;Ig_Other];
                                                pEDFM_Af = [pEDFM_Af;0];
                                                pEDFM_AvgDist = [pEDFM_AvgDist;0];
                                                counter = counter+1;
                                                Index = counter;
                                            else
                                                Index = find(pEDFM_List==Ig_Other);
                                            end

                                            % Calculating the projected area fraction Af_proj:
                                            % The Cosine of the projection angle can be calculated as dot product between the
                                            % normal vectors of fracture length and projection line.
                                            switch Dimension
                                                case(1)
                                                    Aff_proj = (Aff(m,2) / 2) * abs(dot( Fracture(g).AD_vec/norm(Fracture(g).AD_vec) , IntersectLine.unitVec) );
                                                case(2)
                                                    Aff_proj = (Aff(m,2) / 2) * abs(dot( Fracture(g).AB_vec/norm(Fracture(g).AB_vec) , IntersectLine.unitVec) );
                                            end
                                            pEDFM_Af(Index) = pEDFM_Af(Index) + Aff_proj;

                                            % Calculating average distance AvgDist_proj:
                                            PointA = [ Fracture(g).GridCoords( i_Other   , j_Other+1 , 1 )
                                                       Fracture(g).GridCoords( i_Other   , j_Other+1 , 2 )
                                                       Fracture(g).GridCoords( i_Other   , j_Other+1 , 3 ) ]';
                                            PointB = [ Fracture(g).GridCoords( i_Other   , j_Other   , 1 )
                                                       Fracture(g).GridCoords( i_Other   , j_Other   , 2 )
                                                       Fracture(g).GridCoords( i_Other   , j_Other   , 3 ) ]';
                                            PointC = [ Fracture(g).GridCoords( i_Other+1 , j_Other   , 1 )
                                                       Fracture(g).GridCoords( i_Other+1 , j_Other   , 2 )
                                                       Fracture(g).GridCoords( i_Other+1 , j_Other   , 3 ) ]';
                                            PointD = [ Fracture(g).GridCoords( i_Other+1 , j_Other+1 , 1 )
                                                       Fracture(g).GridCoords( i_Other+1 , j_Other+1 , 2 )
                                                       Fracture(g).GridCoords( i_Other+1 , j_Other+1 , 3 ) ]';

                                            NeighborCell = tetragon_DARSim(PointA,PointB,PointC,PointD);
                                            Refinement = 6;
                                            AvgDist_proj = NeighborCell.Obtain_Average_Distance_Line_From_Tetragon(IntersectLine, Refinement);      
                                            pEDFM_AvgDist(Index) = AvgDist_proj;

                                            % Measuring the alpha:
                                            % alpha is a the fraction of blockage (from 0 to 1) between the
                                            % neighboring matrix cells due to existence of fracture plates when
                                            % considering pEDFM method. the alpha=0 means no blockage and
                                            % alpha=1 means total blockage.
                                            switch Dimension
                                                case(1)
                                                    alpha = (Aff_proj) / Fracture(g).D_Width_AD;
                                                    Fracture(g).alpha_Tx( max(i_This,i_Other) , j_This ) =      Fracture(g).alpha_Tx( max(i_This,i_Other) , j_This ) + alpha;
                                                    Fracture(g).alpha_Tx( max(i_This,i_Other) , j_This ) = min( Fracture(g).alpha_Tx( max(i_This,i_Other) , j_This ) , 1 );
                                                case(2)
                                                    alpha = (Aff_proj) / Fracture(g).D_Length_AB;
                                                    Fracture(g).alpha_Ty( i_This , max(j_This,j_Other) ) =      Fracture(g).alpha_Ty( i_This , max(j_This,j_Other) ) + alpha;
                                                    Fracture(g).alpha_Ty( i_This , max(j_This,j_Other) ) = min( Fracture(g).alpha_Ty( i_This , max(j_This,j_Other) ) , 1 );
                                            end
                                        end
                                    end
                                end

                                % Adding the pEDFM information to the fracture
                                pEDFM_Temp = [pEDFM_List,pEDFM_Af,pEDFM_AvgDist,pEDFM_Af./pEDFM_AvgDist];
                                if ~isempty(pEDFM_Temp)
                                    pEDFM_Temp = sortrows(pEDFM_Temp,1);
                                    Fracture(f).Af_frac_pEDFM{g}{If}      = pEDFM_Temp(:,[1,2]);
                                    Fracture(f).avgDist_frac_pEDFM{g}{If} = pEDFM_Temp(:,[1,3]);
                                    Fracture(f).CI_frac_pEDFM{g}{If}      = pEDFM_Temp(:,[1,4]);
                                    Fracture(f).NumOf_fracConn_pEDFM(If)  = Fracture(f).NumOf_fracConn_pEDFM(If) + size(pEDFM_Temp,1);
                                end
                            end
                        end
                    end
                end
            end
            
            % Each two none-neighboring cells have separate calculation for
            % connectivity index (CI). However, they need to share one 
            % unique CI. Therefore both CIs (if they co-exist) will be
            % first harmonically averaged and then will be added to the dataset.
            for f = 1 : length(Fracture)
                for If = 1 : Fracture(f).N_Total
                    for g = 1 : length(Fracture)
                        if f ~= g && ~isempty(Fracture(f).CI_frac_pEDFM{g})
                            if ~isempty(Fracture(f).CI_frac_pEDFM{g}{If})
                                Ig_List = Fracture(f).CI_frac_pEDFM{g}{If}(:,1);
                                for n = 1 : length(Ig_List)
                                    bothExist = false;
                                    Ig_This = Ig_List(n);
                                    CI_frac1 = Fracture(f).CI_frac_pEDFM{g}{If}(n,2);
                                    if ~isempty( Fracture(g).CI_frac_pEDFM{f}{Ig_This} )
                                        ind_List = Fracture(g).CI_frac_pEDFM{f}{Ig_This}(:,1);
                                        if ismember(If,ind_List)
                                            ind = find(ind_List == If);
                                            CI_frac2 = Fracture(g).CI_frac_pEDFM{f}{Ig_This}(ind,2);
                                            bothExist = true;
                                        else
                                            CI_frac2 = CI_frac1;
                                        end
                                    else
                                        CI_frac2 = CI_frac1;
                                    end
                                    CI_Harmonic = 2 * CI_frac1 * CI_frac2 / (CI_frac1 + CI_frac2);
                                    Fracture(f).CI_frac_pEDFM{g}{If}(n,2) = CI_Harmonic;
                                    if ~bothExist
                                        Fracture(g).CI_frac_pEDFM{f}{Ig_This} = [ Fracture(g).CI_frac_pEDFM{f}{Ig_This} ; [If , CI_Harmonic] ];
                                        Fracture(g).CI_frac_pEDFM{f}{Ig_This} = sortrows(Fracture(g).CI_frac_pEDFM{f}{Ig_This});
                                        Fracture(g).NumOf_fracConn_pEDFM(Ig_This) = Fracture(g).NumOf_fracConn_pEDFM(Ig_This) + 1;
                                    else
                                        Fracture(g).CI_frac_pEDFM{f}{Ig_This}(ind,2) = CI_Harmonic;
                                    end
                                end
                            end
                        end
                    end
                end
            end
            
            % Improving pEDFM alpha corrections
            obj.ReservoirGrid.pEDFM_alpha_Tx( obj.ReservoirGrid.pEDFM_alpha_Tx > 0 ) = 1;
            obj.ReservoirGrid.pEDFM_alpha_Ty( obj.ReservoirGrid.pEDFM_alpha_Ty > 0 ) = 1;
            obj.ReservoirGrid.pEDFM_alpha_Tz( obj.ReservoirGrid.pEDFM_alpha_Tz > 0 ) = 1;
            for f = 1 : length(Fracture)
                Fracture(f).alpha_Tx( Fracture(f).alpha_Tx > 0 ) = 1;
                Fracture(f).alpha_Ty( Fracture(f).alpha_Ty > 0 ) = 1;
            end
            
            fprintf('\n---------------------------------------------------------\n');
        end
    end
end