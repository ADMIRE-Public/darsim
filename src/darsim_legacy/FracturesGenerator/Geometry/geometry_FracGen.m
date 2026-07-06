% Simulation class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Mousa HosseiniMehr
%TU Delft
%Created: 12 July 2016
%Last modified: 17 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef geometry_FracGen < handle
    properties
        FractureModelType % EDFM or pEDFM
        DiscretizationType % Cartesian or CornerPointGrid
        NumOfDimensions   % 2D or 3D
        ReservoirGrid
        FracturesGrid
    end
    methods
        function Compute_EDFM_Connectivities(obj)
            %%
            warning("more import class - EDFM AND PEDFM GENERATION: take a look here")
            Fracture = obj.FracturesGrid.Fracture;
            Epsilon = obj.ReservoirGrid.Epsilon;

            %% Intializing necessary variables
            for f = 1 : length(Fracture)
                % Coordinates of intersections between each two fracture plates:
                Fracture(f).intersectCoord_fracturePlate = cell( length(Fracture) , 1 );
                % Overlap status of fracture f cells due to intersection with fracture g:
                Fracture(f).Overlap_frac = cell( length(Fracture) , 1 );
                
                Fracture(f).intersectCoord_rock = cell( Fracture(f).N_Total , 1 );   % Coordinates of intersections between each fracture cell and each matrix cell
                Fracture(f).      Af_rock_EDFM  = cell( Fracture(f).N_Total , 1 );   % Area fraction of each fracture cell inside each matrix cube
                Fracture(f). avgDist_rock_EDFM  = cell( Fracture(f).N_Total , 1 );   % Average distance between each fracture cell and each matrix cube
                Fracture(f).      CI_rock_EDFM  = cell( Fracture(f).N_Total , 1 );   % Geometrical Transmissibility between each fracture cell and each matrix cube
                Fracture(f).  Af_rock_EDFM_sum  = 0;                                 % Summation of all area fractions (only to validate)
                Fracture(f).intersectCoord_frac = cell( length(Fracture) , 1 );      % Coordinates of intersections between each two fracture cells of distinct fracture plates
                Fracture(f).       Af_frac_EDFM = cell( length(Fracture) , 1 );      % Area fraction of of intersection line between each two fracture cells of distinct fracture plates
                Fracture(f).  avgDist_frac_EDFM = cell( length(Fracture) , 1 );      % Average distance between each fracture cell and the intersection line between each two fracture cells of distinct fracture plates
                Fracture(f).       CI_frac_EDFM = cell( length(Fracture) , 1 );      % Geometrical Transmissibility between each two non-neighboring fracture cells of distinct fracture plates
                Fracture(f).NumOf_fracConn_EDFM = zeros( Fracture(f).N_Total , 1);   % Number of non-neighboring fracture cells connections
                Fracture(f).      Af_rock_pEDFM = cell(Fracture(f).N_Total,1);
                Fracture(f). avgDist_rock_pEDFM = cell(Fracture(f).N_Total,1);
                Fracture(f).      CI_rock_pEDFM = cell(Fracture(f).N_Total,1);
                Fracture(f).      Af_frac_pEDFM = cell(length(Fracture),1);
                Fracture(f). avgDist_frac_pEDFM = cell(length(Fracture),1);
                Fracture(f).      CI_frac_pEDFM = cell(length(Fracture),1);
                Fracture(f).NumOf_fracConn_pEDFM = zeros( Fracture(f).N_Total , 1);
            end
            
            %% Intersections of fracture plates with each other (if any)
            % fracture fracture intersection
            fprintf('Obtaining the intersection of fracture plates ...\n');
            for f = 1 : length(Fracture)
                for g = f+1 : length(Fracture)
                    
                    % Obtaining the intersection Points
                    [Geostatus, intersectPoints] = Fracture(f).Obtain_Tetragon_Tetragon_Intersection( Fracture(g), Epsilon );
                    intersectCoord_fracturePlate_final = [];
                    
                    if ~isempty(intersectPoints) && (Geostatus.areCoplanar ~= 1 )
                        % Processing the intersection Points
                        for nr = 1 : size(intersectPoints,1)
                            
                            % Removing the Points that are not inside the fracture rectangles
                            isInside1 = Fracture(f).Is_Point_Inside_Tetragon(intersectPoints(nr,:), Epsilon);
                            isInside2 = Fracture(g).Is_Point_Inside_Tetragon(intersectPoints(nr,:), Epsilon);
                            if ( isInside1 == 0 ) || ( isInside2 == 0 )
                                intersectPoints(nr,:) = [ NaN , NaN , NaN ];
                                continue;
                            end
                            
                            % Removing the points that are repeated
                            subtracted = ( intersectPoints - intersectPoints(nr,:) );
                            subtracted = sqrt( subtracted(:,1).^2 + subtracted(:,2).^2 + subtracted(:,3).^2 );
                            if ~isempty( find( subtracted([1:nr-1 , nr+1:end]) <= Epsilon , 1 ) )
                                intersectPoints(nr,:) = [ NaN , NaN , NaN ];
                                continue;
                            end
                            
                            if ~isnan(intersectPoints(nr,1))
                                intersectCoord_fracturePlate_final = [ intersectCoord_fracturePlate_final ; intersectPoints(nr,:) ];
                            end
                        end
                        
                        % If there is no or only one final intersection Point, do not plot
                        if size(intersectCoord_fracturePlate_final,1) < 2,  continue;  end
                        
                        % Writing the final intersection Points into the fracture data structure
                        Fracture(f).intersectCoord_fracturePlate{g} = intersectCoord_fracturePlate_final;
                        Fracture(g).intersectCoord_fracturePlate{f} = intersectCoord_fracturePlate_final;
                        Fracture(f).Overlap_frac{g} = zeros( Fracture(f).N_Length_AB * Fracture(f).N_Width_AD , 1 );
                        Fracture(g).Overlap_frac{f} = zeros( Fracture(g).N_Length_AB * Fracture(g).N_Width_AD , 1 );
                    end
                end
            end
            fprintf('---------------------------------------------------------\n');
            
            %% Obtaining fracture - matrix connectivities
            % fracture matrix intersection properties
            obj.Compute_EDFM_Fracture_Matrix_Connectivities();
            %% Obtaining fracture - fracture connectivities
            % fracture fracture intersection properties
            obj.Compute_EDFM_Fracture_Fracture_Connectivities();
        end
        function Compute_EDFM_Fracture_Fracture_Connectivities(obj)
            Fracture = obj.FracturesGrid.Fracture;
            Epsilon = obj.ReservoirGrid.Epsilon;
            %% Intersections of fracture cells of distinct fractures with each other (if any)
            fprintf('Obtaining fracture - fracture connectivities:\n');
            for f = 1 : length(Fracture)
                for g = f+1 : length(Fracture)
                    % Finding the intersections between the cells of distinct fractures
                    if ~isempty(Fracture(f).intersectCoord_fracturePlate{g})
                        fprintf('Fracture %02d <---> %02d\n',f,g);
                        Fracture(f).intersectCoord_frac{g} = cell( Fracture(f).N_Total , 1 );
                        Fracture(g).intersectCoord_frac{f} = cell( Fracture(g).N_Total , 1 );
                        
                        Fracture(f).     Af_frac_EDFM{g} = cell( Fracture(f).N_Total , 1 );
                        Fracture(f).avgDist_frac_EDFM{g} = cell( Fracture(f).N_Total , 1 );
                        Fracture(f).     CI_frac_EDFM{g} = cell( Fracture(f).N_Total , 1 );
                        Fracture(f). Af_frac_EDFM_sum{g} = 0;
                        
                        Fracture(g).     Af_frac_EDFM{f} = cell( Fracture(g).N_Total , 1 );
                        Fracture(g).avgDist_frac_EDFM{f} = cell( Fracture(g).N_Total , 1 );
                        Fracture(g).     CI_frac_EDFM{f} = cell( Fracture(g).N_Total , 1 );
                        Fracture(g). Af_frac_EDFM_sum{f} = 0;
                        
                        % Loop for fracture f
                        If_List = find(Fracture(f).Overlap_frac{g});
                        for nf = 1:length(If_List)
                            If = If_List(nf);
                            [i_f,j_f,~] = Index_Global_To_Local_3D(Fracture(f).N_Length_AB, Fracture(f).N_Width_AD, 1, If);
                            
                            % Corner Points of the 1st fracture cell
                            PointA = [ Fracture(f).GridCoords(i_f  ,j_f  ,1) , Fracture(f).GridCoords(i_f  ,j_f  ,2) , Fracture(f).GridCoords(i_f  ,j_f  ,3) ];
                            PointB = [ Fracture(f).GridCoords(i_f  ,j_f+1,1) , Fracture(f).GridCoords(i_f  ,j_f+1,2) , Fracture(f).GridCoords(i_f  ,j_f+1,3) ];
                            PointC = [ Fracture(f).GridCoords(i_f+1,j_f+1,1) , Fracture(f).GridCoords(i_f+1,j_f+1,2) , Fracture(f).GridCoords(i_f+1,j_f+1,3) ];
                            PointD = [ Fracture(f).GridCoords(i_f+1,j_f  ,1) , Fracture(f).GridCoords(i_f+1,j_f  ,2) , Fracture(f).GridCoords(i_f+1,j_f  ,3) ];
                            Plane_fracCell_1 = tetragon_DARSim(PointA,PointB,PointC,PointD);
                            
                            % Loop for fracture g
                            Ind_f_Intersect  = 0;
                            Ig_List = find(Fracture(g).Overlap_frac{f});
                            for ng = 1:length(Ig_List)
                                Ig = Ig_List(ng);
                                [i_g,j_g,~] = Index_Global_To_Local_3D(Fracture(g).N_Length_AB, Fracture(g).N_Width_AD, 1, Ig);
                                
                                % Corner Points of the 2nd fracture cell
                                PointA = [ Fracture(g).GridCoords(i_g  ,j_g  ,1) , Fracture(g).GridCoords(i_g  ,j_g  ,2) , Fracture(g).GridCoords(i_g  ,j_g  ,3) ];
                                PointB = [ Fracture(g).GridCoords(i_g  ,j_g+1,1) , Fracture(g).GridCoords(i_g  ,j_g+1,2) , Fracture(g).GridCoords(i_g  ,j_g+1,3) ];
                                PointC = [ Fracture(g).GridCoords(i_g+1,j_g+1,1) , Fracture(g).GridCoords(i_g+1,j_g+1,2) , Fracture(g).GridCoords(i_g+1,j_g+1,3) ];
                                PointD = [ Fracture(g).GridCoords(i_g+1,j_g  ,1) , Fracture(g).GridCoords(i_g+1,j_g  ,2) , Fracture(g).GridCoords(i_g+1,j_g  ,3) ];
                                Plane_fracCell_2 = tetragon_DARSim(PointA,PointB,PointC,PointD);
                                
                                [~, intersectPoints] = Plane_fracCell_1.Obtain_Tetragon_Tetragon_Intersection( Plane_fracCell_2, Epsilon );
                                
                                intersectCoord_frac_final = [];
                                
                                if ~isempty(intersectPoints)
                                    for nr = 1 : size(intersectPoints,1)
                                        
                                        % Removing the Points that are not inside the fracture tetragon
                                        isInside1 = Plane_fracCell_1.Is_Point_Inside_Tetragon(intersectPoints(nr,:), Epsilon);
                                        isInside2 = Plane_fracCell_2.Is_Point_Inside_Tetragon(intersectPoints(nr,:), Epsilon);
                                        if ( isInside1 == 0 ) || ( isInside2 == 0 )
                                            intersectPoints(nr,:) = [ NaN , NaN , NaN ];
                                            continue;
                                        end
                                        
                                        % Removing the Points that are repeated
                                        subtracted = ( intersectPoints - intersectPoints(nr,:) );
                                        subtracted = sqrt( subtracted(:,1).^2 + subtracted(:,2).^2 + subtracted(:,3).^2 );
                                        if ~isempty( find( subtracted([1:nr-1 , nr+1:end]) < Epsilon , 1 ) )
                                            intersectPoints(nr,:) = [ NaN , NaN , NaN ];
                                            continue;
                                        end
                                        
                                        if ~isnan(intersectPoints(nr,1))
                                            intersectCoord_frac_final = [ intersectCoord_frac_final ; intersectPoints(nr,:) ];
                                        end
                                    end
                                    
                                    if ( size( intersectCoord_frac_final , 1 ) > 2 )
                                        warning('More than two intersection Points exist between two intersecting fracture cells! The implementation might be wrong!');
                                    end
                                    
                                    % Obtaining the area fractions and average distances and writing them into structure
                                    if ( size( intersectCoord_frac_final , 1 ) == 2 )
                                        Ind_f_Intersect = Ind_f_Intersect + 1;
                                        Ind_g_Intersect = size(Fracture(g).intersectCoord_frac{f}{Ig},1) + 1;
                                        Fracture(f).intersectCoord_frac{g}{If}{Ind_f_Intersect,1} = Ig;
                                        Fracture(f).intersectCoord_frac{g}{If}{Ind_f_Intersect,2} = intersectCoord_frac_final;
                                        Fracture(g).intersectCoord_frac{f}{Ig}{Ind_g_Intersect,1} = If;
                                        Fracture(g).intersectCoord_frac{f}{Ig}{Ind_g_Intersect,2} = intersectCoord_frac_final;
                                        
                                        Fracture(f).NumOf_fracConn_EDFM(If) = Fracture(f).NumOf_fracConn_EDFM(If) + 1;
                                        Fracture(g).NumOf_fracConn_EDFM(Ig) = Fracture(g).NumOf_fracConn_EDFM(Ig) + 1;
                                        
                                        Fracture(f).     Af_frac_EDFM{g}{If}(Ind_f_Intersect,1) = Ig;
                                        Fracture(f).avgDist_frac_EDFM{g}{If}(Ind_f_Intersect,1) = Ig;
                                        Fracture(f).     CI_frac_EDFM{g}{If}(Ind_f_Intersect,1) = Ig;
                                        Fracture(g).     Af_frac_EDFM{f}{Ig}(Ind_g_Intersect,1) = If;
                                        Fracture(g).avgDist_frac_EDFM{f}{Ig}(Ind_g_Intersect,1) = If;
                                        Fracture(g).     CI_frac_EDFM{f}{Ig}(Ind_g_Intersect,1) = If;
                                        
                                        intersectLine = lineSegment_DARSim( intersectCoord_frac_final(1,:) , intersectCoord_frac_final(2,:) );
                                        
                                        [Af_frac1 , avgDist_frac1, Collinearity1 ] = obj.Line_Plane_Connectivity_3D( Plane_fracCell_1 , intersectLine );
                                        [Af_frac2 , avgDist_frac2, Collinearity2 ] = obj.Line_Plane_Connectivity_3D( Plane_fracCell_2 , intersectLine );

                                        if ( Collinearity1 == 1 ),  Af_frac2 = Af_frac2 /2;  end
                                        if ( Collinearity2 == 1 ),  Af_frac1 = Af_frac1 /2;  end
                                        
                                        Fracture(f).     Af_frac_EDFM{g}{If}(Ind_f_Intersect,2) = Af_frac1;
                                        Fracture(f).avgDist_frac_EDFM{g}{If}(Ind_f_Intersect,2) = avgDist_frac1;
                                        Fracture(g).     Af_frac_EDFM{f}{Ig}(Ind_g_Intersect,2) = Af_frac2;
                                        Fracture(g).avgDist_frac_EDFM{f}{Ig}(Ind_g_Intersect,2) = avgDist_frac2;
                                        
                                        % Each two none-neighboring cells have separate calculation for
                                        % connectivity index (CI). However, they need to share one unique CI.
                                        % Therefore both CIs will be first harmonically averaged and then will be added to the dataset.
                                        CI_frac1 = Af_frac1 / avgDist_frac1;
                                        CI_frac2 = Af_frac2 / avgDist_frac2;
                                        CI_Harmonic = 2 * CI_frac1 * CI_frac2 / (CI_frac1 + CI_frac2);
                                        
                                        Fracture(f).CI_frac_EDFM{g}{If}(Ind_f_Intersect,2) = CI_Harmonic;
                                        Fracture(g).CI_frac_EDFM{f}{Ig}(Ind_g_Intersect,2) = CI_Harmonic;
                                        
                                        Fracture(f).Af_frac_EDFM_sum{g} = Fracture(f).Af_frac_EDFM_sum{g} + Af_frac1;
                                        Fracture(g).Af_frac_EDFM_sum{f} = Fracture(g).Af_frac_EDFM_sum{f} + Af_frac2;
                                    end
                                end
                            end
                        end
                    end
                end
            end
            fprintf('\n---------------------------------------------------------\n');
        end
        function [ Af , AvgDist , Collinearity] = Line_Plane_Connectivity_3D( obj, Plane , Line )
            Epsilon = obj.ReservoirGrid.Epsilon;
            Collinearity = 0;
            
            %% Calculating Af of the line segment inside the plane segment
            Af = 2 * norm ( Line.AB_vec );
            AB_Inf = lineInfinite_DARSim(Plane.AB_vec);  AB_Inf.AddPoint(Plane.PointA);  AB_Inf.AddEquation(Plane.PointA);
            BC_Inf = lineInfinite_DARSim(Plane.BC_vec);  BC_Inf.AddPoint(Plane.PointB);  BC_Inf.AddEquation(Plane.PointA);
            CD_Inf = lineInfinite_DARSim(Plane.CD_vec);  CD_Inf.AddPoint(Plane.PointC);  CD_Inf.AddEquation(Plane.PointA);
            DA_Inf = lineInfinite_DARSim(Plane.DA_vec);  DA_Inf.AddPoint(Plane.PointD);  DA_Inf.AddEquation(Plane.PointA);
            
            [Geostatus_AB, ~] = Line.Obtain_LineInfinite_LineInfinite_Intersection( AB_Inf, Epsilon );
            [Geostatus_BC, ~] = Line.Obtain_LineInfinite_LineInfinite_Intersection( BC_Inf, Epsilon );
            [Geostatus_CD, ~] = Line.Obtain_LineInfinite_LineInfinite_Intersection( CD_Inf, Epsilon );
            [Geostatus_DA, ~] = Line.Obtain_LineInfinite_LineInfinite_Intersection( DA_Inf, Epsilon );
            
            if ( Geostatus_AB.areCollinear == 1 ) || ( Geostatus_BC.areCollinear == 1 ) || ...
               ( Geostatus_CD.areCollinear == 1 ) || ( Geostatus_DA.areCollinear == 1 )
                Collinearity = 1;
                Af = Af /2;
            end
            
            %% Obtaining aveDist of the plane segment from the line segment
            AvgDist  = Plane.Obtain_Average_Distance_Line_From_Tetragon(Line, 6);
        end
        function PrintInfo(obj)
            disp( ['************* ', num2str(obj.NumOfDimensions), 'D-', obj.FractureModelType, ' FractureGenerator for ', obj.DiscretizationType, ' Discretization **************']);
        end
        function PlotFractures(obj)
            %% Plotting fracture plates
            Fracture = obj.FracturesGrid.Fracture;
            fprintf('Plotting fractures ...\n');
%             % Apenas um teste: --------------------------------------------
%             minK = 1; maxK = 0;
%             for f = 1 : length(Fracture)
%                 if Fracture(f).Permeability<minK
%                    minK = Fracture(f).Permeability; 
%                 end
%                 if Fracture(f).Permeability>maxK
%                    maxK = Fracture(f).Permeability; 
%                 end
%             end
%             %--------------------------------------------------------------
            figure();
            for f = 1 : length(Fracture)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Plotting the boundaries of each fracture plate
                
%                 % Apenas um teste: ----------------------------------------
                 % if Fracture(f).Permeability==maxK, color = 'r'; elseif Fracture(f).Permeability==minK, color = 'b';end
%                 %----------------------------------------------------------
                color = 'r'
                plot3(Fracture(f).Vertex([1:4,1],1),Fracture(f).Vertex([1:4,1],2),Fracture(f).Vertex([1:4,1],3),color,'LineWidth',3);
                xlim([0,obj.ReservoirGrid.Lx]); ylim([0,obj.ReservoirGrid.Ly]); zlim([0,obj.ReservoirGrid.Lz]);
                xlabel('X [cm]'); ylabel('Y [cm]'); zlabel('Z [cm]');
                %grid ON;
                
                % Plotting the normal vector of each fracture plate
                hold on;
                quiver3(Fracture(f).Centroid(1),Fracture(f).Centroid(2),Fracture(f).Centroid(3),Fracture(f).nVec_plot(1),Fracture(f).nVec_plot(2),Fracture(f).nVec_plot(3),'b','LineWidth',2);
                
                % Plotting cell edges of each fracture plate
                for i = 1 : Fracture(f).N_Length_AB+1
                    plot3(Fracture(f).GridCoords(i,:,1),Fracture(f).GridCoords(i,:,2),Fracture(f).GridCoords(i,:,3),color,'LineWidth',1);
                end
                for j = 1 : Fracture(f).N_Width_AD+1
                    plot3(Fracture(f).GridCoords(:,j,1),Fracture(f).GridCoords(:,j,2),Fracture(f).GridCoords(:,j,3),color,'LineWidth',1);
                end
                
                % Writing the letters of A,B,C,D to the corners
                text( Fracture(f).PointA(1) - 0.05*Fracture(f).AD_vec(1) - 0.05*Fracture(f).AB_vec(1) , ...
                      Fracture(f).PointA(2) - 0.05*Fracture(f).AD_vec(2) - 0.05*Fracture(f).AB_vec(2) , ...
                      Fracture(f).PointA(3) - 0.05*Fracture(f).AD_vec(3) - 0.05*Fracture(f).AB_vec(3) , ...
                      'A' , 'Color' , 'black' , 'FontSize' , 14);
                text( Fracture(f).PointB(1) - 0.05*Fracture(f).AD_vec(1) + 0.05*Fracture(f).AB_vec(1) , ...
                      Fracture(f).PointB(2) - 0.05*Fracture(f).AD_vec(2) + 0.05*Fracture(f).AB_vec(2) , ...
                      Fracture(f).PointB(3) - 0.05*Fracture(f).AD_vec(3) + 0.05*Fracture(f).AB_vec(3) , ...
                      'B' , 'Color' , 'black' , 'FontSize' , 14);
                text( Fracture(f).PointC(1) + 0.05*Fracture(f).AD_vec(1) + 0.05*Fracture(f).AB_vec(1) , ...
                      Fracture(f).PointC(2) + 0.05*Fracture(f).AD_vec(2) + 0.05*Fracture(f).AB_vec(2) , ...
                      Fracture(f).PointC(3) + 0.05*Fracture(f).AD_vec(3) + 0.05*Fracture(f).AB_vec(3) , ...
                      'C' , 'Color' , 'black' , 'FontSize' , 14);
                text( Fracture(f).PointD(1) + 0.05*Fracture(f).AD_vec(1) - 0.05*Fracture(f).AB_vec(1) , ...
                      Fracture(f).PointD(2) + 0.05*Fracture(f).AD_vec(2) - 0.05*Fracture(f).AB_vec(2) , ...
                      Fracture(f).PointD(3) + 0.05*Fracture(f).AD_vec(3) - 0.05*Fracture(f).AB_vec(3) , ...
                      'D' , 'Color' , 'black' , 'FontSize' , 14);
                text( Fracture(f).Centroid(1), Fracture(f).Centroid(2), Fracture(f).Centroid(3), ...
                      num2str(f) , 'Color' , 'magenta' , 'FontSize' , 20, 'FontWeight','bold');
                
                % Plotting the intersection line of each two fracture plates
                for g = f+1 : length(Fracture)
                    if ~isempty(Fracture(f).intersectCoord_fracturePlate{g})
                        plot3(Fracture(f).intersectCoord_fracturePlate{g}(:,1),...
                              Fracture(f).intersectCoord_fracturePlate{g}(:,2),...
                              Fracture(f).intersectCoord_fracturePlate{g}(:,3),...
                              'Color',[0,0.5,0],'LineWidth',5);
                    end
                end
            end
            
            rotate3d on;
            view([0 90]);
            fprintf('---------------------------------------------------------\n');
        end
    end
end