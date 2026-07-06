% class fracturer for DARSim2FracGen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Mousa HosseiniMehr
%TU Delft
%Created: 2017-07-12
%Last modified: 2019-02-08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef fracture_FracGen < tetragon_DARSim
    properties
        Aperture
        Porosity
        Permeability
        Area_Total
        nVec_plot
        Width_AD
        N_Length_AB
        N_Width_AD
        N_Total
        ADM
        CoarseNodeIndex
        D_Length_AB
        D_Width_AD
        D_Length_AB_vec
        D_Width_AD_vec
        CellCenterCoordsV1
        CellCenterCoordsV2
        GridCoords
        intersectCoord_rock
        Af_rock_EDFM
        Af_rock_EDFM_sum
        avgDist_rock_EDFM
        Coplanar_EDFM
        CI_rock_EDFM
        Af_rock_pEDFM
        avgDist_rock_pEDFM
        CI_rock_pEDFM
        MatrixOverlapList = [];
        intersectCoord_fracturePlate
        intersectCoord_frac
        Overlap_frac
        Af_frac_EDFM
        Af_frac_EDFM_sum
        avgDist_frac_EDFM
        CI_frac_EDFM
        NumOf_fracConn_EDFM
        Af_frac_pEDFM
        avgDist_frac_pEDFM
        CI_frac_pEDFM
        NumOf_fracConn_pEDFM
        Neighbours
        alpha_Tx
        alpha_Ty
        FitsInsideReservoir
        % alpha is a the fraction of blockage (from 0 to 1) between
        % neighboring fracture cells due o existence of other intersecting
        % fracture planes when considering pEDFM method. The alpha=0 means 
        % no blockage and alpha=1 means total blockage.
    end
    methods
        
        function obj = fracture_FracGen(PointA,PointB,PointC,PointD)
            if nargin==4
                obj.Initializetetragon(obj,PointA,PointB,PointC,PointD);
            end
        end
        function BuildFracture(obj, Frac_Input, ReservoirGrid, Epsilon, f)
            % Dimensions of Reservoir in each direction
            Lx = ReservoirGrid.Lx;
            Ly = ReservoirGrid.Ly;
            Lz = ReservoirGrid.Lz;
            
            if isnan( Frac_Input.CenterCoord(1) )
                % Fracture construction with three Points input type:
                % The three vectors connecting 3 corners of the fracture plate (forming a triangle)
                Vec12 = Frac_Input.CornerCoords(:,2) - Frac_Input.CornerCoords(:,1);
                Vec23 = Frac_Input.CornerCoords(:,3) - Frac_Input.CornerCoords(:,2);
                Vec31 = Frac_Input.CornerCoords(:,3) - Frac_Input.CornerCoords(:,1);
                
                % Finding which of the three Points is the corner with the 90 deg angle
                dotProduct = [ dot( Vec31 , Vec12 ) , dot( Vec12 , Vec23 ) , dot( Vec23 , Vec31 ) ];
                PointRightAngle   = find( abs(dotProduct) <= Epsilon , 1 );
                PointSalientAngle = find( abs(dotProduct) >  Epsilon     );
                if isempty(PointRightAngle),  error('Fracture # %02d is not a rectangle! Please check the input file.',f);  end
                
                % Setting attributes to the corners (A is the corner with 90 deg angle)
                PointA = Frac_Input.CornerCoords( : , PointRightAngle      )';
                PointB = Frac_Input.CornerCoords( : , PointSalientAngle(1) )';
                PointD = Frac_Input.CornerCoords( : , PointSalientAngle(2) )';
                
                % The central Point of the plate
                PointM = ( PointB + PointD ) /2;
                
                % Obtaining the 4th corner (by using the centeral Point)
                PointC = PointA + 2* ( PointM - PointA );
            else
                
                % Fracture construction with central Point and angles input type
                PointM = Frac_Input.CenterCoord;  PointM = reshape(PointM,1, []);
                
                % Obtaining the corners before rotation
                PointA_noRot       = [  PointM(1) - Frac_Input.Length /2  ,  PointM(2) + Frac_Input.Width /2  ,  PointM(3)  ];
                PointB_noRot       = [  PointM(1) + Frac_Input.Length /2  ,  PointM(2) + Frac_Input.Width /2  ,  PointM(3)  ];
                PointC_noRot       = [  PointM(1) + Frac_Input.Length /2  ,  PointM(2) - Frac_Input.Width /2  ,  PointM(3)  ];
                PointD_noRot       = [  PointM(1) - Frac_Input.Length /2  ,  PointM(2) - Frac_Input.Width /2  ,  PointM(3)  ];
                
                % Rotating around Z axis with horizontal rotation angle
                RotAxis_Hor        = [ 0 ; 0 ; 1 ];
                RotAngle_Hor_Rad   = Frac_Input.RotAngleAlongZ * pi / 180;
                PointA_Rot_Hor     = Rotate_Point_Around_Line_3D( PointA_noRot' , RotAxis_Hor , PointM' , RotAngle_Hor_Rad )';
                PointB_Rot_Hor     = Rotate_Point_Around_Line_3D( PointB_noRot' , RotAxis_Hor , PointM' , RotAngle_Hor_Rad )';
                PointC_Rot_Hor     = Rotate_Point_Around_Line_3D( PointC_noRot' , RotAxis_Hor , PointM' , RotAngle_Hor_Rad )';
                PointD_Rot_Hor     = Rotate_Point_Around_Line_3D( PointD_noRot' , RotAxis_Hor , PointM' , RotAngle_Hor_Rad )';
                
                % Rotating around "BMC - DMA" (along central length)
                PointBMC           = ( PointB_Rot_Hor + PointC_Rot_Hor ) /2;
                PointDMA           = ( PointD_Rot_Hor + PointA_Rot_Hor ) /2;
                RotAxis_Length     = PointDMA - PointBMC;
                RotAxis_Length_Rad  = Frac_Input.RotAngleAlongL * pi / 180;
                PointA_Rot_Length   = Rotate_Point_Around_Line_3D( PointA_Rot_Hor' , RotAxis_Length' , PointM' , RotAxis_Length_Rad )';
                PointB_Rot_Length   = Rotate_Point_Around_Line_3D( PointB_Rot_Hor' , RotAxis_Length' , PointM' , RotAxis_Length_Rad )';
                PointC_Rot_Length   = Rotate_Point_Around_Line_3D( PointC_Rot_Hor' , RotAxis_Length' , PointM' , RotAxis_Length_Rad )';
                PointD_Rot_Length   = Rotate_Point_Around_Line_3D( PointD_Rot_Hor' , RotAxis_Length' , PointM' , RotAxis_Length_Rad )';
                
                % Rotating around "AMB - CMD" (along central width)
                PointAMB           = ( PointA_Rot_Length + PointB_Rot_Length ) /2;
                PointCMD           = ( PointC_Rot_Length + PointD_Rot_Length ) /2;
                RotAxis_Width      = PointCMD - PointAMB;
                RotAxis_Width_Rad = Frac_Input.RotAngleAlongW * pi / 180;
                PointA     = Rotate_Point_Around_Line_3D( PointA_Rot_Length' , RotAxis_Width' , PointM' , RotAxis_Width_Rad )';
                PointB     = Rotate_Point_Around_Line_3D( PointB_Rot_Length' , RotAxis_Width' , PointM' , RotAxis_Width_Rad )';
                PointC     = Rotate_Point_Around_Line_3D( PointC_Rot_Length' , RotAxis_Width' , PointM' , RotAxis_Width_Rad )';
                PointD     = Rotate_Point_Around_Line_3D( PointD_Rot_Length' , RotAxis_Width' , PointM' , RotAxis_Width_Rad )';
            end
            
            %% Starting with discretization
            % Initializing the Plane Segment properties of the fracture
            obj.InitializeTetragon(PointA,PointB,PointC,PointD);
            
            % Checking if fracture fits inside the matrix box
            if ( min(obj.Vertex(:,1)) < 0 - Epsilon ) || ( max(obj.Vertex(:,1)) > Lx + Epsilon ) || ...
               ( min(obj.Vertex(:,2)) < 0 - Epsilon ) || ( max(obj.Vertex(:,2)) > Ly + Epsilon ) || ...
               ( min(obj.Vertex(:,3)) < 0 - Epsilon ) || ( max(obj.Vertex(:,3)) > Lz + Epsilon )
               obj.FitsInsideReservoir = 0;
            else
               obj.FitsInsideReservoir = 1;
            end

            % The area of the fracture plane (each plane has two sides)
            obj.Area_Total = Triangle_Area_3D( obj.PointA , obj.PointB , obj.PointC ) * 2 * 2;
            
            % The normal vector of plate with its size normalized to 1/10 of the reservoir length (Lx) (just for plotting)
            obj.nVec_plot = obj.nVec * Lz / 10;
            
            % Check if the fracture plate is paralel to any of X,Y,Z axes
            obj.CheckIfParallelToReferenceSurface(Epsilon);
            
            % Length (Length_AB) and width (Width_AD) of the fracture plate
            if isnan( Frac_Input.CenterCoord(1) )
                obj.Length_AB = norm( obj.AB_vec );
                obj.Width_AD  = norm( obj.AD_vec );
            else
                obj.Length_AB = Frac_Input.Length;
                obj.Width_AD  = Frac_Input.Width;
                if ReservoirGrid.Nz == 1 && ~obj.FitsInsideReservoir
                    if obj.Width_AD ~= ReservoirGrid.Lz
                        warning('The length of the fracture # %02d was not modified. Make sure that you have set the length properly in the input file', f);
                        warning('The width of the fracture # %02d did not match with the thickness of the 2D reservoir and is now modified!', f);
                        obj.Width_AD  = ReservoirGrid.Lz;
                    else
                        % All is fine. No action should be taken.
                    end
                end
                if ReservoirGrid.Ny == 1 && ~obj.FitsInsideReservoir
                    if obj.Width_AD ~= ReservoirGrid.Ly
                        warning('The length of the fracture # %02d was not modified. Make sure that you have set the length properly in the input file', f);
                        warning('The width of the fracture # %02d did not match with the thickness of the 2D reservoir and is now modified!', f);
                        obj.Width_AD  = ReservoirGrid.Ly;
                    else
                        % All is fine. No action should be taken.
                    end
                end
                if ReservoirGrid.Nx == 1 && ~obj.FitsInsideReservoir
                    if obj.Width_AD ~= ReservoirGrid.Lx
                        warning('The length of the fracture # %02d was not modified. Make sure that you have set the length properly in the input file', f);
                        warning('The width of the fracture # %02d did not match with the thickness of the 2D reservoir and is now modified!', f);
                        obj.Width_AD  = ReservoirGrid.Lx;
                    else
                        % All is fine. No action should be taken.
                    end
                end
            end
            
            % Adding ADM configuration
            obj.ADM = Frac_Input.ADM;
            
            % Discretizing the fracture plane
            obj.Discretize(Frac_Input, ReservoirGrid, f);
            
            % Aperture of the fracture
            obj.Aperture     = Frac_Input.Aperture;
            if isnan(obj.Aperture)
                obj.Aperture = mean(ReservoirGrid.CellVolume)/10;
            end
            
            % Porosity and Permeability of the fracture
            if ~isempty(ReservoirGrid.Stress)
                Frac_Input = obj.ModifyFracturesPropertiesBasedOnStress(ReservoirGrid.Stress, Frac_Input);
            end
            obj.Porosity     = Frac_Input.Porosity;
            obj.Permeability = Frac_Input.Permeability;
            
            
            
            %% Assigning Coordinates (x,y,z) of cell centers and cell interfaces of the fracture plane
            
            obj.AssignCellCoordinates();
        end
        function Discretize(obj, Frac_Input, ReservoirGrid, f)
            %%
            % Number of reservoir grid cells in each direction
            Nx = ReservoirGrid.Nx;
            Ny = ReservoirGrid.Ny;
            Nz = ReservoirGrid.Nz;
            dV = mean( ReservoirGrid.CellVolume );
            
            % Obtaining the number of grid cells at each direction
            if ~isnan( Frac_Input.GridResRatio )
                % Discretizing with "GridResRatio"
                obj.N_Length_AB = round( obj.Length_AB / dV^(1/3) * Frac_Input.GridResRatio );
                obj.N_Width_AD  = round( obj.Width_AD / dV^(1/3) * Frac_Input.GridResRatio );
                if ReservoirGrid.ADM.Active && obj.ADM(2) > 0
                    % There are two approaches for coarse grid constuctio:
                    % 1- Original, with identical primal coarse grids. In
                    % this case, fine grid number must be odd.
                    % 2- Coarse nodes on the corners and edges, resulting
                    % in unequal primal coarse grids (coarse grids on the
                    % corners and edges have smaller size as they are cut).
                    % In this case,  fine grid number must be even.
                    if mod(Nx,ReservoirGrid.ADM.CoarseRatio(1)^ReservoirGrid.ADM.CoarseLevel) == 0 && ...
                       mod(Ny,ReservoirGrid.ADM.CoarseRatio(2)^ReservoirGrid.ADM.CoarseLevel) == 0
                        % Coarse grids will be constructed with original method
                        % Correcting N_Length_AB to match MMs level and coarsening ratio
                        if ( mod( obj.N_Length_AB , obj.ADM(3)^obj.ADM(2) )~=0 )
                            obj.N_Length_AB = obj.N_Length_AB - mod( obj.N_Length_AB , obj.ADM(3)^obj.ADM(2) );
                        end
                        % Correcting N_Width_AD to match MMs level and coarsening ratio
                        if ( mod( obj.N_Width_AD  , obj.ADM(4)^obj.ADM(2) )~=0 )
                            obj.N_Width_AD  = obj.N_Width_AD  - mod( obj.N_Width_AD  , obj.ADM(4)^obj.ADM(2) );
                        end
                        % At least 2 coarse grid cells at the highest coarsening level in each direction
                        obj.N_Length_AB = max( obj.N_Length_AB , 2*obj.ADM(3)^obj.ADM(2) );
                        obj.N_Width_AD  = max( obj.N_Width_AD  , 2*obj.ADM(4)^obj.ADM(2) );
                    elseif mod(Nx-1,ReservoirGrid.ADM.CoarseRatio(1)^ReservoirGrid.ADM.CoarseLevel) == 0 && ...
                           mod(Ny-1,ReservoirGrid.ADM.CoarseRatio(2)^ReservoirGrid.ADM.CoarseLevel) == 0
                        % Coarse nodes will be put on the corners as well.
                        % Correcting N_Length_AB to match MMs level and coarsening ratio
                        if ( mod( obj.N_Length_AB-1 , obj.ADM(3)^obj.ADM(2) )~=0 )
                            obj.N_Length_AB = obj.N_Length_AB - mod( obj.N_Length_AB-1 , obj.ADM(3)^obj.ADM(2) );
                        end
                        % Correcting N_Width_AD to match MMs level and coarsening ratio
                        if ( mod( obj.N_Width_AD-1  , obj.ADM(4)^obj.ADM(2) )~=0 )
                            obj.N_Width_AD  = obj.N_Width_AD  - mod( obj.N_Width_AD  , obj.ADM(4)^obj.ADM(2) );
                        end
                        % At least 2 coarse grid cells at the highest coarsening level in each direction
                        obj.N_Length_AB = max( obj.N_Length_AB , obj.ADM(3)^obj.ADM(2)+1 );
                        obj.N_Width_AD  = max( obj.N_Width_AD  , obj.ADM(4)^obj.ADM(2)+1 );
                    else
                        error('The number of grid cells in reservoir do not match its ADM settings. Please check the input file!');
                    end
                end

                % At least 2 grid cells in each direction
                % obj.N_Length_AB = max( obj.N_Length_AB , 2);
                % obj.N_Width_AD  = max( obj.N_Width_AD  , 2);
                
                % If reservoir has only 1 grid cell in Z direction 
                % (i.e., 2D reservoir), fracture must also have only
                % 1 grid cell along its width (i.e., 1D fracture)
                if Nz == 1,  obj.N_Width_AD = 1;  end
                
                % If reservoir has only 1 grid cell in Y direction 
                % (i.e., 2D reservoir), fracture must also have only
                % 1 grid cell along its width (i.e., 1D fracture)
                if Ny == 1,  obj.N_Width_AD = 1;  end
                
                % If reservoir has only 1 grid cell in X direction 
                % (i.e., 2D reservoir), fracture must also have only
                % 1 grid cell along its width (i.e., 1D fracture)
                if Nx == 1,  obj.N_Width_AD = 1;  end
                
            else
                
                % Discretizing with "GridNumAlongL" and "GridNumAlongW"
                obj.N_Length_AB = Frac_Input.GridNumAlongL;
                obj.N_Width_AD  = Frac_Input.GridNumAlongW;
                
                if ReservoirGrid.ADM.Active && obj.ADM(2) > 0
                    if mod(Nx,2)==1 && mod(Ny,2)==1
                        % In case original coarse grid construction, with identical primal coarse grids.
                        % In this case, fine grid number must be odd.
                        if mod( obj.N_Length_AB -1 , obj.ADM(3)^obj.ADM(2) ) ~= 0
                            error('In fracture #%02d, the number of grid cells along length is not acceptable to the ADM properties!\nPlease check the input file.',f);
                        end
                        if mod( obj.N_Width_AD-1 , obj.ADM(4)^obj.ADM(2) ) ~= 0
                            error('In fracture #%02d, the number of grid cells along width is not acceptable to the ADM properties!\nPlease check the input file.',f);
                        end
                    else
                        % In case of coarse nodes on the corners and edges, resulting
                        % in unequal primal coarse grids (coarse grids on the
                        % corners and edges have smaller size as they are cut).
                        % In this case, fine grid number must be even.
                        if mod( obj.N_Length_AB , obj.ADM(3)^obj.ADM(2) ) ~= 0
                            error('In fracture #%02d, the number of grid cells along length is not acceptable to the ADM properties!\nPlease check the input file.',f);
                        end
                        if mod( obj.N_Width_AD  , obj.ADM(4)^obj.ADM(2) ) ~= 0
                            error('In fracture #%02d, the number of grid cells along width is not acceptable to the ADM properties!\nPlease check the input file.',f);
                        end
                    end
                    
                end
            end
            
            % Total number of grid cells in the fracture
            obj.N_Total = obj.N_Length_AB * obj.N_Width_AD;
            
            % The size of each fracture grid cell in AB and AD directions with their vectors
            obj.D_Length_AB     = obj.Length_AB / obj.N_Length_AB;
            obj.D_Width_AD      = obj.Width_AD  / obj.N_Width_AD;
            obj.D_Length_AB_vec = obj.AB_vec    / obj.N_Length_AB;
            obj.D_Width_AD_vec  = obj.AD_vec    / obj.N_Width_AD;
            
            % Add the neighbors of each fracture cell
            obj.AssignNeighbours();
            
            % alpha is a the fraction of blockage (from 0 to 1) between
            % neighboring fracture cells due to existence of other intersecting
            % fracture planes when considering pEDFM method. The alpha=0 means
            % no blockage and alpha=1 means total blockage.
            obj.alpha_Tx = zeros( obj.N_Length_AB+1, obj.N_Width_AD   );
            obj.alpha_Ty = zeros( obj.N_Length_AB  , obj.N_Width_AD+1 );
        end
        function AssignCellCoordinates(obj)
            %% Coordinates (x,y,z) of cell centers and cell interfaces of the fracture plate
            obj.CellCenterCoordsV1  = zeros( obj.N_Length_AB   , obj.N_Width_AD   , 3 );
            obj.CellCenterCoordsV2  = zeros( obj.N_Length_AB * obj.N_Width_AD , 3 );
            obj.GridCoords = zeros( obj.N_Length_AB+1 , obj.N_Width_AD+1 , 3 );
            
            %% Coordinates of fracture cell centers
            for i = 1 : obj.N_Length_AB
                obj.CellCenterCoordsV1(i,:,1) = linspace( obj.PointA(1) + (i-1/2)*obj.D_Length_AB_vec(1) + obj.D_Width_AD_vec(1)/2 , ...
                                                          obj.PointD(1) + (i-1/2)*obj.D_Length_AB_vec(1) - obj.D_Width_AD_vec(1)/2 , obj.N_Width_AD );
                
                obj.CellCenterCoordsV1(i,:,2) = linspace( obj.PointA(2) + (i-1/2)*obj.D_Length_AB_vec(2) + obj.D_Width_AD_vec(2)/2 , ...
                                                          obj.PointD(2) + (i-1/2)*obj.D_Length_AB_vec(2) - obj.D_Width_AD_vec(2)/2 , obj.N_Width_AD );
                
                obj.CellCenterCoordsV1(i,:,3) = linspace( obj.PointA(3) + (i-1/2)*obj.D_Length_AB_vec(3) + obj.D_Width_AD_vec(3)/2 , ...
                                                          obj.PointD(3) + (i-1/2)*obj.D_Length_AB_vec(3) - obj.D_Width_AD_vec(3)/2 , obj.N_Width_AD );
            end
            
            %% Coordinates of fracture cell centers in vector form
            for j = 1 : obj.N_Width_AD
                obj.CellCenterCoordsV2( (j-1)*obj.N_Length_AB+1 : j*obj.N_Length_AB , : ) = obj.CellCenterCoordsV1(:,j,:);
            end
            
            %% Coordinates of fracture cell interfaces
            for i = 1 : obj.N_Length_AB+1
                obj.GridCoords(i,:,1) = linspace( obj.PointA(1) + (i-1)*obj.D_Length_AB_vec(1) , ...
                                                  obj.PointD(1) + (i-1)*obj.D_Length_AB_vec(1) , obj.N_Width_AD+1 );
                
                obj.GridCoords(i,:,2) = linspace( obj.PointA(2) + (i-1)*obj.D_Length_AB_vec(2) , ...
                                                  obj.PointD(2) + (i-1)*obj.D_Length_AB_vec(2) , obj.N_Width_AD+1 );
                
                obj.GridCoords(i,:,3) = linspace( obj.PointA(3) + (i-1)*obj.D_Length_AB_vec(3) , ...
                                                  obj.PointD(3) + (i-1)*obj.D_Length_AB_vec(3) , obj.N_Width_AD+1 );
            end
        end
        function AssignNeighbours(obj)
            Nx = obj.N_Length_AB;
            Ny = obj.N_Width_AD;
            Nz = 1;
            N = obj.N_Total;
            % West Neighbors
            [i,j,k] = meshgrid( 0:Nx-1 , 1:Ny , 1:Nz );  i(i<1) = nan;
            IW = (k-1).*Nx.*Ny + (j-1).*Nx + i;
            IW = reshape(permute(IW,[2 1 3]),N,1);
            
            % East Neighbors
            [i,j,k] = meshgrid( 2:Nx+1 , 1:Ny , 1:Nz );  i(i>Nx) = nan;
            IE = (k-1).*Nx.*Ny + (j-1).*Nx + i;
            IE = reshape(permute(IE,[2 1 3]),N,1);
            
            % South Neighbors
            [i,j,k] = meshgrid( 1:Nx , 0:Ny-1 , 1:Nz );  j(j<1) = nan;
            IS = (k-1).*Nx.*Ny + (j-1).*Nx + i;
            IS = reshape(permute(IS,[2 1 3]),N,1);
            
            % North Neighbors
            [i,j,k] = meshgrid( 1:Nx , 2:Ny+1 , 1:Nz );  j(j>Ny) = nan;
            IN = (k-1).*Nx.*Ny + (j-1).*Nx + i;
            IN = reshape(permute(IN,[2 1 3]),N,1);
            
            % Bottom Neighbors
            [i,j,k] = meshgrid( 1:Nx , 1:Ny , 0:Nz-1 );  k(k<1) = nan;
            IB = (k-1).*Nx.*Ny + (j-1).*Nx + i;
            IB = reshape(permute(IB,[2 1 3]),N,1);
            
            % Top Neighbors
            [i,j,k] = meshgrid( 1:Nx , 1:Ny , 2:Nz+1 );  k(k>Nz) = nan;
            IT = (k-1).*Nx.*Ny + (j-1).*Nx + i;
            IT = reshape(permute(IT,[2 1 3]),N,1);
            
            % Assembling the array for all the neighbors
            obj.Neighbours = [IW,IE,IS,IN,IB,IT];
            obj.Neighbours = num2cell(obj.Neighbours,2);
            obj.Neighbours = cellfun(@(x) x(~isnan(x)),obj.Neighbours ,'UniformOutput' ,false);
        end
        function Frac_Input = ModifyFracturesPropertiesBasedOnStress(obj,Stress,Frac_Input)
            %% This function Modifies the properties of the fractures based on the stress field of the reservoir.
            % Written by MSc student Daan Ris.
            s1_e = Stress.Sigma1 - Stress.Pressure;
            s3_e = Stress.Sigma3 - Stress.Pressure;
            
            %% ------------ Obtaining mohr circle radius and centres ------------
            Diam1 = s1_e-s3_e;                   %Diameter of the Sigma 1, Sigma 3 circle (1 refers to circle #1)
            r1 = 0.5*Diam1;                      %Radius of the Sigma 1, Sigma 3 circle
            Centre1 = r1 + s3_e ;
            a = Centre1;
            
            %% ------------ Plotting Mohr-Coulomb circle ------------
            x = s3_e:0.1:s1_e;                 %x values for circle ranging from 0 to sigma 1 with steps of 0.1
            y1 = (r1.^2-(x-a).^2).^0.5;     %y values for circle for sigma 3 and sigma 1
            
            figure; plot (x,y1)             %Plotting mohr circle for sigma 1 and 3
            hold on
            ylim([0 s1_e])
            ylim([0 s1_e])
            
            %% ------------ Plotting Friction coefficient of 0.6 ------------
            f_c = 0.6;                    %Plot the friction coefficient line mu=1
            x2 = 0:0.1:s1_e;
            y2 = f_c * x2;
            plot(x2,y2)
            
            %% ------------ Make abc formula to find intersection friction coefficient and mohr circle ------------
            
            a_abc = 1+(f_c)^2;                %a value for abc formula
            b_abc = -2*a;                       %b value for abc formula
            c_abc = -(r1^2-a^2);                %c value for abc formula
            D_abc = b_abc^2 - (4*a_abc*c_abc);      %D value for abc formula
            if D_abc < 0
                message = sprintf('no intersection');
                return
            end
            
            x_i1 = (-b_abc + D_abc^0.5)/(2*a_abc);  %intersection point 1 (right) (x value)
            x_i2 = (-b_abc - D_abc^0.5)/(2*a_abc);  %intersection point 2 (left) (x value)
            
            y_i1 = x_i1*f_c;                  %intersection point 1 (right)(y value)
            y_i2 = x_i2*f_c;                  %intersection point 2 (left) (y value)
            
            plot(x_i1,y_i1,'d')                 %right intersection
            plot(x_i2,y_i2,'d')                 %left intersection
            
            %% ------------ Plotting lines to centre of mohr circle left intersection ------------
            xslope = x_i2:0.01:Centre1;              %x values for the slope of the line connecting intersection with the centre of the circle
            slope1 = y_i2/(x_i2-Centre1);           %slope of the linear line connecting intersection with the centre of the circle
            y_slope = (slope1*(xslope-x_i2))+y_i2;  %y values for the slope of the line connecting intersection with the centre of the circle
            plot(xslope,y_slope)
            
            %% ------------ Plotting lines to centre of mohr circle right intersection ------------
            if x_i1 < Centre1
                slope2 = y_i1/(x_i1-Centre1);
                xslope2 = x_i1:0.01:Centre1+1;
                y_slope2 = (slope2*(xslope2-x_i1))+y_i1;
            elseif x_i1 > Centre1
                slope2 = y_i1/(x_i1-Centre1);
                xslope2 = Centre1:0.05:x_i1+0.04;
                y_slope2 = (slope2*(xslope2-x_i1))+y_i1;
                
            end
            
            plot(xslope2,y_slope2);
            hold off
            %% ------------ Calculate angle between the two linear lines ------------
            beta1 = atan(slope1)* 180/pi;
            if beta1>0
                beta1 = beta1;
            elseif beta1<0
                beta1 = beta1+180;
            end
            beta2 = atan(slope2)* 180/pi;
            theta1 = .5*beta1;
            theta2 = .5*beta2;
            
            %% ------------ Calculate ranges for crittically stressed fractures ------------
            %Because theta1 and theta2 represent the poles of the fractures, this must
            %be converted to the orientation of the fractures. This is done by adding
            %90 degrees to theta1 and theta2. Also move according to direction of S_1
            gamma1 = 90-theta1;
            gamma2 = 90-theta2;
            gamma3 = 180 - gamma2;
            gamma4 = 180 - gamma1;
            
            %% ------------ Check if the fracture orientation is within the range ------------
            fracture_orientation = Frac_Input.RotAngleAlongZ;
            
            if gamma1<fracture_orientation && fracture_orientation<gamma2
                frac = 1;
            elseif gamma3<fracture_orientation && fracture_orientation<gamma4
                frac = 1;
            else
                frac = 0;
            end
            if frac == 1
                message = sprintf('fracture is hydraulically active');
                Frac_Input.Permeability = 1.0e-06;
            elseif frac == 0
                message = sprintf('fracture is not hydraulically active');
                Frac_Input.Permeability = 1.0e-22;
            end
        end
        function PrintInfo(obj,f)
            fprintf('Frac %2d: Size= %4.2f x %4.2f [m2], Grid= %2.0f x %2.0f, ADM Switch=%1.0f, Coarse Level=%1.0f, Coarse Ratio=%2.0f x %2.0f\n', ...
                    f, obj.Length_AB, obj.Width_AD, obj.N_Length_AB, obj.N_Width_AD, obj.ADM(1), obj.ADM(2), obj.ADM(3), obj.ADM(4) );
        end
    end
end