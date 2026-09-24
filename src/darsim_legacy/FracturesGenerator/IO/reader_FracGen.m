% Output writer base class
%--------------------------------------------------------------------------
% DARSim 2 Reservoir Simulator
% Author: Matteo Cusini
% TU Delft
% Created: 13 July 2016
% Last modified: 14 July 2022 by Túlio de M. Cavalcante (UFPE)
%--------------------------------------------------------------------------
classdef reader_FracGen < handle
    properties
        Directory
        File
        InputMatrix
        CornerPointGridMatrix
        TetrahedralGridMatrix
    end
    methods
        function obj = reader_FracGen(dir, file)
            obj.Directory = dir;
            obj.File = fullfile(dir,file);
        end
        function ReadInputFile(obj,Builder)
            %ReadInputFile
            fileID = fopen(obj.File, 'r');
            %// Read lines from input file
            matrix = textscan(fileID, '%s', 'Delimiter', '\n');
            obj.InputMatrix = matrix{1};
            fclose(fileID);
            
            Builder.MainInputData = obj.ReadInformation();
        end
        function MainInputData = ReadInformation(obj)
            %% INPUT FILE
            % Passinbg the main path of the input directory
            MainInputData.Directory = obj.Directory;
            
            % The type of the fracture model (EDFM or pEDFM)
            temp = strfind(obj.InputMatrix, 'TYPE');
            index = find(~cellfun('isempty', temp));
            MainInputData.Type = obj.InputMatrix{index+1};
            if ~strcmp(MainInputData.Type,'EDFM') && ~strcmp(MainInputData.Type,'pEDFM')
                error('The type mentioned in the input file should either be "EDFM" or "pEDFM"!');
            end
            
            % Reading reservoir properties from the input file
            MainInputData.ReservoirProperties = obj.ReadReservoirProperties();
            
            if strcmp(MainInputData.ReservoirProperties.Discretization,'TetrahedralGrid')
                MainInputData.NumOfDimensions = 3;
            else        
                if MainInputData.ReservoirProperties.Grid.N(3) > 1
                    MainInputData.NumOfDimensions = 3;
                else
                    MainInputData.NumOfDimensions = 2;
                end
            end
            
            % Reading fractures properties
            temp = strfind(obj.InputMatrix, 'FRAC_START');
            MainInputData.FracStartIndex = find(~cellfun('isempty', temp));
            temp = strfind(obj.InputMatrix, 'FRAC_END');
            MainInputData.FracEndIndex = find(~cellfun('isempty', temp));
            MainInputData.NrOfAllFracs = MainInputData.FracEndIndex - MainInputData.FracStartIndex - 1;
            MainInputData.FracturesProperties = obj.ReadFracturesProperties(MainInputData);
            
            % The precision of the comparison of two values (near zero aproximation)
            temp = strfind(obj.InputMatrix, 'CompareAccuracy');
            index = find(~cellfun('isempty', temp));
            MainInputData.compareAccuracy = str2double( obj.InputMatrix{index+1} );
        end
        function ReservoirProperties = ReadReservoirProperties(obj)
            % The type of the discretization of the reservoir (CartesianGrid or CornerPointGrid)
            temp = strfind(obj.InputMatrix, 'GRID_DISCRETIZATION');
            index = find(~cellfun('isempty', temp));
            if ~isempty(index)
                ReservoirProperties.Discretization = obj.InputMatrix{index+1};
            else
                warning('The keyword "GRID_DISCRETIZATION" is missing. Cartesian discretization is set by default.\n');
                ReservoirProperties.Discretization = 'CartesianGrid';
            end
            flag_cartesianbox = 0;
            % Check the discretization approach
            switch ReservoirProperties.Discretization
                case('CornerPointGrid')
                    % Read the cornerpointgrid data input file
                    FileName = obj.InputMatrix{index+2};
                    if isfile(FileName)
                        CornerPointGridFile = FileName;
                    elseif isfile(strcat(obj.Directory, '/', FileName))
                        CornerPointGridFile = strcat(obj.Directory, '/', FileName);
                    else %isfile(strcat(obj.Directory, '/','CornerPointGridData.mat'))
                        CornerPointGridFile = "CartesianBox";
                        flag_cartesianbox = 1;
                    % else
                    %     error('There is no CornerPointGrid input file given in the main input file. Please specify it correctly in the main input file.');
                    end
                    % Check if the data has already been loaded and saved, otherwise read the CornerPointGrid input file
                    if isfile(strcat(obj.Directory, '/','CornerPointGridData.mat'))
                        fprintf('"CornerPointGridData.mat" file already exists. No need to load the CornerPointGrid data input file.\n');
                        fprintf('---> Loading the "CornerPointGridData.mat" file ...');
                        load(strcat(obj.Directory, '/','CornerPointGridData.mat'),'ReservoirProperties');
                        fprintf('Done!\n');
                    else                        
                        EclipseReader = reader_eclipse(CornerPointGridFile);
                        if flag_cartesianbox 
                            CartesianBox = struct();
                            ref = find(ismember(obj.InputMatrix,"SPECGRID"));
                            grid_size = str2double(obj.InputMatrix(ref+1:ref+3));
                            ref = find(ismember(obj.InputMatrix,"DIMENS"));
                            dimen_size = str2double(obj.InputMatrix(ref+1:ref+3));                            
                            CartesianBox.Lx = dimen_size(1);
                            CartesianBox.Ly = dimen_size(2);
                            CartesianBox.Lz = dimen_size(3);
                            CartesianBox.Nx = grid_size(1);
                            CartesianBox.Ny = grid_size(2);
                            CartesianBox.Nz = grid_size(3);
                         ReservoirProperties.CornerPointGridData = ...
                                EclipseReader.ConvertCartesianBoxToCornerPointGridData(CartesianBox);                            
                        elseif contains(CornerPointGridFile,'.grdecl') && ~contains(CornerPointGridFile,'.txt')
                            % Reading the grdecl file directly
                            ReservoirProperties.CornerPointGridData = EclipseReader.ReadGRDECL(1);
                            if isfield(ReservoirProperties.CornerPointGridData,'Porosity')
                                ReservoirProperties.CornerPointGridData = rmfield(ReservoirProperties.CornerPointGridData,{'Porosity'});
                            end
                            if isfield(ReservoirProperties.CornerPointGridData,'Permeability')
                                ReservoirProperties.CornerPointGridData = rmfield(ReservoirProperties.CornerPointGridData,{'Permeability'});
                            end

                        elseif contains(CornerPointGridFile,'.txt')
                            % Reading the txt files generated by converting grdecl via DARSim CornerPointGrid data generator
                            ReservoirProperties.CornerPointGridData = EclipseReader.ReadfromTXT();
                        else
                            error('The extension of CornerPointGrid input file is not supported (must be txt or grdecl).');
                        end
                        % fprintf('---> Saving the "CornerPointGridData.mat" file ...');
                        % save(strcat(obj.Directory, '/','CornerPointGridData.mat'),'ReservoirProperties','-v7.3');
                        fprintf('Done!\n');
                    end

                    ReservoirProperties.Grid.N(1) = ReservoirProperties.CornerPointGridData.Nx;
                    ReservoirProperties.Grid.N(2) = ReservoirProperties.CornerPointGridData.Ny;
                    ReservoirProperties.Grid.N(3) = ReservoirProperties.CornerPointGridData.Nz;
                    ReservoirProperties.Grid.N_ActiveCells = ReservoirProperties.CornerPointGridData.N_ActiveCells;
                    
                    ReservoirProperties.Size = [ReservoirProperties.CornerPointGridData.Lx;...
                                                ReservoirProperties.CornerPointGridData.Ly;
                                                ReservoirProperties.CornerPointGridData.Lz];
                
                case('TetrahedralGrid')
                    % Read the tetrahedralgrid data input file
                    FileName = obj.InputMatrix{index+2};
                    if isfile(FileName)
                        TetrahedralGridFile = FileName;
                    elseif isfile(strcat(obj.Directory, '/', FileName))
                        TetrahedralGridFile = fullfile(obj.Directory, FileName);
                    else
                        error('There is no TetrahedralGrid input file given in the main input file. Please specify it correctly in the main input file.');
                    end

                    % Check if the data has already been loaded and saved, otherwise read the TetrahedralGrid input file
                    if isfile(strcat(obj.Directory, '/','TetrahedralGridData.mat'))
                        fprintf('"TetrahedralGridData.mat" file already exists. No need to load the TetrahedralGrid data input file.\n');
                        fprintf('---> Loading the "TetrahedralGridData.mat" file ...');
                        load(strcat(obj.Directory, '/','TetrahedralGridData.mat'),'ReservoirProperties');
                        fprintf('Done!\n');
                    else
                        EclipseReader = reader_eclipse(TetrahedralGridFile);
                        if contains(TetrahedralGridFile,'.grdecl') && ~contains(TetrahedralGridFile,'.txt')
                            % Reading the grdecl file directly
                            ReservoirProperties.TetrahedralGridData = EclipseReader.ReadGRDECL(1);
                            ReservoirProperties.TetrahedralGridData = rmfield(ReservoirProperties.TetrahedralGridData,{'Porosity','Permeability'});

                        elseif contains(TetrahedralGridFile,'.txt')
                            % Reading the txt files generated by converting grdecl via DARSim TetrahedralGrid data generator
                            ReservoirProperties.TetrahedralGridData = EclipseReader.ReadfromTXT();
                        else
                            error('The extension of TetrahedralGrid input file is not supported (must be txt or grdecl).');
                        end
                        fprintf('---> Saving the "TetrahedralGridData.mat" file ...');
                        save(strcat(obj.Directory, '/','TetrahedralGridData.mat'),'ReservoirProperties','-v7.3');
                        fprintf('Done!\n');
                    end

                    ReservoirProperties.Grid.N_ActiveCells = size(ReservoirProperties.TetrahedralGridData.Cell.Volume,1);
                    
                    ReservoirProperties.Size = [max(ReservoirProperties.TetrahedralGridData.Vertex.Coordinates(:,1));...
                                                max(ReservoirProperties.TetrahedralGridData.Vertex.Coordinates(:,2));
                                                max(ReservoirProperties.TetrahedralGridData.Vertex.Coordinates(:,3))];      
                                            
                    ReservoirProperties.Grid.N = ReservoirProperties.Size;
                                            
                case('CartesianGrid')
                    % Assume it is Cartesian
                    % 1. size of the reservoir
                    temp = strfind(obj.InputMatrix, 'DIMENS'); % Search a specific string and find all rows containing matches
                    index = find(~cellfun('isempty', temp));
                    if isempty(index)
                        error('The keyword "DIMENS" is missing. Please check the input file!\n');
                    end
                    ReservoirProperties.Size = [str2double(obj.InputMatrix{index+1});...
                                                str2double(obj.InputMatrix{index+2});
                                                str2double(obj.InputMatrix{index+3})];
                    % 2. GridSize
                    temp = strfind(obj.InputMatrix, 'SPECGRID');
                    index = find(~cellfun('isempty', temp));
                    if isempty(index)
                        error('The keyword "SPECGRID" is missing. Please check the input file!\n');
                    end
                    ReservoirProperties.Grid.N = [str2double(obj.InputMatrix{index+1});...
                                                  str2double(obj.InputMatrix{index+2});
                                                  str2double(obj.InputMatrix{index+3})];
                    ReservoirProperties.Grid.N_ActiveCells = prod(ReservoirProperties.Grid.N);
                otherwise
                    error('The discretization method should either be "CartesianGrid" or "CornerPointGrid". Check the input file!\n');
            end
            
            % The ADM settingss of the reservoir
            temp = strfind(obj.InputMatrix, 'ADM_RESERVOIR');
            index = find(~cellfun('isempty', temp));
            if isempty(index)
                ReservoirProperties.ADM = [0,0,1,1,1];
            else
                ReservoirProperties.ADM = strsplit(obj.InputMatrix{index+1},{' ' ,'','[',',',']'});
                if ReservoirProperties.ADM{1} == '0'
                    ReservoirProperties.ADM = [0,0,1,1,1];
                else                
                    ReservoirProperties.ADM = str2double( ReservoirProperties.ADM(1:5) );
                end
            end
            
            % The stress settings of the reservoir
            temp = strfind(obj.InputMatrix, 'RESERVOIR_STRESS');
            index = find(~cellfun('isempty', temp));
            if ~isempty(index)
                ReservoirProperties.Stress.Sigma1           = str2double( obj.InputMatrix{index+1} );
                ReservoirProperties.Stress.Sigma3           = str2double( obj.InputMatrix{index+2} );
                ReservoirProperties.Stress.Pressure         = str2double( obj.InputMatrix{index+3} );
                ReservoirProperties.Stress.Sigma1_Direction = str2double( obj.InputMatrix{index+4} );
            else
                ReservoirProperties.Stress = [];
            end
        end
        function FracturesProperties = ReadFracturesProperties(obj,MainInputData)
            FracturesProperties = struct();
            for f = 1 : MainInputData.NrOfAllFracs
                frac_temp = strsplit(obj.InputMatrix{MainInputData.FracStartIndex+f},' , ');
                %%% Check if fracture is used or not %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                FracturesProperties(f).isActive = str2double( frac_temp{01} );
                
                %%% Central Coordinate of the fracture plate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                frac_temp{02}       = regexprep(frac_temp{02},' ' ,'');
                CenterCoord_temp    = strsplit(frac_temp{02}, { '[' , ';' , ']' });
                CenterCoord_temp{2} = strsplit(CenterCoord_temp{2}, { '*' , '/' });
                CenterCoord_temp{3} = strsplit(CenterCoord_temp{3}, { '*' , '/' });
                CenterCoord_temp{4} = strsplit(CenterCoord_temp{4}, { '*' , '/' });
                if sum( strcmp('LX' , CenterCoord_temp{2}) ) > 0
                    LX_Index = find( strcmp('LX' , CenterCoord_temp{2}) == 1 );
                    CenterCoord_temp{2}{LX_Index} = num2str( MainInputData.ReservoirProperties.Size(1) );
                end
                if length(CenterCoord_temp{2}) == 1
                    FracturesProperties(f).CenterCoord(1,1) = str2double(CenterCoord_temp{2});
                else
                    FracturesProperties(f).CenterCoord(1,1) = str2double(CenterCoord_temp{2}{1}) * str2double(CenterCoord_temp{2}{2}) / str2double(CenterCoord_temp{2}{3});
                end
                if sum( strcmp('LY' , CenterCoord_temp{3}) ) > 0
                    LY_Index = find( strcmp('LY' , CenterCoord_temp{3}) == 1);
                    CenterCoord_temp{3}{LY_Index} = num2str( MainInputData.ReservoirProperties.Size(2) );
                end
                if length(CenterCoord_temp{3}) == 1
                    FracturesProperties(f).CenterCoord(2,1) = str2double(CenterCoord_temp{3});
                else
                    FracturesProperties(f).CenterCoord(2,1) = str2double(CenterCoord_temp{3}{1}) * str2double(CenterCoord_temp{3}{2}) / str2double(CenterCoord_temp{3}{3});
                end
                if sum( strcmp('LZ' , CenterCoord_temp{4}) ) > 0
                    LZ_Index = find( strcmp('LZ' , CenterCoord_temp{4}) == 1);
                    CenterCoord_temp{4}{LZ_Index} = num2str( MainInputData.ReservoirProperties.Size(3) );
                end
                if length(CenterCoord_temp{4}) == 1
                    FracturesProperties(f).CenterCoord(3,1) = str2double(CenterCoord_temp{4});
                else
                    FracturesProperties(f).CenterCoord(3,1) = str2double(CenterCoord_temp{4}{1}) * str2double(CenterCoord_temp{4}{2}) / str2double(CenterCoord_temp{4}{3});
                end
                
                %%% Length of the fracture plate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Length_temp = strsplit(frac_temp{03}, { '*' , '/' });
                if sum( strcmp('LX' , Length_temp) ) > 0
                    LX_Index = find( strcmp('LX' , Length_temp) == 1);
                    Length_temp{LX_Index} = num2str( MainInputData.ReservoirProperties.Size(1) );
                end
                if length(Length_temp) == 1
                    FracturesProperties(f).Length = str2double(Length_temp);
                else
                    FracturesProperties(f).Length = str2double(Length_temp{1}) * str2double(Length_temp{2}) / str2double(Length_temp{3});
                end
                if sum( strcmp('LY' , Length_temp) ) > 0
                    LY_Index = find( strcmp('LY' , Length_temp) == 1);
                    Length_temp{LY_Index} = num2str( MainInputData.ReservoirProperties.Size(2) );
                end
                if length(Length_temp) == 1
                    FracturesProperties(f).Length = str2double(Length_temp);
                else
                    FracturesProperties(f).Length = str2double(Length_temp{1}) * str2double(Length_temp{2}) / str2double(Length_temp{3});
                end
                if sum( strcmp('LZ' , Length_temp) ) > 0
                    LZ_Index = find( strcmp('LZ' , Length_temp) == 1);
                    Length_temp{LZ_Index} = num2str( MainInputData.ReservoirProperties.Size(3) );
                end
                if length(Length_temp) == 1
                    FracturesProperties(f).Length = str2double(Length_temp);
                else
                    FracturesProperties(f).Length = str2double(Length_temp{1}) * str2double(Length_temp{2}) / str2double(Length_temp{3});
                end
                
                %%% Width of the fracture plate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Width_temp = strsplit(frac_temp{04}, { '*' , '/' });
                if sum( strcmp('LX' , Width_temp) ) > 0
                    LX_Index = find( strcmp('LX' , Width_temp) == 1);
                    Width_temp{LX_Index} = num2str( MainInputData.ReservoirProperties.Size(1) );
                end
                if length(Width_temp) == 1
                    FracturesProperties(f).Width = str2double(Width_temp);
                else
                    FracturesProperties(f).Width = str2double(Width_temp{1}) * str2double(Width_temp{2}) / str2double(Width_temp{3});
                end
                if sum( strcmp('LY' , Width_temp) ) > 0
                    LY_Index = find( strcmp('LY' , Width_temp) == 1);
                    Width_temp{LY_Index} = num2str( MainInputData.ReservoirProperties.Size(2) );
                end
                if length(Width_temp) == 1
                    FracturesProperties(f).Width = str2double(Width_temp);
                else
                    FracturesProperties(f).Width = str2double(Width_temp{1}) * str2double(Width_temp{2}) / str2double(Width_temp{3});
                end
                if sum( strcmp('LZ' , Width_temp) ) > 0
                    LZ_Index = find( strcmp('LZ' , Width_temp) == 1);
                    Width_temp{LZ_Index} = num2str( MainInputData.ReservoirProperties.Size(3) );
                end
                if length(Width_temp) == 1
                    FracturesProperties(f).Width = str2double(Width_temp);
                else
                    FracturesProperties(f).Width = str2double(Width_temp{1}) * str2double(Width_temp{2}) / str2double(Width_temp{3});
                end
                
                %%% Rotation angles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Rotation angle along Z
                FracturesProperties(f).RotAngleAlongZ = str2double( frac_temp{05} );
                
                % Rotation angle along central length
                FracturesProperties(f).RotAngleAlongL = str2double( frac_temp{06} );

                % Rotation angle along central width
                FracturesProperties(f).RotAngleAlongW = str2double( frac_temp{07} );
                
                % If the reservoir is in 2D, force the fracture's plate to be vertical (act as a line)
                if strcmp(MainInputData.ReservoirProperties.Discretization,'TetrahedralGrid')~=1 && ...
                   MainInputData.ReservoirProperties.Grid.N(3) == 1
                    FracturesProperties(f).RotAngleAlongL = 90;
                    FracturesProperties(f).RotAngleAlongW = 0;
                end

                %%% Grid resolution ratio between fracture and matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                FracturesProperties(f).GridResRatio   = str2double( frac_temp{08} );
                
                %%% Grid numbers along the length of the fracture plate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                FracturesProperties(f).GridNumAlongL  = str2double( frac_temp{09} );
                
                %%% Grid numbers along the width of the fracture plate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                FracturesProperties(f).GridNumAlongW  = str2double( frac_temp{10} );
                
                %%% The aperture of the fracture plate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                FracturesProperties(f).Aperture       = str2double( frac_temp{11} );
                
                %%% The porosity of the fracture plate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                FracturesProperties(f).Porosity       = str2double( frac_temp{12} );
                
                %%% The permeability of the fracture plate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                FracturesProperties(f).Permeability   = str2double( frac_temp{13} );

                %%% ADM Configuration of the fracture %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                frac_temp{14}                         = regexprep(frac_temp{14},' ' ,'');
                ADM_temp                              = strsplit(frac_temp{14}, { '[' , ',' , ']' });
                FracturesProperties(f).ADM            = [ str2double(ADM_temp{2}) , str2double(ADM_temp{3}) , str2double(ADM_temp{4}) , str2double(ADM_temp{5}) ];
                %%% The coordinates of 3 corners of fracture plate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                frac_temp{15}                         = regexprep(frac_temp{15},' ' ,'');
                CornerCoord_temp                      = strsplit(frac_temp{15}, { '[' , ';' , ',' , ']' });
                FracturesProperties(f).CornerCoords   = [ str2double(CornerCoord_temp{02}) , str2double(CornerCoord_temp{03}) , str2double(CornerCoord_temp{04})
                                                          str2double(CornerCoord_temp{05}) , str2double(CornerCoord_temp{06}) , str2double(CornerCoord_temp{07})
                                                          str2double(CornerCoord_temp{08}) , str2double(CornerCoord_temp{09}) , str2double(CornerCoord_temp{10}) ];
            end
            1
        end
    end
end