function FracturesProperties = legacy_fracture_reader(folder, file_path)
%%%%%%%%%%%%%PROPERTIES OF THE FRACTURES%%%%%%%%%%%%%%%%
% Marking the simulation as fractured
file_path = fullfile(folder,file_path);
FractureMatrix = legacy_fracmatrix_generator(file_path);
FracturesProperties.FractureMatrix = FractureMatrix;
FracturesProperties.isFractured = true;
FracturesProperties.FilePath = file_path;
% Reading the discretization approach used in fracture generator
temp = strfind(FractureMatrix, 'DISCRETIZATION');
index = find(~cellfun('isempty', temp));
if isempty(index)
    error('DARSim Error: The keyword "DISCRETIZATION" is missing. Please check the Fracture_Output file!\n');
end
temp = strsplit(FractureMatrix{index},' ');
FracturesProperties.DiscretizationType = temp{end};

% Reading the fracture model approach (EDFM or pEDFM)
temp = strfind(FractureMatrix, 'FRACTURE_MODEL_TYPE');
index = find(~cellfun('isempty', temp));
if isempty(index)
    error('DARSim Error: The keyword "FRACTURE_MODEL_TYPE" (for fracture model approach EDFM or pEDFM) is missing. Please check the Fracture_Output file!\n');
end
temp = strsplit(FractureMatrix{index},' ');
FracturesProperties.FractureModelType = temp{end};

% Reading the reservoir size (Lx,Ly,Lz [m3]) from the fracture file
temp = strfind(FractureMatrix, 'DIMENSION');
index = find(~cellfun('isempty', temp));
if isempty(index)
    error('DARSim Error: The keyword "DIMENSION" is missing. Please check the Fracture_Output file!\n');
end
temp = strsplit(FractureMatrix{index},' ');
FracturesProperties.Reservoir.Lx = str2double(temp{2});
FracturesProperties.Reservoir.Ly = str2double(temp{4});
FracturesProperties.Reservoir.Lz = str2double(temp{6});

% Reading the reservoir grid (Nx,Ny,Nz) from the fracture file
temp = strfind(FractureMatrix, 'RESERVOIR_GRID');
index = find(~cellfun('isempty', temp));
if isempty(index)
    error('DARSim Error: The keyword "RESERVOIR_GRID" is missing. Please check the Fracture_Output file!\n');
end
temp = strsplit(FractureMatrix{index},' ');
FracturesProperties.Reservoir.Nx = str2double(temp{2});
FracturesProperties.Reservoir.Ny = str2double(temp{4});
FracturesProperties.Reservoir.Nz = str2double(temp{6});

% Reading the number of active grids reservoir grids from the fracture file
temp = strfind(FractureMatrix, 'RESERVOIR_ACTIVE_GRID');
index = find(~cellfun('isempty', temp));
if isempty(index)
    error('DARSim Error: The keyword "RESERVOIR_ACTIVE_GRID" is missing. Please check the Fracture_Output file!\n');
end
temp = strsplit(FractureMatrix{index},' ');
FracturesProperties.Reservoir.N = str2double(temp{2});

% Reading the ADM settings of the reservoir from the fracture file
temp = strfind(FractureMatrix, 'RESERVOIR_ADM');
index = find(~cellfun('isempty', temp));
if isempty(index)
    warning('DARSim Warning: The keyword "RESERVOIR_ADM" is missing in the Fracture_Output file.\nThe default values will be set.');
    FracturesProperties.ReservoirADM.isActive = 0;
    FracturesProperties.ReservoirADM.CoarseLevel = 1;
    FracturesProperties.ReservoirADM.CoarseRatio = [ 1 ; 1 ; 1];
end
temp = strsplit(FractureMatrix{index},{' ',',','[',']'});
FracturesProperties.ReservoirADM.isActive = str2double(temp{2});
FracturesProperties.ReservoirADM.CoarseLevel = str2double(temp{3});
FracturesProperties.ReservoirADM.CoarseRatio = [ str2double(temp{4}) ; str2double(temp{5}) ; str2double(temp{6}) ];

% Reading the number of fractures
temp = strfind(FractureMatrix, 'NUM_FRACS');
index = find(~cellfun('isempty', temp));
if isempty(index)
    error('DARSim Error: The keyword "NUM_FRACS" is missing. Please check the Fracture_Output file!\n');
end
temp = strsplit(FractureMatrix{index},' ');
FracturesProperties.NumOfFrac = str2double( temp{2} );

% Reading the number of all fractures grids
temp = strfind(FractureMatrix, 'NUM_TOTAL_FRACS_GRIDS');
index = find(~cellfun('isempty', temp));
if isempty(index)
    error('DARSim Error: The keyword "NUM_TOTAL_FRACS_GRIDS" is missing. Please check the Fracture_Output file!\n');
end
temp = strsplit(FractureMatrix{index},' ');
FracturesProperties.NumOfAllFracGrids = str2double( temp{2} );
end