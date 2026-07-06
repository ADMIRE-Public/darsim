function FracturesProperties = ReadFracturesPropertiesUns(obj)
%%%%%%%%%%%%%PROPERTIES OF THE FRACTURES%%%%%%%%%%%%%%%%
% Marking the simulation as fractured
FracturesProperties.isFractured = 1;
fracture_geometry_blueprint = FracReader(obj.FractureMatrix);
FracturesProperties.FractureModelType = fracture_geometry_blueprint.blueprint.Type;
FracturesProperties.DiscretizationType = fracture_geometry_blueprint.blueprint.Discretization;
FracturesProperties.Reservoir.Lx = fracture_geometry_blueprint.blueprint.Size(1);
FracturesProperties.Reservoir.Ly = fracture_geometry_blueprint.blueprint.Size(2);
FracturesProperties.Reservoir.Lz = fracture_geometry_blueprint.blueprint.Size(3);
FracturesProperties.Reservoir.Nx = fracture_geometry_blueprint.blueprint.N(1);
FracturesProperties.Reservoir.Ny = fracture_geometry_blueprint.blueprint.N(2);
FracturesProperties.Reservoir.Nz = fracture_geometry_blueprint.blueprint.N(3);
temp = strsplit(fracture_geometry_blueprint.blueprint.ADM,{' ',',','[',']'});
FracturesProperties.ReservoirADM.isActive = str2double(temp{2});
FracturesProperties.ReservoirADM.CoarseLevel = str2double(temp{3});
FracturesProperties.ReservoirADM.CoarseRatio = [ str2double(temp{4}) ; str2double(temp{5}) ; str2double(temp{6}) ];
FracturesProperties.Reservoir.N = FracturesProperties.Reservoir.Nx * FracturesProperties.Reservoir.Ny *FracturesProperties.Reservoir.Nz;
FracturesProperties.NumOfFrac = sum([fracture_geometry_blueprint.blueprint.FracturesProperties.isActive]);
FracturesProperties.NumOfAllFracGrids = length([fracture_geometry_blueprint.blueprint.FracturesProperties.isActive]);

% Reading the number of fractures
% temp = strfind(obj.FractureMatrix, 'NUM_FRACS');
% index = find(~cellfun('isempty', temp));
% if isempty(index)
%     error('DARSim Error: The keyword "NUM_FRACS" is missing. Please check the Fracture_Output file!\n');
% end
% temp = strsplit(obj.FractureMatrix{index},' ');
% FracturesProperties.NumOfFrac = str2double( temp{2} );

% Reading the number of all fractures grids
% temp = strfind(obj.FractureMatrix, 'NUM_TOTAL_FRACS_GRIDS');
% index = find(~cellfun('isempty', temp));
% if isempty(index)
%     error('DARSim Error: The keyword "NUM_TOTAL_FRACS_GRIDS" is missing. Please check the Fracture_Output file!\n');
% end
% temp = strsplit(obj.FractureMatrix{index},' ');
% FracturesProperties.NumOfAllFracGrids = str2double( temp{2} );
end