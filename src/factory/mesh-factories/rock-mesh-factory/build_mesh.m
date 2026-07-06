function ReservoirProperties = ...
                build_mesh(ReservoirProperties, blueprint, folder_path)
ReservoirProperties.Discretization =  ...
    get_first(blueprint.GRID.GRID_DISCRETIZATION);

switch ReservoirProperties.Discretization 
    case "CartesianGrid"
       ReservoirProperties.size = blueprint.GENERAL.DIMENSION;
       Grid = build_cartesian(blueprint.GRID);
       ReservoirProperties.Grid = Grid;
    case "CornerPointGrid"
       ReservoirProperties.CornerPointGridData  = ...
           build_cornerpoint(blueprint, folder_path);
       ReservoirProperties.Grid.N_ActiveCells = ...
           ReservoirProperties.CornerPointGridData.N_ActiveCells;
       ReservoirProperties.Grid.N(1,1) = ReservoirProperties.CornerPointGridData.Nx;
       ReservoirProperties.Grid.N(2,1) = ReservoirProperties.CornerPointGridData.Ny;
       ReservoirProperties.Grid.N(3,1) = ReservoirProperties.CornerPointGridData.Nz;
       ReservoirProperties.size(1) = ...
                        ReservoirProperties.CornerPointGridData.Lx;
       ReservoirProperties.size(2) = ...
                        ReservoirProperties.CornerPointGridData.Ly;
       ReservoirProperties.size(3) = ...
                        ReservoirProperties.CornerPointGridData.Lz;              
    case "UnstructuredGrid"
        Grid = build_cartesian(blueprint.GRID);
        ReservoirProperties.size = blueprint.GENERAL.DIMENSION;
        [ReservoirProperties.UnstructuredGridData, ...
         ReservoirProperties.Grid, ReservoirProperties.size, ...
         ReservoirProperties.Unsbox] = build_unstructured(blueprint.GRID,...
         folder_path, Grid, ReservoirProperties.size);
    case "TetrahedralGrid"
       ReservoirProperties.TetrahedralGridData = ...
           build_tetrahedron(blueprint.GRID, folder_path);
       ReservoirProperties.size(:) = nan;
end
end


function grid = build_cartesian(blueprint)
    grid = struct();
    grid.N = blueprint.SPECGRID;
    grid.N_ActiveCells = prod(grid.N);
end


function [UnstructuredGridData,Grid, Size,UnsBox] = ...
    build_unstructured(mesh_bp, folder_path, Grid, Size)
if ischar(mesh_bp.GRID_DISCRETIZATION)
    UnsBox = true;
    meshcore = create_unsbox(Grid.N, Size);
else
    FileName = mesh_bp.GRID_DISCRETIZATION{2};
    FilePath = fullfile(folder_path,FileName);
    [meshcore] = manage_extension(FilePath);
    UnsBox = false;
    Grid.N(:) = NaN;
end
Grid.N_ActiveCells = length(meshcore.element.volume);
UnstructuredGridData = meshcore;
size_total = max(UnstructuredGridData.vertex.coord) - ...
    min(UnstructuredGridData.vertex.coord);
Size = size_total;
end

function [TetrahedralGridData] = build_tetrahedron(mesh_bp, folder_path)
matfile = fullfile(folder_path,'TetrahedralGridData.mat');
if isfile(matfile)
    fprintf('"TetrahedralGridData.mat" already exists. Reloading file.\n');
    fprintf('---> Loading the "TetrahedralGridData.mat" file... ');
    load(matfile, 'TetrahedralGridData');
    fprintf('Done!\n');
else
    % Read the cornerpointgrid data input file
    FileName = mesh_bp.GRID_DISCRETIZATION{2};
    tetrahedron_file = fullfile(folder_path, FileName);
    if ~isfile(tetrahedron_file)
        error("ERROR: %s not found", tetrahedron_file);
    end
    TetraReader = reader_tetrahedral_grid_file(tetrahedron_file);
    TetrahedralGridData = TetraReader.ReadMSH;
    fprintf('---> Saving the "TetrahedralGridData.mat" file... ');
    save(fullfile(folder_path, 'TetrahedralGridData.mat'),'TetrahedralGridData','-v7.3');
    fprintf('Done!\n');
end
end

function out = get_first(cell_data)
       out = convertCharsToStrings(horzcat(cell_data));
       out = out(1);
end
