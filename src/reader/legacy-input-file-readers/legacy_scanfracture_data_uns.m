function [FracturesGrid, CrossConnections] = legacy_scanfracture_data_uns(obj,ReservoirGrid)
%% simulating loading of basic settings
settings = struct();
settings.parallel = false;
settings.useHalf = false;
%% verifying availability
if ~exist('half','class')
    settings.useHalf = settings.useHalf & false;
end
%% simulating creation of the rock frac mesh intersection settings
rockfrac_int_blueprint = struct();
rockfrac_int_blueprint.num_area_points = 70;
rockfrac_int_blueprint.num_frac_boundbox = 7;
rockfrac_int_blueprint.num_frac_boundbox = 2;
rockfrac_int_blueprint.parallel = settings.parallel;
rockfrac_int_blueprint.useHalf = settings.useHalf;
rockfrac_int_blueprint.num_area_points = 120;
rockfrac_int_blueprint.num_frac_boundbox = 7;
rockfrac_int_blueprint.num_frac_boundbox = 1;

rockfrac_int_blueprint.edfm = true;
rockfrac_int_blueprint.pedfm = true;
%% intialize parallel pool
if settings.parallel
    p = gcp('nocreate');
    if isempty(p)
        p = parpool(maxNumCompThreads);
    end
end
%% loading input file
% new fracture generator
%addpath(genpath(pwd))
%case_input = fullfile(rFolder,"Fracture_Input_Tetrahedral3.txt");
case_input = FractureMatrix;
%% generating fracture geometry blueprint from DARSim file
[fPath, ~, ~] = fileparts(FractureMatrix);
if isfile(fullfile(fPath,"fracture_process.mat"))
    fracture_geometry_blueprint = FracReader(case_input);
    fractures = FractureFactory(fracture_geometry_blueprint, obj.SimulationInput.ReservoirProperties.UnstructuredGridData, ...
        rockfrac_int_blueprint);
    save(fullfile(fPath,"fracture_process.mat"),"fractures")
else
    disp("Fractures and Intersections reloaded.")
    load(fullfile(fPath,"fracture_process.mat"),"fractures")
end
NumOfFrac = obj.SimulationInput.FracturesProperties.NumOfFrac;
fprintf('Extracting data from %02d fractures ...\n', NumOfFrac);
FracturesGrid = fractures_grid(NumOfFrac);
FracturesGrid.dfrac = fractures;
trans_obj = translate_newDarsim_to_oldDarsim(fractures.frac_core, ...
    obj.SimulationInput.ReservoirProperties.UnstructuredGridData, fractures.fracmesh);
%% updating the number of fracture grid since on unstructured grids this cannot be load a priori
obj.SimulationInput.FracturesProperties.NumOfAllFracGrids = sum(trans_obj.num_frac_el);
%CrossConnections = trans_obj.CrossConnections;
CrossConnections = cross_connections();
for index = 1:length(trans_obj.CrossConnections)
    CrossConnections(index).Cells = trans_obj.CrossConnections(index).Cells;
    CrossConnections(index).CI = trans_obj.CrossConnections(index).CI;
    CrossConnections(index).T_Geo = trans_obj.CrossConnections(index).T_Geo;
    CrossConnections(index).T_Geo_Cond = trans_obj.CrossConnections(index).T_Geo_Cond;
    CrossConnections(index).UpWind = trans_obj.CrossConnections(index).UpWind;
    CrossConnections(index).U_Geo = trans_obj.CrossConnections(index).U_Geo;
end

Nf = zeros(NumOfFrac,1);

warning("creation of cartesian grid in each fracture is done here")
for f = 1 : NumOfFrac
    FractureGrid = unstructured_fracture_grid(trans_obj.fracmesh(f));
    FractureGrid.GridCoords = FractureGrid.fracmesh.vertex.coords;
    FracturesGrid.AddGrid(FractureGrid, f);
    Nf(f) = FractureGrid.N;
end
%% on darsim only inner faces are considered,therefore we are retrieving 
% only the inner faces. DARSIM 3 data structured uses [inner; bound];
trans_alpha = fractures.frac_core.alpha;
trans_alpha = trans_alpha(1:obj.SimulationInput.ReservoirProperties.UnstructuredGridData.face.inner.num);
ReservoirGrid.AddpEDFMCorrections(trans_alpha);

fprintf(' ---> Complete!\n');
end

