%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel;
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FracturesGrid, CrossConnections] = scanfracture_data_uns(obj,ReservoirGrid)
    % Loading of basic settings
    fracture_file = fullfile(obj.SimulationInput.Directory, obj.SimulationInput.FracturesProperties.FractureMatrix);
    
    if strcmpi(obj.SimulationInput.FracturesProperties.FractureModelType,"EDFM")
        flag_edfm = true;
    else
        flag_edfm = false;
    end
    
    if strcmpi(obj.SimulationInput.FracturesProperties.FractureModelType,"pEDFM")
        flag_pedfm = true;
    else
        flag_pedfm = false;
    end
    
    % Creation of the rock frac mesh intersection settings
    rockfrac_int_blueprint = struct( ...
        "edfm_method",...
        obj.SimulationInput.FracturesProperties.all_settings.EDFM_METHOD, ...
        "num_area_points", ...
        obj.SimulationInput.FracturesProperties.all_settings.FRACTURE_NUM_POINTS,...
        "num_frac_boundbox", ...
        obj.SimulationInput.FracturesProperties.all_settings.FRACTURE_BOUNDING_BOX,...
        "parallel", ...
        obj.SimulatorSettings.parallel, ...
        "useHalf",...
        obj.SimulatorSettings.half_precision, ...
        "edfm",...
        flag_edfm, ...
        "pedfm", ...
        flag_pedfm, ...
        "d_avg_refinement", ...
        obj.SimulationInput.FracturesProperties.all_settings.FRACTURE_AVERAGE_DISTANCE_REFINEMENT...
        );

    % Translating fracture blueprint
    UnstructureData = obj.SimulationInput.ReservoirProperties.UnstructuredGridData;
    if isempty(UnstructureData.props)
        UnstructureData.props = struct("perm",obj.SimulationInput.ReservoirProperties.Perm);
    end
    
    % Generating fracture geometry blueprint from DARSim file
    [fPath, ~, ~] = fileparts(fracture_file);
    procFile = fullfile(fPath, "fracture_process.mat");
    
    if isfile(procFile)
        disp("Fractures and Intersections reloaded.");
        load(fullfile(fPath,"fracture_process.mat"),"fractures")
    else
        fractures = FractureFactory( ...
            obj.SimulationInput.FracturesProperties, ...
            UnstructureData, ...
            rockfrac_int_blueprint);
        save(fullfile(fPath,"fracture_process.mat"),"fractures")
    end
    
    NumOfFrac = obj.SimulationInput.FracturesProperties.NumOfFrac;
    fprintf('Extracting data from %02d fractures ...\n', NumOfFrac);

    FracturesGrid = fractures_grid(NumOfFrac);
    FracturesGrid.dfrac = fractures;
    
    % Translating between DARSim versions
    trans_obj = translate_newDarsim_to_oldDarsim(fractures.frac_core, ...
        obj.SimulationInput.ReservoirProperties.UnstructuredGridData, fractures.fracmesh);
    
    % Updating the number of fracture grid since on unstructured grids this cannot be load a priori
    obj.SimulationInput.FracturesProperties.NumOfAllFracGrids = sum(trans_obj.num_frac_el);
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
    
    % Creation of cartesian grid in each fracture
    for f = 1 : NumOfFrac
        FractureGrid = unstructured_fracture_grid(trans_obj.fracmesh(f));
        FractureGrid.GridCoords = FractureGrid.fracmesh.vertex.coords;
        FracturesGrid.AddGrid(FractureGrid, f);
        Nf(f) = FractureGrid.N;
    end

    % In darsim only inner faces are considered,therefore we are retrieving
    % only the inner faces. DARSIM 3 data structured uses [inner; bound];
    trans_alpha = fractures.frac_core.alpha;
    trans_alpha = trans_alpha(1:obj.SimulationInput.ReservoirProperties.UnstructuredGridData.face.inner.num);
    ReservoirGrid.AddpEDFMCorrections(trans_alpha);
    fprintf(' ---> Complete!\n');

end


