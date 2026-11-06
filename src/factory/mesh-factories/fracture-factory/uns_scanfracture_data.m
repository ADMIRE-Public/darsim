%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Ryan Haagenson; Artur Castiel
%TU Delft
%
% Scan for previous FractureFactory output file. If none exists, run
% FractureFactory() to generate (p)EDFM connections.
%
% Input:
%   ReservoirGrid:  struct
%                   Contains reservoir grid information
%
% Output:
%   FracturesGrid:    struct
%                     Contains fracture grid information
%   CrossConnections: struct
%                     Contains (p)EDFM cross connection information for 
%                     each fracture cell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FracturesGrid, CrossConnections] = ...
         uns_scanfracture_data(obj,ReservoirGrid)
    
    %%%%%%%%%%%%%%%%%%
    % Initialization %
    %%%%%%%%%%%%%%%%%%

    % Loading of basic settings
    FractureFile = fullfile(obj.SimulationInput.Directory, ...
        obj.SimulationInput.FracturesProperties.FractureMatrix);
    
    % Set the EDFM/pEDFM flags
    if strcmpi(...
            obj.SimulationInput.FracturesProperties.FractureModelType,...
            "EDFM")
        FlagEDFM = true;
    else
        FlagEDFM = false;
    end
    
    if strcmpi(...
            obj.SimulationInput.FracturesProperties.FractureModelType,...
            "pEDFM")
        FlagPEDFM = true;
    else
        FlagPEDFM = false;
    end
    
    % Create rock-frac mesh intersection settings blueprint
    % These are set by user-input (SimulationInput or SimulatorSettings)
    RockFracIntBlueprint = struct( ...
        "edfm",...
        FlagEDFM, ...
        "pedfm", ...
        FlagPEDFM);

    % Unstructure rock mesh properties
    UnstructureData = ...
        obj.SimulationInput.ReservoirProperties.UnstructuredGridData;
    if isempty(UnstructureData.props)
        UnstructureData.props = ...
            struct("perm",obj.SimulationInput.ReservoirProperties.Perm);
    end
    
    % MATLAB file to save FractureFactory output
    [fPath, ~, ~] = fileparts(FractureFile);
    ProcFile = fullfile(fPath, "fracture_process.mat");

    nFrac = obj.SimulationInput.FracturesProperties.NumOfFrac;
    fprintf('Extracting data from %02d fractures ...\n', nFrac);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Creating CrossConnections %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Read in saved FractureFactory output if available
    if isfile(ProcFile)
        disp("Fractures and Intersections reloaded.");
        load(fullfile(fPath,"fracture_process.mat"),"Fractures")
    % Else call FractureFacotry and save output
    else
        Fractures = FractureFactory( ...
            obj.SimulationInput.FracturesProperties, ...
            UnstructureData, ...
            RockFracIntBlueprint);
        save(fullfile(fPath,"fracture_process.mat"),"Fractures")
    end
    
    % RH: Artur created the Fractures object, but this isn't actually used
    % by DARSim. We need to translate this to the CrossConnections object
    % that is actually used by the flow and transport code of DARSim.
    TransObj = ...
        translate_Fractures_to_CrossConnections(Fractures.FracCore, ...
        obj.SimulationInput.ReservoirProperties.UnstructuredGridData, ...
        Fractures.FracMesh);

    % Creating the CrossConnections object
    % This is used to build Jacobian 
    % Same object exists for Cartesian- and CPG-based pEDFM method
    CrossConnections = cross_connections();
    for index = 1:length(TransObj.CrossConnections)
        CrossConnections(index).Cells = ...
            TransObj.CrossConnections(index).Cells;
        CrossConnections(index).CI = ...
            TransObj.CrossConnections(index).CI;
        CrossConnections(index).T_Geo = ...
            TransObj.CrossConnections(index).T_Geo;
        CrossConnections(index).T_Geo_Cond = ...
            TransObj.CrossConnections(index).T_Geo_Cond;
        CrossConnections(index).UpWind = ...
            TransObj.CrossConnections(index).UpWind;
        CrossConnections(index).U_Geo = ...
            TransObj.CrossConnections(index).U_Geo;
    end 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Creating FracturesGrid %
    %%%%%%%%%%%%%%%%%%%%%%%%%%

    FracturesGrid = fractures_grid(nFrac);
    FracturesGrid.dfrac = Fractures;
    
    % Creating a Cartesian grid for each fracture
    for iFrac = 1 : nFrac
        FractureGrid = unstructured_fracture_grid(TransObj.FracMesh(iFrac));
        FractureGrid.GridCoords = FractureGrid.fracmesh.vertex.coords;
        FracturesGrid.AddGrid(FractureGrid, iFrac);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Apply Alpha values to rock and frac meshes %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % DARim considers only inner faces; retrieving only the inner faces.
    TransAlpha = Fractures.FracCore.RockAlpha;
    TransAlpha = TransAlpha(1:Fractures.RockMesh.face.inner.num);
    ReservoirGrid.AddpEDFMCorrections(TransAlpha);

    for iFrac = 1:nFrac
        TransAlpha = Fractures.FracCore.FracAlpha{iFrac};
        TransAlpha = ...
            TransAlpha(1:Fractures.FracMesh(iFrac).faces.inner.num);
        FracturesGrid.Grids(iFrac).AddpEDFMCorrections(TransAlpha);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating total number of frac cells in system %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating the number of fracture grid cells tracked in SimulationInput,
    % since this cannot be load a priori with unstructrued grids
    obj.SimulationInput.FracturesProperties.NumOfAllFracGrids = ...
        sum(TransObj.nFracElem);

    fprintf(' ---> Complete!\n');

end


