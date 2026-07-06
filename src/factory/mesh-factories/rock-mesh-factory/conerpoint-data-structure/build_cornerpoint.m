function CornerPointGridData = build_cornerpoint(blueprint, folder)
% Read the cornerpointgrid data input file
if ischar(blueprint.GRID.GRID_DISCRETIZATION)
    CornerPointGridFile ='CartesianBox';
    includeCornerPointGridRockProperties = false;
else
    CornerPointGridFile = blueprint.GRID.GRID_DISCRETIZATION{2};
    if ~isfile(CornerPointGridFile)
        CornerPointGridFile = fullfile(folder, CornerPointGridFile);
    end
    if size(CornerPointGridFile,1) > 2
        if contains(blueprint.GRID.GRID_DISCRETIZATION{end}, ...
                "INCLUDE_ROCK_PROPERTIES")
            includeCornerPointGridRockProperties = true;
        else
            includeCornerPointGridRockProperties = false;
        end
    end
end

if ~isfile(fullfile(folder,"CornerPointGridData.mat"))
    EclipseReader = reader_eclipse(CornerPointGridFile);
    if contains(CornerPointGridFile, 'CartesianBox', "IgnoreCase",true)
        CartesianBox = struct();
        CartesianBox.Lx = blueprint.GENERAL.DIMENSION(1)';
        CartesianBox.Ly = blueprint.GENERAL.DIMENSION(2)';
        CartesianBox.Lz = blueprint.GENERAL.DIMENSION(3)';
        CartesianBox.Nx = blueprint.GRID.SPECGRID(1);
        CartesianBox.Ny = blueprint.GRID.SPECGRID(2);
        CartesianBox.Nz = blueprint.GRID.SPECGRID(3);
        CornerPointGridData = ...
            EclipseReader.ConvertCartesianBoxToCornerPointGridData(CartesianBox);

    elseif contains(CornerPointGridFile, '.txt', "IgnoreCase",true)
        CornerPointGridData =...
            EclipseReader.ReadfromTXT(CornerPointGridFile);
    elseif contains(CornerPointGridFile, '.grdecl', "IgnoreCase",true)
        CornerPointGridData = EclipseReader.ReadGRDECL(1);
    else
        error("CPG inputfile not recognized.")
    end
    save(fullfile(folder,"CornerPointGridData.mat"), "CornerPointGridData");
else
    load(fullfile(folder,"CornerPointGridData.mat"), "CornerPointGridData");
end

% if ~isempty(blueprint.GRID.GRID_DISCRETIZATION{3,1})
if size(blueprint.GRID.GRID_DISCRETIZATION, 1) == 3
    includeCornerPointGridRockProperties = 1;
else
    includeCornerPointGridRockProperties = 0;
end

% Deleting the permeability data if simulation is set to be homogeneous
% if ~includeCornerPointGridRockProperties
%     try
%         CornerPointGridData = rmfield(CornerPointGridData,{'Porosity', ...
%             'Permeability','PermUnit','PermScale'});
%     catch
%         error("Error while trying to remove CPG data.")
%     end
% end

end