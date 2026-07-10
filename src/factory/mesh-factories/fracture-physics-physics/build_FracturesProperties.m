function FracturesProperties = build_FracturesProperties(blueprint, ...
                                        fracfile_path)
    if isempty(blueprint)
        FracturesProperties = struct("isFractured", 0, "NumOfFrac", 0);
    elseif blueprint.legacy
        FracturesProperties = blueprint;
        FracturesProperties.legacy = blueprint.legacy;
    elseif ~blueprint.legacy
        FracturesProperties = struct();
        FracturesProperties.legacy = blueprint.legacy;
        FracturesProperties.isFractured = 1;
        FracturesProperties.FractureModelType = blueprint.TYPE;
        FracturesProperties.UnsFracmesh = ...
            blueprint.UNSTRUCTURED_FRACTURE_MESH;        
        discretizationType = blueprint.GRID_DISCRETIZATION;
        if ischar(discretizationType) || isstring(discretizationType)
            FracturesProperties.DiscretizationType = ...
                blueprint.GRID_DISCRETIZATION;
        else
            FracturesProperties.DiscretizationType = ...
                discretizationType{1};
        end        
        FracturesProperties.Reservoir.Lx = blueprint.DIMENS(1);
        FracturesProperties.Reservoir.Ly = blueprint.DIMENS(2);
        FracturesProperties.Reservoir.Lz = blueprint.DIMENS(3);
        FracturesProperties.Reservoir.Nx = blueprint.SPECGRID(1);
        FracturesProperties.Reservoir.Ny = blueprint.SPECGRID(2);
        FracturesProperties.Reservoir.Nz = blueprint.SPECGRID(3);
        [FracturesProperties.ReservoirADM.isActive, ...
        FracturesProperties.ReservoirADM.CoarseLevel, ...
        FracturesProperties.ReservoirADM.CoarseRatio] = ...
                                split_adm_input(blueprint.ADM_RESERVOIR);
        frac_active = str2double(blueprint.FRAC(:,1));
        FracturesProperties.NumOfFrac =numel(frac_active);
        FracturesProperties.NumOfAllFracGrids = length(frac_active);
        FracturesProperties.FractureMatrix = fracfile_path{2};
        FracturesProperties.all_settings = blueprint;
        FracturesProperties.blueprint = blueprint;
        FracturesProperties.blueprint = ...
            generate_blueprint_description(FracturesProperties);

    end

end

function fractures_description = ...
        generate_fracture_description(FRAC_input,Lsize, Nsize)    
fractures_description  = [];
num_frac = size(FRAC_input,1);
for index = 1:num_frac
        local_frac = FRAC_input(index,:);
        fractures_description = [fractures_description; ...
                           generate_frac_struct(local_frac,Lsize, Nsize)];

end
end

function fracture_geometry_blueprint = generate_blueprint_description(FracturesProperties)
fracture_geometry_blueprint = struct();
fracture_geometry_blueprint.Type = FracturesProperties.blueprint.TYPE;
fracture_geometry_blueprint.Discretization  = convertCharsToStrings(FracturesProperties.blueprint.GRID_DISCRETIZATION);
fracture_geometry_blueprint.NumOfDimensions = str2double(FracturesProperties.blueprint.DOMAIN(1));     
fracture_geometry_blueprint.ADM = FracturesProperties.ReservoirADM;
fracture_geometry_blueprint.Size = [FracturesProperties.Reservoir.Lx  ...
    FracturesProperties.Reservoir.Ly ...
   FracturesProperties.Reservoir.Lz]';
fracture_geometry_blueprint.N = [FracturesProperties.Reservoir.Nx  ...
    FracturesProperties.Reservoir.Ny ...
   FracturesProperties.Reservoir.Nz]';
fracture_geometry_blueprint.compareAccuracy = FracturesProperties.blueprint.CompareAccuracy;
fracture_geometry_blueprint.FracturesProperties = ...
    generate_fracture_description(FracturesProperties.blueprint.FRAC, ...
    fracture_geometry_blueprint.Size, fracture_geometry_blueprint.N);
end



function fracture_struct = generate_frac_struct(frac, Lsize, Nsize)
    LX = Lsize(1);
    LY = Lsize(2);
    LZ = Lsize(3);
    NX = Nsize(1);
    NY = Nsize(2);
    NZS = Nsize(3);
    fracture_struct = struct();
    fracture_struct.isActive = str2double(frac{1});
    fracture_struct.CenterCoord = eval(frac{2})';
    fracture_struct.Length = str2double(frac{3});
    fracture_struct.Width = str2double(frac{4});
    fracture_struct.RotAngleAlongZ = str2double(frac{5});
    fracture_struct.RotAngleAlongL = str2double(frac{6});
    fracture_struct.RotAngleAlongW = str2double(frac{7});
    fracture_struct.GridResRatio = str2double(frac{8});
    fracture_struct.GridNumAlongL = str2double(frac{9});
    fracture_struct.GridNumAlongW = str2double(frac{10});
    fracture_struct.Aperture = str2double(frac{11});
    fracture_struct.Porosity = str2double(frac{12});
    fracture_struct.Permeability = str2double(frac{13});
    adm_expr = join(convertCharsToStrings(frac(:,14:17)),",");
    fracture_struct.ADM = eval(adm_expr);
    corner_expression = join(convertCharsToStrings(frac(:,18:24)),",");
    fracture_struct.CornerCoords = eval(corner_expression);
end


function [isActive, CoarseLevel, CoarseRatio] = split_adm_input(adm_text)
    if numel(adm_text) == 1
        split_text = strsplit(adm_text,{' ',',','[',']'});
        ref = ismember(split_text,'');
        split_text = str2double(split_text(~ref));
        isActive = split_text(1);
        CoarseLevel = split_text(2);
        CoarseRatio = split_text(3:end);
    elseif numel(adm_text) == 3
        isActive = adm_text{1};
        CoarseLevel = adm_text{2};
        CoarseRatio = eval(adm_text{3});
    end
end