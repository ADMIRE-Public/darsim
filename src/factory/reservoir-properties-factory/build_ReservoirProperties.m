%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Artur Castiel;
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ReservoirProperties = build_ReservoirProperties(blueprint)
ReservoirProperties = struct();
ReservoirProperties = build_mesh(ReservoirProperties, blueprint, ...
                                blueprint.get_folder_file());
ReservoirProperties = build_media_information(blueprint, ...
    ReservoirProperties);
end



function ReservoirProperties = build_media_information(blueprint,...
                                                ReservoirProperties) 

ReservoirProperties.Temperature = blueprint.BC.INITIAL_TEMPERATURE;
ReservoirProperties.PermUnit = blueprint.ROCK.PERMEABILITY_UNIT;
ReservoirProperties.PermScale = blueprint.ROCK.PERMEABILITY_SCALE;
perm_contrast_reduction = blueprint.ROCK.PERMEABILITY_CONTRAST_REDUCTION;
if isempty(perm_contrast_reduction)
    ReservoirProperties.PermContrastReduction = true;
    perm_mean = blueprint.ROCK.PERMEABILITY_CONTRAST_MEAN;
    ReservoirProperties.PermContrastMean = perm_mean;
    perm_order = blueprint.ROCK.PERMEABILITY_CONTRAST_ORDER;
    ReservoirProperties.PermContrastOrder = perm_order;
elseif ~perm_contrast_reduction
    ReservoirProperties.PermContrastReduction = false;
end
[ReservoirProperties, flag] = read_from_props_mesh_data(...
        blueprint.GRID.GRID_DISCRETIZATION, ReservoirProperties);
if strcmp(ReservoirProperties.Discretization,'CornerPointGrid') && sum(strcmp(fieldnames(ReservoirProperties.CornerPointGridData), 'Porosity')) && ~isempty(ReservoirProperties.CornerPointGridData.Porosity)
    % Do nothing, beacuse the porosity data already exist in CornerPointGridData.
    ReservoirProperties.phi = ReservoirProperties.CornerPointGridData.Porosity;
    ReservoirProperties.PorosityInclude = 0;
else
ReservoirProperties.phi = blueprint.ROCK.POROSITY;
end
if strcmp(ReservoirProperties.Discretization,'CornerPointGrid') && sum(strcmp(fieldnames(ReservoirProperties.CornerPointGridData), 'Permeability')) && ~isempty(ReservoirProperties.CornerPointGridData.Permeability)
                ReservoirProperties.Perm = ReservoirProperties.CornerPointGridData.Permeability;
else
    %%if Permeability Porosity are not found in the grid generate perm field.
    if ~flag
        rock_perm = {blueprint.ROCK.PERMEABILITY_X,  ...
        blueprint.ROCK.PERMEABILITY_Y, blueprint.ROCK.PERMEABILITY_Z};
        ReservoirProperties.Perm = generate_permeability_field(rock_perm,  ...
        blueprint.get_folder_file, ReservoirProperties);
    end
end
ReservoirProperties.KRotAngle = [blueprint.ROCK.ROTANGLE_X ...
                    blueprint.ROCK.ROTANGLE_Y blueprint.ROCK.ROTANGLE_Z];

ReservoirProperties.Perm;
ReservoirProperties = generate_porfield(ReservoirProperties,blueprint);
if isfield(blueprint.OTHERS, "FACIES_INDEX")
    ReservoirProperties.FaciesIndex = load(fullfile(blueprint.get_folder_file, blueprint.OTHERS.FACIES_INDEX));
end

end

function [ReservoirProperties] = generate_porfield(ReservoirProperties,...
    blueprint) % fixed by Mengjie
if strcmp(ReservoirProperties.Discretization,'CornerPointGrid') && ...
        sum(strcmp(fieldnames(ReservoirProperties.CornerPointGridData), ...
        'Porosity')) && ~isempty(ReservoirProperties.CornerPointGridData.Porosity)
    ReservoirProperties.phi = ...
        ReservoirProperties.CornerPointGridData.Porosity;
    ReservoirProperties.PorosityInclude = false;
else
    por = blueprint.ROCK.POROSITY;
    if  isstring(por) || ischar(por)
        ReservoirProperties.PorosityInclude = true;
        ReservoirProperties.PorosityFile = fullfile(...
            blueprint.get_folder_file, por);
        N_active_cells = ReservoirProperties.Grid.N_ActiveCells;
        grid_size = ReservoirProperties.Grid.N';
        ReservoirProperties.phi = read_unique_input_file_por(ReservoirProperties.PorosityFile, grid_size, ...
        N_active_cells);

    else
        ReservoirProperties.PorosityInclude = false;
        ReservoirProperties.phi = ...
            por .* ones(ReservoirProperties.Grid.N_ActiveCells, 1);        
    end
end
end

function [ReservoirProperties, flag] = ...
    read_from_props_mesh_data(GRID_DISCRETIZATION, ReservoirProperties)
    flag = false;
if strcmp(GRID_DISCRETIZATION, "CornerPointGrid")
    Permeability = ...
        fetch_label_mesh(ReservoirProperties,"Permeability");
    Porosity = ...
        fetch_label_mesh(ReservoirProperties,"Porosity");
    PermUnit = ...
        fetch_label_mesh(ReservoirProperties,"PermUnit");
    PermScale = ...
        fetch_label_mesh(ReservoirProperties,"PermScale");
    if ~isempty(PermScale)
        if ~strcmp(PermUnit, blueprint.ROCK.PERMEABILITY_UNIT)
            warning("Permeability Unit inputfile and CPG file are " + ...
                "different.")
        end
    end
    if ~isempty(Permeability)
        ReservoirProperties.Perm = Permeability;
        flag = true;
    end
    if ~isempty(PermUnit)
        ReservoirProperties.PermUnit = PermUnit;
        flag = true;
    end
    if ~isempty(Porosity)
        ReservoirProperties.Porosity = Porosity;
        flag = true;
    end
    if ~isempty(PermScale)
        ReservoirProperties.Scale = PermScale;
        flag = true;
    end
elseif strcmp(GRID_DISCRETIZATION, "UnstructuredGrid")
    Perm = ...
        fetch_label_mesh(ReservoirProperties,"perm");
    PermUnit = ...
        fetch_label_mesh(ReservoirProperties,"PermUnit");
    if ~isempty(Perm)
        ReservoirProperties.Perm = Perm;
        ReservoirProperties.PermUnit = PermUnit;
        flag = true;
    end
end

end




function prop = fetch_label_mesh(struct_data, label)
prop  = [];
cpg_flag =  strcmp(struct_data.Discretization, "CornerPointGrid");
uns_flag =  strcmp(struct_data.Discretization, "UnstructuredGrid");
if cpg_flag
    if isfield(struct_data.CornerPointGridData, label)
        prop = struct_data.CornerPointGridData.(label);
    end
end
if uns_flag
    if isfield(struct_data.UnstructuredGridData.props, label)
        prop = struct_data.UnstructuredGridData.props.(label);
    end
end
end

function phi = read_unique_input_file_por(unique_input_file, grid_size, n_total)
     raw_data = load(unique_input_file);
     file_grid_info = raw_data(1:3);
     if ~isnan(grid_size)
        grid_flag = isequal(file_grid_info,grid_size');
        uns_flag_header = false;
        uns_flag_headerless = false;
        structured_flag = true;
     else
        uns_flag_header = legnth(raw_data) + 1  == n_total;
        uns_flag_headerless = legnth(raw_data) == n_total;
        grid_flag = uns_flag_header | uns_flag_headerless;
        structured_flag = false;
     end  
     
     if ~grid_flag 
         error("The size of the permeability file is not consistent with" + ...
             "the mesh.")
     end

     if structured_flag
         phi = [raw_data(4:end)];
     elseif uns_flag_header
         phi = [raw_data(2:end)];
     elseif uns_flag_headerless
         phi = [raw_data(1:end)];
     end
end



            % % Read the permeability
            % if strcmp(ReservoirProperties.Discretization,'CornerPointGrid') && sum(strcmp(fieldnames(ReservoirProperties.CornerPointGridData), 'Permeability'))
            %     % Do nothing, beacuse the permeability data already exist in CornerPointGridData.
            %     ReservoirProperties.Perm = ReservoirProperties.CornerPointGridData.Permeability;
            %     if ~strcmp(ReservoirProperties.PermUnit, ReservoirProperties.CornerPointGridData.PermUnit)
            %         warning('DARSim Warning: The permeability unit reported in the main input file and the CornerPointGrid rock properties input file do not match.\nThe one from CornerPointGrid rock properties input file (%s) is considered.\n', ReservoirProperties.CornerPointGridData.PermUnit);
            %     end
            %     if ~strcmp(ReservoirProperties.PermScale, ReservoirProperties.CornerPointGridData.PermScale)
            %         warning('DARSim Warning: The permeability scale reported in the main input file and the CornerPointGrid rock properties input file do not match.\nThe one from CornerPointGrid rock properties input file (%s) is considered.\n', ReservoirProperties.CornerPointGridData.PermScale);
            %     end
            % else
            %     % The permeability is read from the main input file.
            %     ReservoirProperties.Perm = [];
            %     perm = zeros(3, 1); rotangle = zeros(3, 1);
            % 
            %     temp = strfind(obj.InputMatrix, 'PERMEABILITY_X');
            %     perm(1) = find(~cellfun('isempty', temp));
            %     if isempty(perm(1))
            %         error('DARSim Error: The keyword "PERMEABILITY_X" is missing. Please check the input file!\n');
            %     end
            %     temp = strfind(obj.InputMatrix, 'PERMEABILITY_Y');
            %     perm(2) = find(~cellfun('isempty', temp));
            %     if isempty(perm(2))
            %         error('DARSim Error: The keyword "PERMEABILITY_Y" is missing. Please check the input file!\n');
            %     end
            %     temp = strfind(obj.InputMatrix, 'PERMEABILITY_Z');
            %     perm(3) = find(~cellfun('isempty', temp));
            %     if isempty(perm(3))
            %         error('DARSim Error: The keyword "PERMEABILITY_Z" is missing. Please check the input file!\n');
            %     end
            %     temp = strfind(obj.InputMatrix, 'ROTANGLE_X');
            %     rot = find(~cellfun('isempty', temp));
            %     if isempty(rot)
            %        warning('DARSim Error: The keyword "ROTANGLE_X" is missing. Please check the input file!\n');
            %     else
            %        rotangle(1) = find(~cellfun('isempty', temp));
            %     end
            %     temp = strfind(obj.InputMatrix, 'ROTANGLE_Y');
            %     rot = find(~cellfun('isempty', temp));
            %     if isempty(rot)
            %         warning('DARSim Error: The keyword "ROTANGLE_Y" is missing. Please check the input file!\n');
            %     else
            %        rotangle(2) = find(~cellfun('isempty', temp));
            %     end
            %     temp = strfind(obj.InputMatrix, 'ROTANGLE_Z');
            %     rot = find(~cellfun('isempty', temp));
            %     if isempty(rot)
            %         warning('DARSim Error: The keyword "ROTANGLE_Z" is missing. Please check the input file!\n');
            %     else
            %        rotangle(3) = find(~cellfun('isempty', temp));
            %     end
            %     DimChar = {'x' , 'y' , 'z'};
            %     for i=1:3
            %         if contains(obj.InputMatrix(perm(i) - 1), 'INCLUDE')
            %             ReservoirProperties.PermInclude(i) = 1;
            %             if isfile( char(obj.InputMatrix(perm(i)+1)) )
            %                 ReservoirProperties.PermFile{i} = char(obj.InputMatrix(perm(i)+1));
            %             elseif isfile( strcat( obj.Directory, '/', char(obj.InputMatrix(perm(i)+1)) ) )
            %                 ReservoirProperties.PermFile{i} = strcat( obj.Directory, '/', char(obj.InputMatrix(perm(i)+1)) );
            %             else
            %                 error('DARSim Error: The path to permeability input file is not valid. Please specify it correctly in the main input file.');
            %             end
            %             ReservoirProperties.KRotAngleInclude(i) = 0;
            %             ReservoirProperties.KRotAngle(i) = 0;
            %         else
            %             ReservoirProperties.PermInclude(i) = 0;
            %             ReservoirProperties.KRotAngleInclude(i) = 0;
            %             ReservoirProperties.Perm(i) = str2double(obj.InputMatrix(perm(i) +1));
            %             ReservoirProperties.KRotAngle(i) = str2double(obj.InputMatrix(rotangle(i) +1));
            %             if isnan(ReservoirProperties.Perm(i))
            %                 error('DARSim Error: The permeability value in %s direction is not valid. Please check the input file.\n',DimChar{i});
            %             end
            %             if isnan(ReservoirProperties.KRotAngle(i))
            %                 warning('DARSim Error: The permeability rotation angle value in %s direction is not valid. Please check the input file.\n',DimChar{i});
            %                 ReservoirProperties.KRotAngle(i) = 0;
            %             end
            %         end
            %     end
            % 
            %     %%%% Homogenized/upscaled permeabilities for each coarsening level %%%%%%%
            %     temp = strfind(obj.SettingsMatrix, 'DLGR');
            %     x = find(~cellfun('isempty', temp));
            %     if isempty(temp)
            %         DLGR = 0;
            %     else
            %         DLGR = str2double(obj.SettingsMatrix(x+1));
            %     end
            %     if DLGR
            %         temp = strfind(obj.SettingsMatrix, 'ADM');
            %         adm = find(~cellfun('isempty', temp));
            %         if str2double(obj.SettingsMatrix(adm + 1)) == 1
            %             temp = strfind(obj.SettingsMatrix, 'LEVELS');
            %             x = find(~cellfun('isempty', temp));
            %             CoarseningLevel = str2double(obj.SettingsMatrix(x+1));
            %             for L = 1:CoarseningLevel
            %                 temp = strfind( obj.InputMatrix , strcat('PERM_COARSE_L',num2str(L),'_X') );
            %                 permx = find(~cellfun('isempty', temp));
            %                 temp = strfind( obj.InputMatrix , strcat('PERM_COARSE_L',num2str(L),'_Y') );
            %                 permy = find(~cellfun('isempty', temp));
            %                 temp = strfind( obj.InputMatrix , strcat('PERM_COARSE_L',num2str(L),'_Z') );
            %                 permz = find(~cellfun('isempty', temp));
            %                 if ~isempty(permx)
            %                     ReservoirProperties.CoarsePermFile{L,1} = strcat(obj.PermDirectory, char(obj.InputMatrix(permx+1)));
            %                 end
            %                 if ~isempty(permy)
            %                     ReservoirProperties.CoarsePermFile{L,2} = strcat(obj.PermDirectory, char(obj.InputMatrix(permy+1)));
            %                 end
            %                 if ~isempty(permz)
            %                     ReservoirProperties.CoarsePermFile{L,3} = strcat(obj.PermDirectory, char(obj.InputMatrix(permz+1)));
            %                 end
            %             end
            %         end
            %     end
            % end
