function [K] = generate_permeability_field(var, folder,...
                                        ReservoirProperties)
perm_unit_mult = dictionary('m2',1,'D', 1e-12,'mD',1e-15);
perm_scale = dictionary('Linear',1, 'Logarithmic', 0);
perm_scale_flag = perm_scale(ReservoirProperties.PermScale);
perm_unit_constant = perm_unit_mult(ReservoirProperties.PermUnit);
[flag] = input_type_evaluate(var);
flag_unique_input_file = strcmp(var{1}, var{2}) && strcmp(var{1}, var{3}); %all(flag); Fixed by Amin, Case of reading permeabilities from the same file
flag_numerical_input_file = all(~flag);
N_active_cells = ReservoirProperties.Grid.N_ActiveCells;
grid_size = ReservoirProperties.Grid.N';
if flag_unique_input_file
    unique_input_file = fullfile(folder,findIdenticalValue(var));
    K = read_unique_input_file(unique_input_file, grid_size, ...
        N_active_cells);
    K = calculate_perm(perm_scale_flag, perm_unit_constant, K);
elseif flag_numerical_input_file
    K = ones(N_active_cells, 3);
    K(:,1) = K(:,1) * calculate_perm(perm_scale_flag, ...
        perm_unit_constant, var{1});
    K(:,2) = K(:,2) * calculate_perm(perm_scale_flag, ...
        perm_unit_constant, var{2});
    K(:,3) = K(:,3) * calculate_perm(perm_scale_flag, ...
        perm_unit_constant, var{3});  
else
    input_files = createFullPaths(var, folder);
    K = read_from_different_files(input_files,grid_size,N_active_cells, ...
        flag, var, perm_scale_flag, perm_unit_constant);
end
if ReservoirProperties.PermContrastReduction && ~flag_numerical_input_file
    K =  apply_perm_contrast_filter(K, ...
        ReservoirProperties.PermContrastMean, ...
        ReservoirProperties.PermContrastOrder, perm_unit_constant);
end
end

function K_filtered =  apply_perm_contrast_filter(K, k_mean, k_order, ...
    perm_multiplier)
warning("needs testing")
K_Log10 = log10(K);
if isstring(k_mean) || ischar(k_mean)
    if strcmpi(k_mean, "default")
        k_mean = mean(K,1);
    end
else
    k_mean = k_mean*perm_multiplier;
end
k_mean = log10(k_mean);
k_mean = repmat(k_mean,size(K,1),1);
maxContrast = k_order;
ratio = (max(K_Log10) - min(K_Log10)) / maxContrast;
ratio = repmat(ratio,size(K,1),1);
if ratio ~= 0
    K_filtered = ((K_Log10 - k_mean)./ratio)  +  k_mean;
    K_filtered = 10.^(K_filtered);
else
    K_filtered  = K;
end
end

function K = read_from_different_files(files,grid_size,n_total, flag,...
                                var, perm_scale_flag, multiplier)
    flag_numerical = ~flag;
    K = ones(n_total, 3);
    for index = find(flag_numerical)
        K(:,index)  = K(:,index) * calculate_perm(perm_scale_flag, ...
                        multiplier, var{index});
    end    
    for index = find(flag)
        local_var = read_single_file(files{index}, n_total); %Fixed by Amin; grid_size(index) > n_total
        K(:,index) = calculate_perm(perm_scale_flag, ...
                        multiplier, local_var); %Fixed by Amin; K(:,index) * calculate_perm...
    end    
end

function k_axis = read_single_file(local_file, n_total)
    raw_data = load(local_file);
    flag_headerless = length(raw_data) == n_total;
    flag_header = length(raw_data) == n_total + 1; %Fixed by Amin; (length(raw_data) +1) == n_total
    flag_longheader = length(raw_data) == n_total + 3;
    if flag_headerless
        k_axis = raw_data(1:end);
    elseif flag_header
        k_axis = raw_data(2:end);
    elseif flag_longheader
        k_axis = raw_data(4:end);
    else
        error("The size of the dimension of the permeability file is" + ...
            " not consistent with " + "the mesh.")
    end
end


function K = read_unique_input_file(unique_input_file, grid_size, n_total)
     raw_data = load(unique_input_file);
     file_grid_info = raw_data(1:3);
     if ~isnan(grid_size)
        grid_flag = isequal(file_grid_info,grid_size');
        uns_flag_header = false;
        uns_flag_headerless = false;
        structured_flag = true;
     else
        uns_flag_header = length(raw_data) + 1  == n_total;
        uns_flag_headerless = length(raw_data) == n_total;
        grid_flag = uns_flag_header | uns_flag_headerless;
        structured_flag = false;
     end  
     
     if ~grid_flag 
         error("The size of the permeability file is not consistent with" + ...
             "the mesh.")
     end

     if structured_flag
         K = [raw_data(4:end) raw_data(4:end) raw_data(4:end)];
     elseif uns_flag_header
         K = [raw_data(2:end) raw_data(2:end) raw_data(2:end)];
     elseif uns_flag_headerless
         K = [raw_data(1:end) raw_data(1:end) raw_data(1:end)];
     end

end

function perm = calculate_perm(flag_type, multiplier, var)
    if flag_type
        perm = multiplier.*var;
    else
        perm = multiplier.*(10.^var);
    end
end

function fullPathCellArray = createFullPaths(fileNames, folder)
    % Check if the input cell array is empty
    if isempty(fileNames) || isempty(folder)
        error('Input cell array or folder path is empty.');
    end
    
    % Initialize the cell array for full paths
    fullPathCellArray = cell(size(fileNames));    
    % Iterate through the file names and create full paths
    for i = 1:numel(fileNames)
        fullPathCellArray{i} = fullfile(folder, fileNames{i});
    end
end

function identicalValue = findIdenticalValue(cellArray)
    % Check if the input cell array is empty
    if isempty(cellArray)
        identicalValue = cellArray;
        return;
    end
    
    % Check if the cell array has at least two elements to compare
    if numel(cellArray) < 2
        identicalValue = cellArray;
        return;
    end
    
    % Iterate through the elements of the cell array
    for i = 1:numel(cellArray)-1
        for j = i+1:numel(cellArray)
            % Compare the current pair of elements for equality
            if isequal(cellArray{i}, cellArray{j})
                identicalValue = cellArray{i};
                return;
            end
        end
    end
    
    % If no identical values are found, return the original cell array
    identicalValue = cellArray;
end