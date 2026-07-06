function resultCell = cellArrayToMixed(inputCell)
    % Initialize a cell array to store the results
    resultCell = cell(size(inputCell));

    % Iterate through the input cell array
    for i = 1:numel(inputCell)
        % Use str2double to attempt conversion
        convertedValue = str2double(inputCell{i});
        
        % Check if the conversion is successful (not NaN)
        if ~isnan(convertedValue)
            resultCell{i} = convertedValue; % Store as a number
        else
            resultCell{i} = inputCell{i}; % Store as a string
        end
    end
end