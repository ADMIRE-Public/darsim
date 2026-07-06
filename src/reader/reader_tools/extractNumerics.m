function var = extractNumerics(var)
size_var = length(var);
if iscell(var)
    for index = 1:size_var
        if isstring(var{index})
            var{index} = extractNumerics_el(var{index});
        end
    end
    if isCellArrayNumerical(var)
        var = cell2mat(var);
    end
elseif ischar(var) || isstring(var)
   var = extractNumerics_el(var);
end
end

function canBeConverted = isCellArrayNumerical(cellArray)
    % Initialize a flag to indicate if the cell array can be converted
    % to a matrix
    canBeConverted = true;

    % Loop through each element of the cell array
    for i = 1:numel(cellArray)
        % Check if the current element is not numeric
        if ~isnumeric(cellArray{i})
            % If any element is not numeric, set the flag to false and break
            canBeConverted = false;
            break;
        end
    end
end

function numericValues = extractNumerics_el(inputStr)
    % Regular expression to match numerical values (including decimals)
    % pattern = '\d+(\.\d+)?';
    %pattern = '[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?';
    special_words = ["m2" , "D", "mD", "2D", "3D"];    
    pattern = '[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?';

    % Find all matches in the input string
    matches = regexp(inputStr, pattern, 'match');    
    flag_not_numerical = isnan(str2double(inputStr));
    % Convert the matched strings to numeric values
    numericValues = str2double(matches);    
    % Check if there are any non-NaN values in the result
    if ismember(inputStr,special_words) || flag_not_numerical 
        numericValues = inputStr;
    else
        % No numeric values found, return the original string
        numericValues = numericValues(~isnan(numericValues));
    end        
    % 
    % 
    % 
    % if any(~isnan(numericValues)) && ~ismember(inputStr,special_words) ...
    %         && ~contains_file_extension(inputStr)
    %     % Remove NaN values (non-numeric parts)
    %     numericValues = numericValues(~isnan(numericValues));
    % else
    %     % No numeric values found, return the original string
    %     numericValues = inputStr;
    % end
end


function isFileExtension = contains_file_extension(str)
    % Define a regular expression pattern to match file extensions
    pattern = '^\.[a-zA-Z]+$'; % Alphabetic characters (upper and lower case) are allowed    
    % Use the regular expression function to check if the input string matches the pattern
    isMatch = regexp(str, pattern, 'once');    
    % If there is a match, isMatch will be a non-empty array
    % If there is a match, it's a valid file extension
    % Otherwise, it's not a valid file extension
    if ~isempty(isMatch)
        isFileExtension = true;
    else
        isFileExtension = false;
    end
end
