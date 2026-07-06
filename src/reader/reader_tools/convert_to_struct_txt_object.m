function [blueprint] = convert_to_struct_txt_object(InputMatrix)
% verify and checks if all required keywords are available. converts the
% object to a matlab structure;
blueprint = struct();
for label = InputMatrix.mapped_data_label'
    value = replaceKeywordsWithValues(InputMatrix.read_mapped_data(label));
    value = extractNumerics(value);
    blueprint.(label)  = value;
end

for index = 1:length(InputMatrix.range_data_label')
    label = InputMatrix.range_data_label(index);
    count = InputMatrix.range_data_count(index);
    value = cell(count,1);
    for ii  = 1:count
        local_value = convert_to_struct_txt_object((...
                            InputMatrix.read_range_data(label,ii)));
        value{ii} = local_value;
    end
    blueprint.(label)  = value;
end
if InputMatrix.isfrac
    blueprint.("FRAC") = InputMatrix.read_frac_range_data;
end
blueprint = convert_string_to_matlab_cells_matrix(blueprint);
end

function blueprint = convert_string_to_matlab_cells_matrix(blueprint)
fields = fieldnames(blueprint);
for index = 1:length(fields)
    field = fields{index};
    value = blueprint.(field);
    if isstring(value) || ischar(value)
        if isEnclosedByBrackets(value)
            value = eval(value);
            blueprint.(field) = value;
        end
    end
end
end

function result = isEnclosedByBrackets(inputString)
% Check if the input string is enclosed by {} or []

% Remove leading and trailing whitespace from the input string
inputString = strtrim(inputString);

if isempty(inputString)
    result = 0; % An empty string is not enclosed by brackets
else
    % Check if the first and last characters are { } or [ ]
    firstChar = inputString(1);
    lastChar = inputString(end);

    if (firstChar == '{' && lastChar == '}') ...
            || (firstChar == '[' && lastChar == ']')
        result = 1; % Enclosed by brackets
    else
        result = 0; % Not enclosed by brackets
    end
end
end

function modifiedArray = replaceKeywordsWithValues(inputArray)
% Initialize the output cell array
keywords = ["ON", "TRUE", "OFF", "FALSE", "NAN"];
keyValues = [true, true, false, false, nan];
modifiedArray = inputArray;
% Ensure that keywords and keyValues have the same length
if numel(keywords) ~= numel(keyValues)
    error('Number of keywords must be equal to the number of key values.');
end

% Loop through the elements of the input cell array
for i = 1:numel(inputArray)
    % Find the index of the current element in the keywords array
    keywordIndex = find(strcmpi(inputArray{i}, keywords));
    % If a matching keyword is found, replace with the corresponding value
    if ~isempty(keywordIndex)
        modifiedArray{i} = keyValues(keywordIndex);
    else
        value = convertStringToNumber(inputArray{i});
        modifiedArray{i} = splitAndConvertToNumbers(value);
    end
end
modifiedArray = modifiedArray';
if length(modifiedArray) == 1
    modifiedArray = modifiedArray{1};
end
end



function result = convertStringToNumber(str)
% Try to convert the string to a numerical value
numValue = str2double(str);
% Check if the conversion was successful and the result is not NaN
if ~isnan(numValue)
    result = numValue;
else
    % Return the original string value if not numerical
    result = str;
end
end


function numbers = splitAndConvertToNumbers(inputStr)
% Split the input string by spaces
if isstring(inputStr)
    splitValues = strsplit(inputStr);

    % Initialize an array to store the converted numbers
    numbers = [];

    % Loop through the split values
    for i = 1:numel(splitValues)
        % Convert each split value to a number
        numValue = str2double(splitValues{i});

        % Check if the conversion was successful and the result is a number
        if ~isnan(numValue)
            numbers = [numbers, numValue];
        else
            % If any value is not a number, return an empty array
            numbers = [];
            return;
        end
    end
else
    numbers = inputStr;
end
end
