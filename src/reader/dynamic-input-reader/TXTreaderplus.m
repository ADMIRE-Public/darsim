%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DARSim
% Author: Artur Castiel Reis de Souza
% TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This class reads raw text files in the DARSim input file formats automatically.
% New labels and keys can be automatically added directlly on the Input files without altering
% this routine. 
% raw input file is store based on the label, and size of keys.  Depending on the size some labels are 
% flaged as mapped or ranged data.
% The classs tries to consider all previous DARSim's input file formats. 
classdef TXTreaderplus
    properties
        mapped_data_label
        range_data_label
        range_data_count
        isfrac
    end

    properties(Access=private)
        range_data
        raw_input
        search_vector
        keywords = ["ON", 'OFF', 'INJ', 'PROD', 'XYZ_LIST', ...
            'VTK', 'COUPLED','BINARY', 'ASCII', 'Multiscale', ...
            'Constant', "2D", "3D", "1D", "EDFM", "pEDFM", ...
            "CartesianGrid", "CornerPointGrid", "UnstructuredGrid",...
            "INCLUDE", "NaN", "INCLUDE_ROCK_PROPERTIES", ...
            "CO2-Brine", "H2-Brine", "H2-CH4-Brine", "srw_frac", "Table"];
        atypical_fields = ["pEDFM_", "Fix"];
        negative_keywords = ["INCLUDE"];
        fracture_keywords = ["FRAC_START" "FRAC_END"];
        blank_maps = [""];
    end

    methods
        function obj = TXTreaderplus(file)
            if isstring(file) || ischar(file)
                raw_input = generate_raw_input_data(file);
            else
                raw_input = file;
            end
            raw_input = remove_negative_keyword(strtrim(...
                removeNestedParentheses(removeDoubleDashCells(...
                remove_coments(raw_input)))),obj.negative_keywords);

           
            [obj.mapped_data_label, obj.range_data_label,...
                obj.range_data, obj.raw_input, obj.search_vector,...
                obj.range_data_count, obj.isfrac] = ...
                                        create_map_label(obj,raw_input);
            
        
        end
        function [mapped_data] = read_mapped_data(obj, field)
            if any(ismember(obj.mapped_data_label,field))
                ref_pos = ismember(obj.raw_input,field);
                ref_value = ...
                    findZerosBetweenOnes(obj.search_vector, find(ref_pos));

                mapped_data = obj.raw_input(ref_value);
                mapped_data = convertCharsToStrings(mapped_data);
                if isstring(mapped_data)
                    mapped_data = arrayfun(@(x)char(mapped_data(x)), ...
                        1:numel(mapped_data),'uni',false);
                end
                isNumeric = checkCellArray(mapped_data);
                mapped_data = convert_to_str(mapped_data,isNumeric);
            else
                mapped_data = NaN;
            end
        end
        function [local_range_data] = read_range_data(obj,field, num)
            if ~strcmp("FRAC",field)
                local_range_data = obj.range_data.(field){num};
            else
                local_range_data = obj.range_data.(field);
            end
        end
        function [local_range_data] = read_frac_range_data(obj)
            if obj.isfrac     
                local_range_data = obj.range_data.("FRAC");
            else
                error("File does not contain FRAC_RANGE");
            end
        end

    end
    methods (Access = private)
        function [mapped_data_label, range_data_label, range_data, ...
                raw_input, search_vector, counts, is_frac] = ...
                create_map_label(obj,raw_input)
            [ref, raw_input] = keyword_pattern_filtering(raw_input, ...
                obj.keywords, obj.atypical_fields,obj.fracture_keywords);
            map_keys = raw_input(ref);
            ref_start = contains(map_keys,"START");
            range_data_label = replace(convertCharsToStrings ...
                (convertCharsToStrings(map_keys(ref_start))),"_START","");
            range_data_label = removeNestedParentheses(range_data_label);
            [range_data_label, ~, idx] = unique(range_data_label);
            counts = histcounts(idx, 1:length(range_data_label)+1);           
            [range_data, raw_input] = ...
                split_range_data(obj,range_data_label,counts,raw_input);
         

            %% filtrar mapped data before returning
            ref = keyword_pattern_filtering(raw_input, obj.keywords, ...
                obj.atypical_fields, obj.fracture_keywords);
            [search_vector] = create_search_vector(ref,raw_input, ...
                                 obj.blank_maps);
            map_keys = raw_input(search_vector);
            mapped_data_label = convertCharsToStrings(map_keys);
            frac_ref = ismember(range_data_label,"FRAC");
            is_frac = any(frac_ref);
            range_data_label(frac_ref) = [];
            counts(frac_ref) = [];
        end

        function [range_data,raw_input] = ...
                split_range_data(~,labels,counts,raw_input)
            range_data = struct("fields",labels);
            for index = 1:length(labels)
                label = labels(index);
                local_counts = counts(index);
                local_data = cell(local_counts,1);
                if ~strcmp(label,"FRAC")
                    for ii = 1:local_counts
                        ref_start = ...
                            contains(raw_input,join([label, "_START"],""));
                        ref_start = find(ref_start);
                        ref_start = ref_start(1);
                        ref_end = ...
                            contains(raw_input,join([label, "_END"],""));
                        ref_end = find(ref_end);
                        ref_end = ref_end(1);
                        local_data{ii} = TXTreaderplus(raw_input( ...
                            ref_start+1:ref_end-1));
                        raw_input(ref_start:ref_end) = [];

                    end
                else
                    ref_start = ...
                        contains(raw_input,join([label, "_START"],""));
                    ref_start = find(ref_start);
                    ref_start = ref_start(1);
                    ref_end = ...
                        contains(raw_input,join([label, "_END"],""));
                    ref_end = find(ref_end);
                    ref_end = ref_end(1);
                    data = raw_input(ref_start+1:ref_end-1);
                    data = cellfun(@(x) strsplit(x, ","), data, ...
                        'UniformOutput', false);                   
                    local_data = vertcat(data{:});
                    raw_input(ref_start:ref_end) = [];
                end
                range_data.(label) = local_data;
            end
        end
    end
end

function clean_raw = remove_negative_keyword(raw_input, bad_keys)
    ref = ismember(raw_input,bad_keys);
    clean_raw = raw_input(~ref);
end

function raw_input = generate_raw_input_data(file)
fileID = fopen(file, 'r');
%// Read lines from input file
matrix = textscan(fileID, '%s', 'Delimiter', '\n');
raw_input = strtrim(matrix{1});
fclose(fileID);
end

function isNumeric = checkCellArray(cellArray)
    % Initialize a logical array to store the results
    isNumeric = false(size(cellArray));

    % Loop through each cell in the cell array
    for i = 1:numel(cellArray)
        % Check if the cell contains only text (no digits)
        if isempty(regexp(cellArray{i}, '\d', 'once'))
            isNumeric(i) = false;  % Cell contains only characters
        else
            % Try to convert the text to a numeric value
            numericValue = str2double(cellArray{i});
            if ~isnan(numericValue)
                isNumeric(i) = true;  % Cell contains a numeric value
            else
                isNumeric(i) = false; % Cell contains non-numeric characters
            end
        end
    end

end

function outputString = remove_coments(inputString)
    % Define a regular expression pattern to match text within parentheses
    pattern = '\([^()]*\)|(--|//)[^\n]*';

    % Replace all occurrences of the pattern with an empty string
    while true
        updatedString = regexprep(inputString, pattern, '');
        if isequal(updatedString, inputString)
            break; % No more matches found, exit the loop
        end
        inputString = updatedString;
    end
    
    outputString = inputString;
end


function outputString = removeNestedParentheses(inputString)
    % Define a regular expression pattern to match text within parentheses
    inputString = strtrim(inputString);
    pattern = '\([^()]*\)';
    % Replace all occurrences of the pattern with an empty string
    while true
        updatedString = regexprep(inputString, pattern, '');
        if isequal(updatedString, inputString)
            break; % No more matches found, exit the loop
        end
        inputString = updatedString;
    end    
    outputString = strtrim(inputString);
    %cleanedStrings = removeUnits(inputString);
    %outputString = cleanedStrings;

end

function filteredCellArray = removeDoubleDashCells(cellArray)
    % Define a regular expression pattern to match cells starting with '--'
    pattern = '^--.*';

    % Use regexprep to remove cells matching the pattern
    filteredCellArray = regexprep(cellArray, pattern, '');

    % Remove any empty cells from the result
    filteredCellArray = filteredCellArray(~cellfun('isempty', ...
                    filteredCellArray));
end

function cell_array = convert_to_str(cell_array, range)
for index = find(range)
    cell_array{index} = str2double(cell_array{index});
end
end


function [search_vector] = create_search_vector(ref, raw_data, blank_map)
ref_blankmap_candidates = ismember(raw_data,blank_map);
ref_blank = false(size(ref));
for index = find(ref_blankmap_candidates)'
    if ref(index+1) == true
        ref_blank(index) = true;
    end
end
ref(ref_blank) = 0;
search_vector = [ref;false];
ref_list = find(ref)';
max_ref = max(ref_list);

for index = ref_list
    if search_vector(index) && search_vector(index+1) && ((index) < max_ref)
        [part1, ~] = splitPattern(raw_data{index});
        [part3, ~] = splitPattern(raw_data{index+1});
        if ~strcmp(part1,part3)
            rf_end= find(search_vector(index+1:end) == 0,1) - 2;
            search_vector(index+1:(index+1+rf_end)) = false;
        end
    else     
    end
end
search_vector(end) = [];
search_vector(end) = false;
search_vector(ref_blank) = true;
end


function [firstPart, secondPart] = splitPattern(inputString)
    % Find the first underscore occurrence
    underscoreIndex = strfind(inputString, '_');
    if ~isempty(underscoreIndex)
        % Extract the first part before the first underscore
        firstPart = inputString(1:underscoreIndex(1)-1);

        % Extract the second part after the first underscore
        secondPart = inputString(underscoreIndex(1)+1:end);
    else
        % Set parts to the input string if no underscore is found
        firstPart = inputString;
        secondPart = '';
    end
end

function zeroPositions = findZerosBetweenOnes(A, POS)
    % Find the next 1 position after POS
    A(end+1) = true;
    nextOneIndex = find(A(POS+1:end) == 1, 1);
    if ~isempty(nextOneIndex)
        % Find the positions of 0s between the current 1 and the next 1
        zeroPositions = find(A(POS+1:POS+nextOneIndex-1) == 0) + POS;
    else
        % If there are no more 1s, return an empty result
        zeroPositions = [];
    end
end

function [ref, raw_input] = keyword_pattern_filtering(raw_input, ...
                                keywords, atypical_fields, fracture_keys)
pat = regexpPattern("^[A-Z].*");
ref = contains(raw_input,pat);

ref_title = contains(raw_input,"TITLE");
ref(find(ref_title) +1) = false;
ref_cammel = hasCamelCase(raw_input);
ref(ref_cammel) = false;
ref_keywords = ismember(raw_input,keywords);
ref(ref_keywords) = false;
ref(contains(raw_input, ".")) = 0;
ref_positive = contains(raw_input,atypical_fields);
ref_var = find(ismember(raw_input,"VARIABLE"))+1;
%% take a look here
sumref = find(ref(ref_var:end) == 0,1) - 3;
ref(ref_var:(ref_var+sumref)) = 0;
ref(ref_positive) = true;

ref_frac = false(size(ref));
[ind_frac_ind] = find(contains(raw_input, fracture_keys));
if ~isempty(ind_frac_ind)
    ref_frac(ind_frac_ind(1)+1:ind_frac_ind(2)-1) = true;
end
raw_input(~ref_frac) = removeUnits(raw_input(~ref_frac));
end



function camelCaseIndicator = hasCamelCase(strings)
    % Initialize a logical array to store the indicators
    camelCaseIndicator = false(size(strings));
    % Loop through the strings
    for i = 1:numel(strings)
        % Use regexp to check if the current string is in CamelCase
        matches = regexp(strings{i}, ...
            '[A-Z][a-z]*(_[A-Z][a-z]*)*', 'match');        
        % Check if there are matches and if the 
        % first match covers the entire string
        if ~isempty(matches) && strcmp(matches{1}, strings{i})
            camelCaseIndicator(i) = true;
        end
    end
end

function cleanedStrings = removeUnits(strings)
    ref1 = contains(strings,"[");
    ref2 = contains(strings,"]");
    ref  = ref1 & ref2;
    refS1 = startsWith(strings,"[");
    refS2 = endsWith(strings,"]");
    refS = refS1 & refS2;   
    ref(refS) = false;
    for index = find(ref)'
        local_string = strings{index};
        [no_comment_string, ~] = splitString(local_string);
        strings{index} = no_comment_string;
    end
    cleanedStrings = strings;
end

function [part1, part2] = splitString(inputString)
    % Split the input string based on white space
    parts = strsplit(inputString);
    % Check if there are exactly two parts
    if numel(parts) == 2
        part1 = parts{1};
        part2 = parts{2};
    else
        part1 = '';
        part2 = '';
    end
end
