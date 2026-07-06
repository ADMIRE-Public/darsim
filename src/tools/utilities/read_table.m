function Table = read_table(filename)

    fid = fopen(filename, 'r');
    if fid == -1
        error('Could not open file: %s', filename);
    end

    Table = struct();
    currentName = '';
    currentData = [];

    while ~feof(fid)
        line = strtrim(fgetl(fid));

        % Skip empty lines
        if isempty(line)
            continue
        end

        % Detect section header, e.g. --wetting_kr or ----Pc_drainage
        if startsWith(line, '-')
            % Save previous section
            if ~isempty(currentName)
                Table.(currentName) = currentData;
            end

            % Remove leading dashes and make a valid MATLAB field name
            name = regexprep(line, '^-+', '');
            currentName = matlab.lang.makeValidName(strtrim(name));

            % Start new data block
            currentData = [];

        else
            % Read numeric data line
            nums = sscanf(line, '%f');

            if numel(nums) == 2
                currentData(end+1, :) = nums(:).';
            end
        end
    end

    % Save last section
    if ~isempty(currentName)
        Table.(currentName) = currentData;
    end

    fclose(fid);
end
% function [table] = read_table(file)
% table = struct();
% file_cell =  read_to_cell(file);
% file_cell = file_cell(~ismember(file_cell,''));
% cleanedText = removeHyphenatedText(file_cell);
% keywords = ["wetting_kr", "nonwetting_drainage_kr", ...
%     "nonwetting_imbibation_kr", "Pc_drainage", "Pc_imbibation", "end"];
% ref_keywords = find(ismember(cleanedText,keywords));
% for word = keywords(1:end-1)
%     ref_start = ref_keywords(ismember(keywords,word));    
%     ref_end = ref_keywords(ref_keywords > ref_start);
%     ref_end = ref_end(1);
%     table.(word) = cellArrayToMatrix(cleanedText(ref_start+1:ref_end-1));
% end
% end
% 
% 
% function outputMatrix = cellArrayToMatrix(inputCellArray)
%     % Get the number of rows in the cell array
%     numRows = numel(inputCellArray);
% 
%     % Initialize the output matrix with zeros
%     outputMatrix = zeros(numRows, 2);
% 
%     % Loop through each cell in the input cell array
%     for i = 1:numRows
%         % Split the string using both tab and space as delimiters
%         parts = strsplit(inputCellArray{i}, {'\t', ' '});
% 
%         % Convert the string parts to numeric values
%         if numel(parts) == 2
%             value1 = str2double(parts{1});
%             value2 = str2double(parts{2});
% 
%             % Check if the conversion was successful
%             if ~isnan(value1) && ~isnan(value2)
%                 outputMatrix(i, 1) = value1;
%                 outputMatrix(i, 2) = value2;
%             else
%                 error('Error: Invalid numeric values in cell %d', i);
%             end
%         else
%             error('Error: Invalid format in cell %d', i);
%         end
%     end
% end
% 
% function cleanedCellArray = removeHyphenatedText(cellArray)
%     % Define the regular expression pattern to match hyphenated text with names
%     pattern = '(-+)+(\w+)';
% 
%     % Initialize a cell array to store the cleaned text
%     cleanedCellArray = cell(size(cellArray));
% 
%     % Loop through each element of the input cell array
%     for i = 1:numel(cellArray)
%         % Find all matches of the pattern in the current cell
%         matches = regexp(cellArray{i}, pattern, 'tokens');
% 
%         % Initialize a copy of the current cell's text
%         cleanedCellArray{i} = cellArray{i};
% 
%         % Process matches and update the cleaned cell array
%         for j = 1:numel(matches)
%             % Extract the name from the match
%             name = matches{j}{2};
% 
%             % Replace the match with the extracted name in the cleaned cell
%             cleanedCellArray{i} = name;
% 
%             % Display the removed text (optional)
%         end
%     end
% end