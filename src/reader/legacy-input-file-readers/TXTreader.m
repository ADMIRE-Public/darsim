classdef TXTreader
    properties
        raw_input
        mapped_data_label
        range_data_label
    end

    methods
        function obj = TXTreader(file)
            obj.raw_input = generate_raw_input_data(file);
            [obj.mapped_data_label, obj.range_data_label] = ...
                                                create_map_label(obj);

        end
        function [mapped_data_label, range_data_label] = ...
                                                    create_map_label(obj)
                pat = regexpPattern("^[A-Z].*");
                ref = contains(obj.raw_input,pat);
                qref = diff(ref);
                qref = [false; qref];
                qref = qref == 1;
                nref = qref & qref;
                map_keys = obj.raw_input(nref);
                ref_start = contains(map_keys,"START");
                ref_end = contains(map_keys,"END");
                ref_map_label = ~ (ref_start | ref_end);
                mapped_data_label = ...
                            convertCharsToStrings(map_keys(ref_map_label));
                range_data_label = replace(convertCharsToStrings ...
                 (convertCharsToStrings(map_keys(ref_start))),"_START","");
        end

        function [mapped_data] = read_mapped_data(obj, field)
            mapped_data = find_section(obj.raw_input, field);
            mapped_data = convertCharsToStrings(mapped_data);
            map_conv = regexp(mapped_data ,'^[\-\,\.e0-9]*$', 'once');
            if size(mapped_data,1) > 1
                if ~isempty(cell2mat(map_conv))
                    mapped_data = str2num(cell2mat(mapped_data));
                end
            else
                if ~isempty(map_conv)
                    mapped_data = str2num(mapped_data);
                end
            end
        end
        
        function [range_data] = read_range_data(obj,field)
            start_field = strcat(field,"_START");
            end_field = strcat(field,"_END");
            range_data = find_range(obj.raw_input, start_field, end_field);
        end
    end
end

function data = find_range(file, start_range, end_range)
% returns the data in file of of the corresponding section
temp = strfind(file, start_range);
index_begin = find(~cellfun('isempty', temp));
temp = strfind(file, end_range);
index_end = find(~cellfun('isempty', temp));
if isempty(index_begin)
    error("ERROR: Start Field %s not found.", start_range)
end
if isempty(index_end)
    error("ERROR: End Field %s not found.", end_range)
end
range = index_begin+1:index_end-1;
data = file(range);
data = cellfun(@(x) strsplit(x, ","), data, 'UniformOutput', false);
for index = 1:size(data,1)
    tmp = data{index};
    data(index,1:size(tmp,2)) = tmp;
end
end


function data = find_section(file, section)
% returns the data in file of of the corresponding section
temp = strfind(file, section);
index = find(~cellfun('isempty', temp));
if isempty(index)
    error("ERROR: Field %s not found.", section)
end
index = index+1;
index_begin = index;
while ~startsWith( file{index}, "--")
    index = index + 1;
end
data = file(index_begin:index-2);
if size(data,1) == 1
    data = data{1};
end
end

function raw_input = generate_raw_input_data(file)
fileID = fopen(file, 'r');
%// Read lines from input file
matrix = textscan(fileID, '%s', 'Delimiter', '\n');
raw_input = matrix{1};
fclose(fileID);
end