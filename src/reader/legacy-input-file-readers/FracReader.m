classdef FracReader
    properties
        blueprint
        folder
        file
        ext
    end
    properties (Access = private)
        fracture_input_raw
    end
    methods
        function obj = FracReader(file)
            [directory, filename, ext] = fileparts(file);
            obj.folder = directory;
            obj.file = filename;
            obj.ext = ext;
            obj.blueprint = struct();
            switch ext
                case ".txt"
                    blueprint_sections = ["Type", "TYPE"; ...
                        "Discretization", "GRID_DISCRETIZATION"; ...
                        "NumOfDimensions", "DOMAIN"; ...
                        "ADM", "ADM_RESERVOIR"; ...
                        "Size", "DIMENS"; ...
                        "N", "SPECGRID"; ...
                        "compareAccuracy", "CompareAccuracy"];
                    blueprint_range = ["FracturesProperties" "FRAC"];
                    obj.fracture_input_raw = TXTreader(file);
                otherwise
                    error("ERROR: Fracture Input File format not supported")
            end
            obj.blueprint = ...
                obj.generate_fracture_info(blueprint_sections);
            dim = convertStringsToChars(obj.blueprint.NumOfDimensions);
            obj.blueprint.NumOfDimensions = str2double(dim(1));
            obj.blueprint = ...
                obj.generate_fracture_structure_description(...
                blueprint_range);
        end
    end

    methods (Access = private)
        function [blueprint] = ...
                generate_fracture_info(obj, blueprint_sections)
            blueprint = obj.blueprint;
            for index = 1:size(blueprint_sections,1)
                blueprint.(blueprint_sections(index,1)) = ...
                    obj.fracture_input_raw.read_mapped_data ...
                    (blueprint_sections(index,2));
            end
        end
        function blueprint = ...
                generate_fracture_structure_description(obj,description)
            blueprint = obj.blueprint;
            data = obj.fracture_input_raw.read_range_data(description(2));
            num_frac = size(data,1);
            frac = struct('isActive',[], 'CenterCoord', [[],[],[]]);
            fracs = repelem(frac,num_frac,1);
            fracs = set_field(fracs,"isActive", ...
                obj.process_scalar(data(:,1)));
            fracs = set_field(fracs,"CenterCoord", ...
                obj.process_vec(data(:,2)));
            fracs = set_field(fracs,"Length", ...
                obj.process_scalar(data(:,3)));
            fracs = set_field(fracs,"Width", ...
                obj.process_scalar(data(:,4)));
            fracs = set_field(fracs,"RotAngleAlongZ", ...
                obj.process_scalar(data(:,5)));
            fracs = set_field(fracs,"RotAngleAlongL", ...
                obj.process_scalar(data(:,6)));
            fracs = set_field(fracs,"RotAngleAlongW", ...
                obj.process_scalar(data(:,7)));
            fracs = set_field(fracs,"GridResRatio", ...
                obj.process_scalar(data(:,8)));
            fracs = set_field(fracs,"GridNumAlongL", ...
                obj.process_scalar(data(:,9)));
            fracs = set_field(fracs,"GridNumAlongW", ...
                obj.process_scalar(data(:,10)));
            fracs = set_field(fracs,"Aperture", ...
                obj.process_scalar(data(:,11)));
            fracs = set_field(fracs,"Porosity", ...
                obj.process_scalar(data(:,12)));
            fracs = set_field(fracs,"Permeability", ...
                obj.process_scalar(data(:,13)));
            adm_expr = join(convertCharsToStrings(data(:,14:17)),",");
            [fracs] = set_field_eval(fracs, "ADM", adm_expr);

            %             fracs = set_field_struct(fracs,"ADM", ...
%                 obj.process_adm(adm_expr,4));

            mat_expr = join(convertCharsToStrings(data(:,18:24)),",");
            % Create a matrix of strings
            [fracs] = set_field_eval(fracs, "CornerCoords", mat_expr);
            blueprint.(description(1)) = fracs;
        end
        function [scalar] = process_scalar(obj,string_array)
            string_array = convertCharsToStrings(string_array);
            LX = obj.blueprint.Size(1);
            LY = obj.blueprint.Size(2);
            LZ = obj.blueprint.Size(3);
            num_frac = size(string_array,1);
            scalar = zeros(num_frac,1);
            for index = 1:num_frac
                tmp = eval(string_array(index));
                scalar(index,:) = tmp';
            end
        end
        function [adm] = process_adm(obj,string_array, varargin)
            txt = replace(string_array, " ", ",");
            txt = erase(txt, "]" | "[" );
            txt = txt.extractBetween(1+1,txt.strlength-1);
            txt = split(txt, ",");
            adm_matrix = cellfun(@str2num,txt);
            adm(size(adm_matrix,1)) = struct();
            for index = 1:size(adm_matrix)
                adm(index).Active = adm_matrix(index,1);
                adm(index).CoarseLevel = adm_matrix(index,2);
                adm(index).CoarseRatio = adm_matrix(index,3:end);
            end
        end
        function [coords] = process_vec(obj,string_array, varargin)
            if isempty(varargin)
                n = 3;
            else
                n = varargin{1};
            end
            string_array = convertCharsToStrings(string_array);
            LX = obj.blueprint.Size(1);
            LY = obj.blueprint.Size(2);
            LZ = obj.blueprint.Size(3);
            num_frac = size(string_array,1);
            coords = zeros(num_frac,n);
            for index = 1:num_frac
                tmp = eval(string_array(index));
                coords(index,:) = tmp';
            end
        end
    end
end



function value = find_str_interval(strg, interval)
n = size(strg,1);
m = size(interval,1);
if n ~= m
    error("The dimmensions of the string and interval are not equal.");
end
value = cell(n,1);
for index = 1:n
    line_strg = convertStringsToChars(strg(index));
    value{index,1} = line_strg(interval(index,1):interval(index,2));
end
value = convertCharsToStrings(value);
end

function [array_struct] =  set_field_struct(array_struct, field, input)
n = size(array_struct,1);
for index = 1:n
    array_struct(index).(field) = input(index);
end
end

function [new_array] =  set_field(array_struct, field, input)
n = size(array_struct,1);
m = size(input,1);
if n ~= m
    error("The dimmensions of the input and the " + ...
        "array struct are not equal.")
end
for index = 1:n
    array_struct(index).(field) = input(index,:);
end
new_array = array_struct;
end

function [new_array] =  set_field_eval(array_struct, field, input)
n = size(array_struct,1);
m = size(input,1);
if n ~= m
    error("The dimmensions of the input and the " + ...
        "array struct are not equal.")
end
for index = 1:n
    matrix = eval(input(index,:));
    array_struct(index).(field) = matrix;
end
new_array = array_struct;
end