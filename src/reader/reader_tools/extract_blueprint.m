function [blueprint] = extract_blueprint(folder,keys,varargin)
%extract_blueprint search input folder building a blueprint based on keys
all_folder_files = properties(folder);
if isempty(varargin)
    folder_files = all_folder_files;
else
    folder_files = varargin;
    is_contained = ismember(folder_files,all_folder_files);
    if ~all(is_contained)
       error("ERROR: %s not found inside this folder.", ...
         strjoin([convertCharsToStrings(folder_files(~is_contained))]," "))
    end
    folder_files = folder_files';
end
blueprint = struct();
for folder_label = folder_files'
    local_folder = folder.(folder_label{1});
    for key  = keys
        if isfield(local_folder,key)
            blueprint.(key) = local_folder.(key);
        end
    end
end
blueprint_fields = convertCharsToStrings(fieldnames(blueprint));
ref_found = ismember(keys,blueprint_fields);
if ~all(ref_found)
    fprintf("WARNING: Fields %s not found.\n", ...
                                        strjoin([keys(~ref_found)], ","))
end
end