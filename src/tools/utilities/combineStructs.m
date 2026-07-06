function combinedStruct = combineStructs(struct1, struct2)
    % Create a new structure to hold the combined data
    combinedStruct = struct();

    % Copy fields from struct1
    fields1 = fieldnames(struct1);
    for i = 1:numel(fields1)
        combinedStruct.(fields1{i}) = struct1.(fields1{i});
    end

    % Copy fields from struct2
    fields2 = fieldnames(struct2);
    for i = 1:numel(fields2)
        combinedStruct.(fields2{i}) = struct2.(fields2{i});
    end
end