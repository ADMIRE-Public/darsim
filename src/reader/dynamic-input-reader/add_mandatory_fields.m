function [blueprint] = add_mandatory_fields(InputStruct,defaultKeys, ...
    defaultValues )
blueprint = InputStruct;
for index = 1:length(defaultValues)
    localKey = defaultKeys(index);
    localValue = defaultValues(index);
    if ~isfield(InputStruct,localKey)
        blueprint.(localKey) = localValue{1};
    end
end
end