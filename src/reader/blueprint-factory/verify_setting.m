function [verifiedSettings] = verify_InputSettings(InputSettings)
verifiedSettings = InputSettings;
defaultKeys = ["GRID_DISCRETIZATION", "FRACTURED" ];
defaultValues = {"CartesianGrid", false};
for index = 1:length(defaultValues)
    localKey = defaultKeys(index);
    localValue = defaultValues(index);
    if ~isfield(InputSettings,localKey)
        verifiedSettings.(localKey) = localValue{1};
    end
end
verifiedSettings.("MINMAXDT")  = process_min(InputSettings.("MINMAXDT"));


end

function out = process_min(str)
    out = str2double(convertCharsToStrings((strsplit(str))));
end