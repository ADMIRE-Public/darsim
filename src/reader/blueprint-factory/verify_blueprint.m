function [blueprint] = verify_blueprint(InputStruct, defaultKeys, ...
                        defaultValues, varargin)
if isempty(varargin)
    timeLabels = [];
else
    timeLabels = varargin{1};
end
%% verifies and imposes default values for the specific structure
%% also converts time labels
blueprint = add_mandatory_fields(InputStruct,defaultKeys, defaultValues);
blueprint = convert_time_labels(blueprint, timeLabels);
end

function blueprint =  convert_time_labels(blueprint,label)
nlabel = length(label);
if (nlabel == 1) && isstring(label)
    blueprint.(label) =  process_time(blueprint.(label));
else 
    for local_label = label
        blueprint.(local_label) =  process_time(blueprint.(local_label));
    end
end
end

function time = process_time(time_var)
 % Initialize variables to store the extracted values
    time = struct();
    [time.days, time.hours, time.minutes, time.seconds] = ...
        parseTimeString(time_var);
    num_time = time.days*24*3600 +time. hours*3600 + time.minutes*60 + ...
        time.seconds;
    time.total_time = num_time;
end


function [days, hours, minutes, seconds] = parseTimeString(inputString)
    % Split the input string into words
    words = strsplit(inputString);

    % Initialize variables to store the values
    days = 0;
    hours = 0;
    minutes = 0;
    seconds = 0;

    % Loop through the words and extract values when present
    for i = 1:length(words)
        if strcmpi(words{i}, 'DAYS')
            if i > 1
                days = str2double(words{i-1});
            end
        elseif strcmpi(words{i}, 'HOURS')
            if i > 1
                hours = str2double(words{i-1});
            end
        elseif strcmpi(words{i}, 'MINUTES')
            if i > 1
                minutes = str2double(words{i-1});
            end
        elseif strcmpi(words{i}, 'SECONDS')
            if i > 1
                seconds = str2double(words{i-1});
            end
        end
    end
end