function output_variable = get_variable(input_variable)
    % Check if the input is already numerical (matrix or scalar)
    if isnumeric(input_variable)
        output_variable = input_variable;
        return;
    end
    
    % Check if the input is a string
    if ischar(input_variable)
        % Split the input string into parts based on spaces or commas
        parts = strsplit(input_variable, {' ', ','});
        
        % Convert the parts to numbers
        nums = str2double(parts);
        
        % Check if all parts were successfully converted to numbers
        if ~any(isnan(nums))
            % If conversion was successful, return the numerical array
            output_variable = nums;
        else
            % If some parts are not numeric, return an error message
            error('Input string cannot be converted to a numerical array.');
        end
    else
        % If the input is neither a string nor a numerical value, return an error message
        error('Input must be a string, matrix, scalar, or a numerical array.');
    end
end