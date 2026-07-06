function isEmpty = isStructEmpty(myStruct)
    % Check if the input is a struct
    if ~isstruct(myStruct)
        error('Input must be a struct');
    end
    
    % Check if the struct is empty
    isEmpty = isempty(fieldnames(myStruct));
end