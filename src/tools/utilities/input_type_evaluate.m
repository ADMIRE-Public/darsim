function [flag] = input_type_evaluate(var)
   %returns 1 if input is text 0 otherwise
   flag =  cellfun(@ischar,var) | cellfun(@isstring,var);
end