function [information] = retrieve_field_from_cell(field,cellarray,varargin)
if isempty(varargin)
    n = 0;
else
    n = varargin{1} - 1;
end
ref_start =  find(ismember(cellarray, field)) + 1;
ref_end = ref_start + n;
information = cellArrayToMixed(cellarray(ref_start:ref_end));
end