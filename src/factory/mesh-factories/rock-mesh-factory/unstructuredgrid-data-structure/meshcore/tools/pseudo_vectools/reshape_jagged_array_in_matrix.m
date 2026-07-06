function [reshape_jagged_array] = reshape_jagged_array_in_matrix(...
                                              adj_matrix,flat_jagged_array)
%reshape_jagged_array_in_matrix Transform a jagged matrix stored in a flat
%jagged array into a not jagged sparse matrix with based on the
%representation of the jagged array in the adj_matrix.
row_sum = sum(adj_matrix,2);
num_entries = sum(row_sum);
if num_entries ~= max(size(flat_jagged_array))
    error("ERROR: The flat jagged array size is not compatible with" + ...
                                                " the adjacencies matrix")
end
[row_wise_ref, colum_wise_ref] = find(adj_matrix');
reshape_jagged_array = ...
                   sparse(row_wise_ref,colum_wise_ref,flat_jagged_array,size(adj_matrix,2), size(adj_matrix,1))';
end