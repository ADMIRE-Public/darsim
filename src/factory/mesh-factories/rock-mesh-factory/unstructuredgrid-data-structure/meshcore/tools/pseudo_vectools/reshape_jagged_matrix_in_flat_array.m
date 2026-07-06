function [flat_jagged_array] = reshape_jagged_matrix_in_flat_array(...
                                              adj_matrix,jagged_matrix)
%reshape_jagged_array_in_matrix Transform a sparse matrix that represents a
% jagged matrix into a flat jagged array according to the adj_matrix.
[ref_global] = find(adj_matrix');
t_jagged_matrix = jagged_matrix';
flat_jagged_array = t_jagged_matrix(ref_global)';
end