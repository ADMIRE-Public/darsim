function [matrix_jagged_array] = mapped_jagged_array(...
                        adj_matrix,ref_fun, varargin)
%mapped_jagged_array performs a mapping operation use the ref_fun using 
%on the matrix_jagged array using as base the adj_matrix
%
[neigh_ref,~] =find(adj_matrix');
matrix_jagged_array = ...
   reshape_jagged_array_in_matrix(adj_matrix,neigh_ref);
row_sum = sum(adj_matrix,2);
total_num_entries  = sum(row_sum);
[row_ind, col_ind] = find(adj_matrix);
flat_val = zeros(total_num_entries,1);
flat_ind = zeros(total_num_entries,1);
flat_pointer = 1;
global_ind_matrix = sparse(col_ind,row_ind, ...
                            sub2ind(size(adj_matrix),row_ind,col_ind))';
for index = unique(row_sum)'
    ref = row_sum == index;
    local_matrix = matrix_jagged_array(ref,:)';
    local_global_ref_matrix = global_ind_matrix(ref,:)';
    local_ref = find(local_matrix);
    ref_fun_arg = full(reshape(local_matrix(local_ref),index,[])'); 
    local_mapped = ref_fun(ref_fun_arg, ref, varargin{:});
    local_val = reshape(local_mapped',1,[])';
    local_ind = local_global_ref_matrix(local_ref);
    n_entries = size(local_val,1);
    flat_val(flat_pointer:(flat_pointer+n_entries-1)) = local_val;
    flat_ind(flat_pointer:(flat_pointer+n_entries-1)) = local_ind;
    flat_pointer = flat_pointer + n_entries;
end
    [ref_row, ref_col] = ind2sub(size(adj_matrix), flat_ind);
    matrix_jagged_array = sparse(ref_row, ref_col, flat_val);
end