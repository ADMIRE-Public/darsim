function [reduced_operation] = reduced_jagged_array(...
                adj_matrix,matrix_jagged_array, num_col, ref_fun, varargin)
%reduced_jagged_array performs a reduction operation use the ref_fun using 
%on the matrix_jagged array using as base the adj_matrix
%   
if num_col < 1
    error("The number of columns of the resulting operation ..." + ...
        "needs to be greater than 0.")
end
row_sum = sum(adj_matrix,2);
reduced_operation = zeros(size(adj_matrix,1),num_col);
for index = unique(row_sum)'
    ref = row_sum == index;
    local_matrix = matrix_jagged_array(ref,:)';
    %local_ref = adj_matrix(ref,:);
    local_ref = find(local_matrix);
    ref_fun_arg = full(reshape(local_matrix(local_ref),index,[])'); %#ok<FNDSB> 
    reduced_operation(ref) = ref_fun(ref_fun_arg, ref, varargin{:});
end
end