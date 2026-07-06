function [fullmat] = sparse_to_full_zeromat(sparse_mat)
%sparse_to_full_zeromat transform bool adj matrix in a full matrix.
[ii,~] = find(sparse_mat');
auxmat = reshape_jagged_array_in_matrix(sparse_mat,ii);
fullmat = sparse_to_zero_ind_mat(auxmat);
end