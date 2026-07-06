function [zero_ind] = sparse_to_zero_ind_mat(sp_mat)
num_colum = sum(sp_mat ~= 0 ,2);
max_colum = max(num_colum);
zero_ind = zeros(size(sp_mat,1),max_colum);
for index = 1:max_colum
    ref = num_colum == index;
    auxmat = sp_mat(ref,:)';
    auxmat = auxmat(auxmat~= 0);      
    zero_ind(ref,1:index) = reshape(auxmat,index,[])';
end
end