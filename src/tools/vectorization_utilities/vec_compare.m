function [full_comp] = vec_compare(op,vec_query,comp_vec, tol)
%vec_compare
%   INPUT:
%       op -> handle of the function operatation needed, i.e. @geq , @leq
%       vec_query -> vector query with n vectors stacked on rows
%       comp_vec -> comparison vector
%   OUTPUT:
%       full_comp -> result of the comparsion column wise between 
%       vec_querry and comp_vec
full_comp = false(size(vec_query));
for index = 1:size(full_comp,2)
    full_comp(:,index) = op(vec_query(:,index),comp_vec(index) + tol);
end
end