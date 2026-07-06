function [out] = create_matedgezero(mat)
out = [mat zeros(size(mat,1),1)];
[~, jj] = min(out,[],2);
ii = [1:size(out,1)]';
linear_ind = sub2ind(size(out),ii,jj);
out(linear_ind) = out(:,1);
end