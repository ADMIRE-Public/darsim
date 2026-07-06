function [out] = matvecprod1dvec(mat,vec)
a = mat(:,1);
b = mat(:,2);
c = mat(:,3);
d = mat(:,4);
e = mat(:,5);
f = mat(:,6);
g = mat(:,7);
h = mat(:,8);
i = mat(:,9);
x = vec(:,1);
y = vec(:,2);
z = vec(:,3);
out = [a.*x + b.*y + c.*z, d.*x + e.*y + f.*z, g.*x + h.*y + i.*z];
end