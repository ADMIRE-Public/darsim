function out = transpose1dvec(mat)
a = mat(:,1);
b = mat(:,2);
c = mat(:,3);
d = mat(:,4);
e = mat(:,5);
f = mat(:,6);
g = mat(:,7);
h = mat(:,8);
i = mat(:,9);
out = [a, d, g, b, e, h, c, f, i];
end