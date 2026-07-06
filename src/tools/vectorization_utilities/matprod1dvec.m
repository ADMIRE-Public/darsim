function out = matprod1dvec(mat1,mat2)
a1 = mat1(:,1);
b1 = mat1(:,2);
c1 = mat1(:,3);
d1 = mat1(:,4);
e1 = mat1(:,5);
f1 = mat1(:,6);
g1 = mat1(:,7);
h1 = mat1(:,8);
i1 = mat1(:,9);
a2 = mat2(:,1);
b2 = mat2(:,2);
c2 = mat2(:,3);
d2 = mat2(:,4);
e2 = mat2(:,5);
f2 = mat2(:,6);
g2 = mat2(:,7);
h2 = mat2(:,8);
i2 = mat2(:,9);
out = [a1.*a2 + b1.*d2 + c1.*g2, ... 
       a1.*b2 + b1.*e2 + c1.*h2, ...
       a1.*c2 + b1.*f2 + c1.*i2, ...
       a2.*d1 + d2.*e1 + f1.*g2, ...
       b2.*d1 + e1.*e2 + f1.*h2, ...
       c2.*d1 + e1.*f2 + f1.*i2, ...
       a2.*g1 + d2.*h1 + g2.*i1, ...
       b2.*g1 + e2.*h1 + h2.*i1, ...
       c2.*g1 + f2.*h1 + i1.*i2];
end