function [out] =  is_above_2d(c1,c2, points,in_point)
tol = 10^-12;
out = false(size(c1,1),1);
den = (c2(:,1) - c1(:,1));
%% m = undefined
ref = abs(den) < tol;

ref_defined = ~ref;

m = (c2(ref_defined,2) - c1(ref_defined,2))./den(ref_defined);
b =  c2(ref_defined,2) - m.*c2(ref_defined,1);

%% analyzing for undefined m
% analyzing Y no inf
refx = (points(ref,2) >= c1(ref,2)) &  (in_point(ref,2) >= c1(ref,2)) | ...
    (points(ref,2) >= c2(ref,2)) &  (in_point(ref,2) >= c2(ref,2));
%% analyzing Y
if ~all(ref)
    R = (m.*points(ref_defined,1) + b) - points(ref_defined,2);
    Rbol = R >= 0;
    Req = (m.*in_point(ref_defined,1) + b) - in_point(ref_defined,2);
    Reqbol = Req >= 0;
    refy = ~xor(Rbol,Reqbol);
    out(ref_defined) = refy;
end
out(ref) = refx;
end