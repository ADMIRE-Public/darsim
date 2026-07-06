function [l2norm] = halfvecnorml2(vecs)
%halfvecnorml2 Calculates the L2 norm of the vectors along the rows of vec
%   Function compatible with half class
l2norm = (dot(vecs,vecs,2)).^.5;
end