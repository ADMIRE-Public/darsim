function [varargout] = unpack_collumn(x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    nout = size(x,2);
    for k = 1:nout
        varargout{k} = x(:,k);
    end
end