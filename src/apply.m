function A = apply(f_cell, X)
% APPLY Apply an array of functions to one or more inputs.
%
%    f_cell := a cell array of k function handles, where each
%              function maps R^n -> R (i.e. is a scalar-valued function of
%              n inputs).  It is assumed here that each function
%              in the array can deal with multiple inputs; the
%              presumed API is:
%
%                        c = f_cell{ii}(X)
%
%              where X is an (m x n) matrix of inputs and c is m x 1.
%
%    A      := an (m x k) matrix containing the function evaluations

k = length(f_cell);
[m,n] = size(X);
A = zeros(m,k);
parfor kk = 1:k
    A(:,kk) = f_cell{kk}(X);
end
