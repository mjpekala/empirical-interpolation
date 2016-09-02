function A = apply(f_cell, X)
% APPLY Apply an array of functions to one or more inputs.
%
%    f_cell := a cell array of k function handles, where each
%              function maps R^n -> R (i.e. is a scalar-valued 
%              function of n inputs).  It is assumed that 
%              each function in the array can support multiple 
%              inputs; the presumed API is:
%
%                        c = f_cell{ii}(X)
%
%              where X is an (m x n) matrix of inputs and c is m x 1.
%
%    X      := the set of points that should be evaluated by each f \in f_cell.
%              a (m x n) matrix of m points each having n dimensions.
%
%    A      := an (m x k) matrix containing all the function evaluations.

%  Example:
%{
    f = {};
    f{1} = @(X) X(:,1) + X(:,2);
    f{2} = @(X) X(:,1).*X(:,2);
    f{3} = @(X) X(:,1).^2 + X(:,2).^2;
    X = [1 1 ; 2 2; 3 3];
    Z = apply(f, X)
%}


k = length(f_cell);
[m,n] = size(X);

A = zeros(m,k);
parfor kk = 1:k
    A(:,kk) = feval(f_cell{kk}, X);
end
