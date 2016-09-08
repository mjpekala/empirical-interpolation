function A = apply(f_cell, X)
% APPLY Apply an array of (basis) functions to one or more data points.
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

% NOTE: a previous version of this software package made more
% extensive use of function handles; I created the apply() in
% anticipation of needing a faster alternative to cellfun.  However, I
% abandoned this fairly quickly in favor of a more purely matrix-basee
% approach.  The apply function is a vestigial remnant of the former
% implementation and can probably be removed/replaced.

% mjp, sept 2016

if 1
    tmp = cellfun(@(f) f(X), f_cell(:)', 'UniformOutput', 0);
    A = cat(2, cell2mat(tmp));
    %assert(size(A,2) == length(f_cell));
    %assert(size(A,1) == size(X,1));
else
    k = length(f_cell);
    [m,n] = size(X);

    % special logic for 1d case
    if m == 1, 
        X = X(:); 
        [m,n] = size(X);
    end;  

    A = zeros(m,k);
    for kk = 1:k
        A(:,kk) = feval(f_cell{kk}, X);
    end
end

