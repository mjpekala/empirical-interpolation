function [Omega, W_n] = make_domain_2d(n_max, delta)
% MAKE_DOMAIN_2D   Creates example spatial and interpolation domains.
%
%  Creates a simple 2d domain (compact subset of R^2) and a polynomial
%  function space for interpolation.  If you want a more interesting
%  domain (e.g. triangular, hexoganal, higher-dimensional) or a
%  different function space (e.g. Legandre polynomials) you'll need
%  different code.
% 
%  The polynomial space we create is the following:
%     W_n(Omega) := { (x^i)*(y^j), (x,y) \in Omega, i+j <= n },  0 <= n <= n_max
%  which has m = (n+1)(n+2)/2 elements.
%
%  Parameters:
%     n_max := the maximum polynomial degree to use when
%              constructing W_n.
%     delta := resolution of the interpolation domain.
%
%  References:
%    Maday et al. "A general multipurpose interpolation procedure:
%                  the magic points," CPAA 2009.
 
if nargin < 2, delta=.1; end
assert(delta > 0);
assert(n_max > 1);

m = (n_max+1)*(n_max+2) / 2;   

x = -1:delta:1;
[X,Y] = meshgrid(x, x);
Omega = [X(:) Y(:)];


% the code below is just a quick way to enumerate the space:
%   ii+jj <= n      for 0 <= n <= n_max
Tmp = tril(ones(n_max+1,n_max+1));
[ii,jj] = ind2sub(size(Tmp), find(Tmp));
jj = jj - 1;
ii = n_max+1-ii;
assert(length(ii) == m);
clear Tmp;


W_n = cell(m,1);  % underscore is my notation indicating this is not a matrix
for kk = 1:m
    W_n{kk} = @(X2d)  X2d(:,1).^ii(kk) .* X2d(:,2).^jj(kk);
end
