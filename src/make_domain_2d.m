function [Omega, U_] = make_domain_2d(n_max, delta)
% MAKE_DOMAIN_2D   Creates example spatial and interpolation domains.
%
%  Creates a simple 2d domain (compact subset of R^2) and a polynomial
%  function space for interpolation.  If you want a more interesting
%  domain (e.g. triangular, hexoganal, higher-dimensional) or a
%  different function space (e.g. Legandre polynomials) you'll need
%  different code.
%
%     n_max := the maximum polynomial degree to use when
%              constructing U_.
%     delta := resolution of the interpolation domain.
 
if nargin < 2, delta=.1; end

x = -1:delta:1;
[X,Y] = meshgrid(x, x);
Omega = [X(:) Y(:)];

% Construct the polynomial basis:
%
%    W_n(Omega) := { (x^i)*(y^j), (x,y) \in Omega, i+j <= n },  0 <= n <= n_max
%
% with m = (n+1)(n+2)/2 elements.
%
m = (n_max+1)*(n_max+2) / 2;   
U_ = cell(m,1);  % underscore is my notation indicating this is not a matrix
for kk = 1:m
    U_{kk} = @(X2d)  X2d(:,1).^ii(kk) .* X2d(:,2).^jj(kk);
end
