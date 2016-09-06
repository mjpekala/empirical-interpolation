function [Omega, W_n] = make_domain_1d(n_max, delta)
% MAKE_DOMAIN_1D   Creates domain and interpolation functions.
%
%  The polynomial space we create is the following:
%     W_n(Omega) := { x^i, x \in Omega, i <= n },  0 <= n <= n_max
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

Omega = -1:delta:1;


W_n = cell(n_max+1,1);
for ii = 0:n_max
    W_n{ii+1} = @(x)  x.^ii;
end
