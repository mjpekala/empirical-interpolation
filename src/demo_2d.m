% DEMO_2D   Simple example of "magic point" interpolation
%
% References:
%   Maday et al. "A general multipurpose interpolation procedure:
%                 the magic points," CPAA 2009.

% mjp, sept 2016


%% construct a toy function on a triangular domain
f = @(X) cos(3*X(:,1)) .* sin(2*X(:,2));
f_rescaled = @(X) f(X*pi);

n_max = 14;
W_n = make_polynomial_basis(2, n_max);
[Omega, domain_info] = make_domain_2d(.02, 'triangle');

compare_methods_2d(f_rescaled, W_n, Omega, domain_info);
