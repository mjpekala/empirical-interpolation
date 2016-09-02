% DEMO_2D   Simple example of "magic point" interpolation
%
% References:
%   Maday et al. "A general multipurpose interpolation procedure:
%                 the magic points," CPAA 2009.

% construct a toy function
f = @(X) cos(3*X(:,1)) .* sin(2*X(:,2));
f_rescaled = @(X) f(X*pi);
[Omega, U_] = make_domain_2d(5, .02);

% visualize the function of interest and a basis function 
Z = f_rescaled(Omega);
n = sqrt(length(Z));
X = reshape(Omega(:,1), n, n);
Y = reshape(Omega(:,2), n, n);
Z = reshape(Z, n, n);
figure; surf(X, Y, Z); title('f_true');

U1 = reshape(U_{1}(Omega), n, n);
figure; surf(X, Y, U1); title('U_1');

% interpolate the function of interest
f_interp = interp_magic(f_rescaled, Omega, U_, 10);

Zhat = reshape(f_interp(Omega), n, n);
figure; surf(X, Y, Zhat); title('$\hat{f}$', 'interpreter', 'latex');
