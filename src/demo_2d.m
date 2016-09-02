% DEMO_2D   Simple example of "magic point" interpolation
%
% References:
%   Maday et al. "A general multipurpose interpolation procedure:
%                 the magic points," CPAA 2009.

% construct a toy function
f = @(x,y) cos(3*x) .* sin(2*y);
f_rescaled = @(x,y) f(x*pi,y*pi);
[Omega, U_] = make_domain_2d(5, .02);
Z = f_rescaled(Omega(:,1), Omega(:,2));

% visualize the function of interest and a basis function 
n = sqrt(length(Z));
X = reshape(Omega(:,1), n, n);
Y = reshape(Omega(:,2), n, n);
Z = reshape(Z, n, n);
figure; surf(X, Y, Z); xlabel('x'); ylabel('y');

U1 = reshape(U_{1}(Omega), n, n);
figure; surf(X, Y, U1); xlabel('x'); ylabel('y'); title('U_1');

% interpolate the function of interest
f_interp = interp_magic(f_rescaled, Omega, U_, 10);
