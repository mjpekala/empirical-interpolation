% DEMO_2D   Simple example of "magic point" interpolation
%
% References:
%   Maday et al. "A general multipurpose interpolation procedure:
%                 the magic points," CPAA 2009.

m = 5;

%% construct a toy function
f = @(X) cos(3*X(:,1)) .* sin(2*X(:,2));
f_rescaled = @(X) f(X*pi);
[Omega, U_] = make_domain_2d(5, .02);

%% visualize the function of interest 
Z = f_rescaled(Omega);
n = sqrt(length(Z));
X = reshape(Omega(:,1), n, n);
Y = reshape(Omega(:,2), n, n);
Z = reshape(Z, n, n);
figure; surf(X, Y, Z); title('f_true');


%% interpolate the function of interest
tic
[f_interp, s] = interp_magic(f_rescaled, Omega, U_, m);
toc

%% visualize result
Zhat = reshape(f_interp(Omega), n, n);
figure; surf(X, Y, Zhat); title('$\hat{f}$', 'interpreter', 'latex');

figure; surf(X, Y, abs(Z-Zhat)), title('error (L2)');

%% Some diagnostics
for k = 1:m
    Ui = reshape(s.u{k}(Omega), n, n);
    figure; surf(X, Y, Ui);  title(sprintf('U_{%d}', k));
end

