% DEMO_2D   Simple example of "magic point" interpolation
%
% References:
%   Maday et al. "A general multipurpose interpolation procedure:
%                 the magic points," CPAA 2009.

n_max = 12;

%% construct a toy function & interpolate
f = @(X) cos(3*X(:,1)) .* sin(2*X(:,2));
f_rescaled = @(X) f(X*pi);
[Omega, U_] = make_domain_2d(n_max, .02);

tic
s = choose_magic(Omega, U_);
toc

tic
f_hat = interp_magic(f_rescaled, s);
toc


%% visualize 
Z = f_rescaled(Omega);
n = sqrt(length(Z));
Z = reshape(Z, n, n);

X = reshape(Omega(:,1), n, n);
Y = reshape(Omega(:,2), n, n);
Z_hat = reshape(f_hat, n, n);

figure; 
surf(X, Y, Z); title('f_true');
hold on;
stem3(s.Omega(s.x,1), s.Omega(s.x,2), f_hat(s.x), 'ro');
hold off;


% the interpolation and error
figure; 
surf(X, Y, Z_hat); title('$\hat{f}$', 'interpreter', 'latex');
hold on;
stem3(s.Omega(s.x,1), s.Omega(s.x,2), f_rescaled(s.Omega(s.x,:)), 'ro');
hold off;

figure; surf(X, Y, abs(Z-Z_hat)), title('error (L2)');
