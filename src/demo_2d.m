% DEMO_2D   Simple example of "magic point" interpolation
%
% References:
%   Maday et al. "A general multipurpose interpolation procedure:
%                 the magic points," CPAA 2009.

m = 15;

%% construct a toy function & interpolate
f = @(X) cos(3*X(:,1)) .* sin(2*X(:,2));
f_rescaled = @(X) f(X*pi);
[Omega, U_] = make_domain_2d(5, .02);

tic
[f_interp, s] = interp_magic(f_rescaled, Omega, U_, m);
toc


%% visualize 
Z = f_rescaled(Omega);
n = sqrt(length(Z));
Z = reshape(Z, n, n);

tic
Zhat = reshape(f_interp(Omega), n, n);
toc

X = reshape(Omega(:,1), n, n);
Y = reshape(Omega(:,2), n, n);

figure; 
surf(X, Y, Z); title('f_true');
hold on;
stem3(s.x(:,1), s.x(:,2), f_interp(s.x), 'ro');
hold off;


% the interpolation and error

figure; 
surf(X, Y, Zhat); title('$\hat{f}$', 'interpreter', 'latex');
hold on;
stem3(s.x(:,1), s.x(:,2), f_interp(s.x), 'ro');
hold off;

figure; surf(X, Y, abs(Z-Zhat)), title('error (L2)');

% some other diagnostics
for k = 1:m
    Ui = reshape(s.u{k}(Omega), n, n);
    figure; surf(X, Y, Ui);  title(sprintf('U_{%d}', k));
end

