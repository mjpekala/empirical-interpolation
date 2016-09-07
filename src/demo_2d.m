% DEMO_2D   Simple example of "magic point" interpolation
%
% References:
%   Maday et al. "A general multipurpose interpolation procedure:
%                 the magic points," CPAA 2009.

n_max = 14;

%% construct a toy function & interpolate
f = @(X) cos(3*X(:,1)) .* sin(2*X(:,2));
f_rescaled = @(X) f(X*pi);

U_ = make_polynomial_basis(2, n_max);
[Omega, x, y, idx] = make_domain_2d(.02, 'shape', 'triangle');

tic
[s, Lambda_M] = choose_magic(Omega, U_);
toc
fprintf('[%s]: Lambda_M = %0.3f\n', mfilename, Lambda_M);

tic
f_hat = interp_magic(f_rescaled, s);
toc

to_square = @(v) reshape(v, sqrt(numel(x)), sqrt(numel(x)));

%% visualize 

z = NaN*ones(size(x));
z(idx) = f_rescaled(Omega);

z_hat = NaN*ones(size(x));
z_hat(idx) = f_hat;

figure; 
mesh(to_square(x), to_square(y), to_square(z));
%hold on;
%stem3(s.Omega(s.x,1), s.Omega(s.x,2), f_hat(s.x), 'ro');
%hold off;
title('f_{true}');


% the interpolation 
figure; 
surf(to_square(x), to_square(y), to_square(z_hat)); 
title('$\hat{f}$', 'interpreter', 'latex');
hold on;
stem3(s.Omega(s.x,1), s.Omega(s.x,2), f_rescaled(s.Omega(s.x,:)), 'ro');
hold off;

% the interpolation error
figure; 
surf(to_square(x), to_square(y), to_square(abs(z-z_hat)));
title('error (L2)');
colorbar();
hold on;
stem3(s.Omega(s.x,1), s.Omega(s.x,2), ...
      abs(f_rescaled(Omega(s.x,:))- f_hat(s.x)), 'ro');
hold off;
