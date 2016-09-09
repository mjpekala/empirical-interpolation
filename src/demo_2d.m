% DEMO_2D   Simple example of "magic point" interpolation
%
% References:
%   Maday et al. "A general multipurpose interpolation procedure:
%                 the magic points," CPAA 2009.

% mjp, sept 2016


%% construct a toy function 
f = @(X) cos(3*X(:,1)) .* sin(2*X(:,2));
f_rescaled = @(X) f(X*pi);

n_max = 14;
U_ = make_polynomial_basis(2, n_max);
[Omega, domain_info] = make_domain_2d(.02, 'triangle');


f_true = f_rescaled(Omega);

figure; 
mesh2(f_true, domain_info);
xlabel('x'); ylabel('y'); 
title('f_{true}');


%% interpolate

tic
[s, Lambda_M] = choose_magic(Omega, U_);
toc
fprintf('[%s]: Lambda_M = %0.3f\n', mfilename, Lambda_M);

tic
f_hat = interp_magic(f_rescaled, s);
toc


figure; 
mesh2(f_hat, domain_info);
xlabel('x'); ylabel('y'); 
title('$\hat{f}$', 'interpreter', 'latex');
hold on;
stem3(s.Omega(s.x,1), s.Omega(s.x,2), f_rescaled(s.Omega(s.x,:)), 'ro');
hold off;


%% the interpolation error
figure; 
mesh2(abs(f_true-f_hat), domain_info);
title('error');
colorbar();
hold on;
stem3(s.Omega(s.x,1), s.Omega(s.x,2), ...
      abs(f_rescaled(Omega(s.x,:))- f_hat(s.x)), 'ro');
hold off;


%% the magic points
figure;
plot(s.Omega(s.x,1), s.Omega(s.x,2), 'x');
grid on;
xlabel('x'); ylabel('y'); title('magic point locations');
