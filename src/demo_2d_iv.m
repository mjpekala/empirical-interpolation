% Another 2d example.

% mjp, sept 2016

%% Setup
c_min = 1.5;  c_max = 3.5;
K_min = 19;   K_max = 21;
delta = .02;

S0 = 21;
r = 0.1;
t = 0.25;

% This scales c, K from [-1 1] to the actual values and
% builds a matrix of values to compute.
iv_scaled = @(X) implied_vol(1+c_min+((c_max - c_min)/2)*X(:,1), ...
                             ones(size(X,1),1)*S0, ...
                             1+K_min+((K_max - K_min)/2)*X(:,2), ...
                             ones(size(X,1),1)*r, ...
                             ones(size(X,1),1)*t);

n_max = 10;
W_n = make_polynomial_basis(2, n_max);
[Omega, x, y, idx] = make_domain_2d(delta, 'square');

to_square = @(v) reshape(v, sqrt(numel(x)), sqrt(numel(x)));


%% Brute force calculation
tic
sigma = iv_scaled(Omega);
toc

figure;
mesh(to_square(x), to_square(y), to_square(sigma));
xlabel('c'); ylabel('k'); 
title('$\sigma_{true}$', 'interpreter', 'latex');

%% Empirical interpolation
tic
[s, Lambda_M] = choose_magic(Omega, W_n);
toc
fprintf('[%s]: Lambda_M = %0.3f\n', mfilename, Lambda_M);

tic
[sigma_hat, sigma_magic] = interp_magic(iv_scaled, s);
toc




%% Matlab's built-in interpolation
tic
sigma_hat_linear = griddata(s.Omega(s.x,1), s.Omega(s.x,2), sigma_magic, ...
                            s.Omega(:,1), s.Omega(:,2), ...
                            'linear');
toc

tic
sigma_hat_cubic = griddata(s.Omega(s.x,1), s.Omega(s.x,2), sigma_magic, ...
                            s.Omega(:,1), s.Omega(:,2), ...
                            'cubic');
toc


%% Visualize fits

figure('Position', [100 100 1200 450]); 

ax1 = subplot(1,3,1);
mesh(to_square(x), to_square(y), to_square(sigma_hat)); 
xlabel('c'); ylabel('k'); 
title('$\hat{\sigma}$', 'interpreter', 'latex');
hold on;
stem3(s.Omega(s.x,1), s.Omega(s.x,2), sigma_magic, 'ro');
hold off;

ax2 = subplot(1,3,2); 
mesh(to_square(x), to_square(y), to_square(sigma_hat_linear)); 
xlabel('c'); ylabel('k'); 
title('$\hat{\sigma}_{linear}$', 'interpreter', 'latex');
hold on;
stem3(s.Omega(s.x,1), s.Omega(s.x,2), sigma_magic, 'ro');
hold off;

ax3 = subplot(1,3,3);
mesh(to_square(x), to_square(y), to_square(sigma_hat_cubic)); 
xlabel('c'); ylabel('k'); 
title('$\hat{\sigma}_{cubic}$', 'interpreter', 'latex');
hold on;
stem3(s.Omega(s.x,1), s.Omega(s.x,2), sigma_magic, 'ro');
hold off;

linkprop([ax1 ax2 ax3], 'CameraPosition');


%% Error analysis

figure('Position', [100 100 1200 450]); 

ax1 = subplot(1,3,1);
mesh(to_square(x), to_square(y), to_square(abs(sigma - sigma_hat))); 
xlabel('c'); ylabel('k'); 
title('err (magic points)');

ax2 = subplot(1,3,2);
mesh(to_square(x), to_square(y), to_square(abs(sigma - sigma_hat_linear))); 
xlabel('c'); ylabel('k'); 
title('err (linear)');

ax3 = subplot(1,3,3);
mesh(to_square(x), to_square(y), to_square(abs(sigma - sigma_hat_cubic))); 
xlabel('c'); ylabel('k'); 
title('err (cubic)');

linkprop([ax1 ax2 ax3], 'CameraPosition');
