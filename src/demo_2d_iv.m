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

tic
[s, Lambda_M] = choose_magic(Omega, W_n);
toc
fprintf('[%s]: Lambda_M = %0.3f\n', mfilename, Lambda_M);

tic
sigma_hat = interp_magic(iv_scaled, s);
toc

figure; 
mesh(to_square(x), to_square(y), to_square(sigma_hat)); 
xlabel('c'); ylabel('k'); 
title('$\hat{\sigma}$', 'interpreter', 'latex');
hold on;
stem3(s.Omega(s.x,1), s.Omega(s.x,2), iv_scaled(s.Omega(s.x,:)), 'ro');
hold off;


figure; 
mesh(to_square(x), to_square(y), to_square(abs(sigma - sigma_hat))); 
xlabel('c'); ylabel('k'); 
title('err');
