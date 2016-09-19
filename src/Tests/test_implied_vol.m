% This example taken from class notes of N. Christou 
%    (UCLA, Stats C183/C283)

% mjp, sept 2016

addpath('..');

%% A single configuration with a known result
c = 1.875;
S0 = 21;
K = 20;
r = 0.1;
t = 0.25;

sigma = implied_vol(c, S0, K, r, t);
assert(abs(sigma-0.2345) < 1e-4);

sigma = implied_vol([c c], [S0 S0], [K K], [r r], [t t]);
assert(length(sigma) == 2);
assert(all(abs(sigma-0.2345) < 1e-4));


%% an example of evaluating multiple points

c_vals = (0:.1:6)';
K_vals = (18:.1:25)';

[C_mat,K_mat] = meshgrid(c_vals, K_vals);

onez = ones(numel(C_mat),1);
S0 = 21 * onez;
r = 0.1 * onez;
t = 0.25 * onez;

tic
sigma = implied_vol(C_mat(:), S0, K_mat(:), r, t);
fprintf('[%s]: time to calculate %d values: %0.2f sec\n', mfilename, numel(onez), toc);

Z = reshape(sigma, size(C_mat,1), size(C_mat,2));
figure; mesh(C_mat,K_mat,Z);
xlabel('c'); ylabel('k');

