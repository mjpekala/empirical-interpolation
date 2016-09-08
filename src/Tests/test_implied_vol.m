% This example taken from class notes of N. Christou 
%    (UCLA, Stats C183/C283)

% mjp, sept 2016

c = 1.875;
S0 = 21;
K = 20;
r = 0.1;
t = 0.25;

sigma = implied_vol(c, S0, K, r, t);
assert(abs(sigma-0.2345) < 1e-4);

