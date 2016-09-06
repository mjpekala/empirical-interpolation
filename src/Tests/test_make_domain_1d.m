% TEST_MAKE_DOMAIN_1d

addpath('..')

[Omega, W_n] = make_domain_1d(10,.01);

figure;
plot(Omega, W_n{1}(Omega), ...
     Omega, W_n{2}(Omega), ...
     Omega, W_n{3}(Omega), ...
     Omega, W_n{4}(Omega));
legend('x^0', 'x^1', 'x^2', 'x^3', 'Location', 'SouthEast');
title('first four polynomial basis functions (1d)');

