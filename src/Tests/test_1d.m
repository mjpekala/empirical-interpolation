% TEST_1d

% mjp, sept 2016

addpath('..')

n_max = 30;
delta = .02;
Omega = -1:delta:1;
W_n = make_polynomial_basis(1, n_max);

figure;
plot(Omega, W_n{1}(Omega), ...
     Omega, W_n{2}(Omega), ...
     Omega, W_n{3}(Omega), ...
     Omega, W_n{4}(Omega));
legend('x^0', 'x^1', 'x^2', 'x^3', 'Location', 'SouthEast');
title('first four polynomial basis functions (1d)');


% see figure 1 in [maday]
lambdas = zeros(n_max,1);
for m = 1:n_max
    [s, lambdas(m)] = choose_magic(Omega, W_n, m);
end

figure;
plot(1:n_max, lambdas, '-o');
xlabel('num. magic points');
ylabel('$\Lambda_M$')

figure;
subplot(2,1,1);
plot(s.Omega(s.x), ones(size(s.x)), 'o');
title('distribution of magic points');
xlabel('x');

subplot(2,1,2);
plot(acos(s.Omega(s.x)), ones(size(s.x)), 'o');
title('distribution of magic points');
xlabel('cos^{-1}(x)');

