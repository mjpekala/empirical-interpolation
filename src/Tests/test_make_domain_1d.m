% TEST_MAKE_DOMAIN_1d

addpath('..')

m_max = 30;
[Omega, W_n] = make_domain_1d(m_max,.01);

figure;
plot(Omega, W_n{1}(Omega), ...
     Omega, W_n{2}(Omega), ...
     Omega, W_n{3}(Omega), ...
     Omega, W_n{4}(Omega));
legend('x^0', 'x^1', 'x^2', 'x^3', 'Location', 'SouthEast');
title('first four polynomial basis functions (1d)');


% see figure 1 in [maday]
lambdas = zeros(m_max,1);
for m = 1:m_max
    [s, lambdas(m)] = choose_magic(Omega, W_n, m);
end

figure;
plot(1:m_max, lambdas, '-o');
xlabel('num. magic points');
ylabel('$\Lambda_M$')

figure;
plot(s.Omega(s.x), ones(size(s.x)), 'o');
title('distribution of magic points');
xlabel('x');

figure;
plot(acos(s.Omega(s.x)), ones(size(s.x)), 'o');
title('distribution of magic points');
xlabel('cos^{-1}(x)');

