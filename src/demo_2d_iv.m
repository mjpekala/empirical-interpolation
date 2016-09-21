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

if 0
    W_n = make_polynomial_basis(2, n_max);
    [Omega, domain_info] = make_domain_2d(delta, 'square');
    compare_methods_2d(iv_scaled, W_n, Omega, domain_info);
else
    [m_vals, err_ell_inf] = error_analysis_2d(iv_scaled, n_max, .01);
    figure;
    plot(m_vals, err_ell_inf(:,1), 'o-', m_vals, err_ell_inf(:,2), 'o-');
    legend('magic points', 'griddata');
    xlabel('m'); ylabel('err ell_inf');
    title('ell_infty error for implied vol example');
end

