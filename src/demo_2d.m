% DEMO_2D   Simple example of "magic point" interpolation
%
% References:
%   Maday et al. "A general multipurpose interpolation procedure:
%                 the magic points," CPAA 2009.

% mjp, sept 2016


%% construct a toy function on a triangular domain
f = @(X) cos(3*X(:,1)) .* sin(2*X(:,2));
f_rescaled = @(X) f(X*pi);

n_max = 14;

if 0
    W_n = make_polynomial_basis(2, n_max);
    [Omega, domain_info] = make_domain_2d(.02, 'triangle');
    compare_methods_2d(f_rescaled, W_n, Omega, domain_info);
else
    [m_vals, err_ell_inf] = error_analysis_2d(iv_scaled, n_max, .01);
    figure;
    plot(m_vals, err_ell_inf(:,1), 'o-', m_vals, err_ell_inf(:,2), 'o-');
    legend('magic points', 'griddata');
    xlabel('m'); ylabel('err ell_inf');
    title('ell_infty error for implied vol example');
end

