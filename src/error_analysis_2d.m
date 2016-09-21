function [m_vals, err_ell_inf] = error_analysis_2d(v, n_max, delta)
% assumes the domain is [-1 1]^2 and a polynomial basis.

if nargin < 3, delta = 0.001; end
if nargin < 2, n_max = 20; end 

W_n = make_polynomial_basis(2, n_max);
[Omega, domain_info] = make_domain_2d(delta, 'square');

fprintf('[%s]: Evaluating true function on finest grid...\n', mfilename);
tic
v_true = v(Omega);
fprintf('[%s]: It took %0.2e seconds to evaluate v() at all points in Omega\n', mfilename, toc);


m_vals = (1:5:length(W_n))';
err_ell_inf = zeros(numel(m_vals),2);

for ii = 1:length(m_vals), m=m_vals(ii);
    % compute error when using m magic points
    [s, Lambda_M] = choose_magic(Omega, W_n, m);
    
    tic
    [v_hat, v_magic] = interp_magic(v, s);
    fprintf('[%s]: It took %0.2e seconds for empirical interpolation\n', mfilename, toc);

    err_ell_inf(ii,1) = max(abs(v_hat(:) - v_true(:)));
    
    % For comparison: matlab's 2d interp. using magic points
    v_hat_cubic = griddata(s.Omega(s.x,1), s.Omega(s.x,2), v_magic, ...
                           s.Omega(:,1), s.Omega(:,2), ...
                           'cubic');
    if isempty(v_hat_cubic)
        err_ell_inf(ii,2) = max(abs(v_true(:)));
    else
        err_ell_inf(ii,2) = max(abs(v_hat_cubic(:) - v_true(:)));
    end
end
