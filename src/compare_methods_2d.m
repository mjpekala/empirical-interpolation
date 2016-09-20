function compare_methods_2d(v, W_n, Omega, domain_info)
% COMPARE_METHODS_2D Visual comparision of interpolation methods
%
%   This assumes a fixed domain/mesh.  A different analysis might
%   investigate error as a function of number of interpolation
%   points...


%% Brute force calculation
tic
v_true = v(Omega);
fprintf('[%s]: It took %0.2e seconds to evaluate v() at all points in Omega\n', mfilename, toc);

plot_mesh = @(v) mesh2(v, domain_info);

figure;
plot_mesh(v_true);
xlabel('x1'); ylabel('x2'); 
title('v_{true}');


%% Empirical interpolation
tic
[s, Lambda_M] = choose_magic(Omega, W_n);
fprintf('[%s]: It took %0.2e seconds to choose magic points\n', mfilename, toc);
fprintf('[%s]: Lambda_M = %0.3f\n', mfilename, Lambda_M);

tic
[v_hat, v_magic] = interp_magic(v, s);
fprintf('[%s]: It took %0.2e seconds for empirical interpolation\n', mfilename, toc);


%% Matlab's built-in interpolation
tic
v_hat_linear = griddata(s.Omega(s.x,1), s.Omega(s.x,2), v_magic, ...
                        s.Omega(:,1), s.Omega(:,2), ...
                        'linear');
fprintf('[%s]: It took %0.2e seconds to run griddata(linear)\n', mfilename, toc);

tic
v_hat_cubic = griddata(s.Omega(s.x,1), s.Omega(s.x,2), v_magic, ...
                       s.Omega(:,1), s.Omega(:,2), ...
                       'cubic');
fprintf('[%s]: It took %0.2e seconds to run griddata(cubic)\n', mfilename, toc);


%% TODO: GPML (Kriging)


%% Chebfun 2d
%
% Note: implicitly assumes domain is [-1,1]^2
%
if exist('chebfun2') && false
    % chebfun2 will provide data *matrices*; requires some reshaping.
    v_chebfun = @(X,Y) reshape(v([X(:) Y(:)]), size(X,1), size(Y,1));
    tic
    cf2 = chebfun2(v_chebfun);
    fprintf('[%s]: It took %0.2e seconds to run chebfun2()\n', mfilename, toc);
end



%% Visualize fits

err_ell1 = @(A,B) abs(A-B);
err_ell_inf = @(A,B) max(abs(A(:)-B(:)));

Err_magic_ell1 = err_ell1(v_true, v_hat);
Err_cubic_ell1 = err_ell1(v_true, v_hat_cubic);
zmax = max([ Err_magic_ell1(:) ; Err_cubic_ell1(:)]);

err_magic_ell_inf = err_ell_inf(v_true, v_hat);
err_cubic_ell_inf = err_ell_inf(v_true, v_hat_cubic);

% show magic points
figure;
scatter(s.Omega(s.x,1), s.Omega(s.x,2), 'ro');
xlabel('x1'); ylabel('x2');
title('interpolation points')


% empirical interpolation
figure;
plot_mesh(v_hat);
xlabel('x1'); ylabel('x2'); 
title('$\hat{v}_{magic}$', 'interpreter', 'latex');
hold on;
stem3(s.Omega(s.x,1), s.Omega(s.x,2), v_magic, 'ro');
hold off;


figure;
plot_mesh(Err_magic_ell1);
caxis([0 zmax]);
xlabel('x1'); ylabel('x2');  zlabel('l_1 err');
title(sprintf('err (magic points); ell_inf=%0.2e', err_magic_ell_inf));

%linkprop([ax1 ax2], 'CameraPosition');


% triangulation + cubic interpolation

figure;
plot_mesh(v_hat_cubic);
xlabel('x1'); ylabel('x2'); 
title('$\hat{v}_{cubic}$', 'interpreter', 'latex');
hold on;
stem3(s.Omega(s.x,1), s.Omega(s.x,2), v_magic, 'ro');
hold off;

figure;
plot_mesh(Err_cubic_ell1);
caxis([0 zmax]);
xlabel('x1'); ylabel('x2'); 
title(sprintf('err (cubic); ell_inf=%0.2e', err_cubic_ell_inf));
