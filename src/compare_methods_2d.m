function compare_methods_2d(v, W_n, Omega, domain_info)
% COMPARE_METHODS_2D Visual comparision of interpolation methods


%% Brute force calculation
tic
v_true = v(Omega);
t = toc;
fprintf('[%s]: It took %0.2e seconds to evaluate v() at all points in Omega\n', mfilename, t);

plot_mesh = @(v) mesh2(v, domain_info);

figure;
plot_mesh(v_true);
xlabel('x1'); ylabel('x2'); 
title('v_{true}');


%% Empirical interpolation
tic
[s, Lambda_M] = choose_magic(Omega, W_n);
t = toc;
fprintf('[%s]: It took %0.2e seconds to choose magic points\n', mfilename, t);
fprintf('[%s]: Lambda_M = %0.3f\n', mfilename, Lambda_M);

tic
[v_hat, v_magic] = interp_magic(v, s);
t = toc;
fprintf('[%s]: It took %0.2e seconds for empirical interpolation\n', mfilename, t);


%% Matlab's built-in interpolation
tic
v_hat_linear = griddata(s.Omega(s.x,1), s.Omega(s.x,2), v_magic, ...
                        s.Omega(:,1), s.Omega(:,2), ...
                        'linear');
t = toc;
fprintf('[%s]: It took %0.2e seconds to run griddata(linear)\n', mfilename, t);

tic
v_hat_cubic = griddata(s.Omega(s.x,1), s.Omega(s.x,2), v_magic, ...
                       s.Omega(:,1), s.Omega(:,2), ...
                       'cubic');
t = toc;
fprintf('[%s]: It took %0.2e seconds to run griddata(cubic)\n', mfilename, t);


%% TODO: GPML (Kriging)


%% Visualize fits

calc_err = @(A,B) abs(A-B);

Err_magic = calc_err(v_true, v_hat);
Err_cubic = calc_err(v_true, v_hat_cubic);
err_max = max([ Err_magic(:) ; Err_cubic(:)]);

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
plot_mesh(Err_magic);
caxis([0 err_max]);
xlabel('x1'); ylabel('x2');
title(sprintf('err (magic points): %0.2e', sum(Err_magic(:))));

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
plot_mesh(Err_cubic);
caxis([0 err_max]);
xlabel('x1'); ylabel('x2'); 
title(sprintf('err (cubic): %0.2e', sum(Err_cubic(:))));
