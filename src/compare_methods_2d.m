function compare_methods_2d(v, W_n, Omega, domain_info)
% COMPARE_METHODS_2D Visual comparision of interpolation methods


%% Brute force calculation
tic
v_true = v(Omega);
t = toc;
fprintf('[%s]: It took %0.2e seconds to evaluate v at all points in Omega\n', mfilename, t);

plot_mesh = @(v) mesh2(v, domain_info);

figure;
plot_mesh(v_true);
xlabel('x1'); ylabel('x2'); 
title('f_{true}');


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


%% TODO: GPML


%% Visualize fits

figure('Position', [100 100 1200 450]); 

ax1 = subplot(1,3,1);
plot_mesh(v_hat);
xlabel('c'); ylabel('k'); 
title('$\hat{\sigma}$', 'interpreter', 'latex');
hold on;
stem3(s.Omega(s.x,1), s.Omega(s.x,2), v_magic, 'ro');
hold off;

ax2 = subplot(1,3,2); 
plot_mesh(v_hat_linear);
xlabel('c'); ylabel('k'); 
title('$\hat{\sigma}_{linear}$', 'interpreter', 'latex');
hold on;
stem3(s.Omega(s.x,1), s.Omega(s.x,2), v_magic, 'ro');
hold off;

ax3 = subplot(1,3,3);
plot_mesh(v_hat_cubic);
xlabel('c'); ylabel('k'); 
title('$\hat{\sigma}_{cubic}$', 'interpreter', 'latex');
hold on;
stem3(s.Omega(s.x,1), s.Omega(s.x,2), v_magic, 'ro');
hold off;

linkprop([ax1 ax2 ax3], 'CameraPosition');


%% Error analysis

figure('Position', [100 100 1200 450]); 

ax4 = subplot(1,3,1);
plot_mesh(abs(v_true - v_hat));
xlabel('c'); ylabel('k'); 
title('err (magic points)');

ax5 = subplot(1,3,2);
plot_mesh(abs(v_true - v_hat_linear));
xlabel('c'); ylabel('k'); 
title('err (linear)');

ax6 = subplot(1,3,3);
plot_mesh(abs(v_true - v_hat_cubic));
xlabel('c'); ylabel('k'); 
title('err (cubic)');

%linkprop([ax4 ax5 ax6], 'CameraPosition');
