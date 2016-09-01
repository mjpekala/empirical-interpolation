function f_interp = interp_magic_2d(f, X, Y, n_max)
% INTERP_MAGIC_2D   Interpolate a 2d function using "magic points"
%
%   X, Y : defines the evaluation domain \bar{Omega}.
%          Specifically, these define a gridding of (a subset of) [-1,1]^2
%          where X, Y are as obtained from meshgrid.
%
% References:
%   Maday et al. "A general multipurpose interpolation procedure:
%                 the magic points," CPAA 2009.
    
% Example:
%{
    % construct a toy function
    f = @(x,y) cos(3*x) .* sin(2*y);
    x = -pi:.1:pi;
    f_rescaled = @(x,y) f(x*pi,y*pi);
    
    [X,Y] = meshgrid(x / pi, x / pi);
    Z = f_rescaled(X,Y);
    surf(X,Y,Z);  xlabel('x'); ylabel('y');
    
    f_interp = interp_magic_2d(f_rescaled, X, Y, 10);
%}

a = min(min(X(:)), min(Y(:)));
b = max(max(X(:)), max(Y(:)));
if (a < -1) || (b > 1)
    error('your problem should be rescaled to [-1,1]^2 before calling');

if (min(X(:)) < -1) || (

% m := number of ("magic") interpolation points
m = (n_max+1)*(n_max+2) / 2;   


%% define basis functions u_i

% the code below is just a quick way to enumerate the space:
%
%   ii+jj <= n      for 0 <= n <= n_max
%
Tmp = tril(ones(n_max+1,n_max+1));
[ii,jj] = ind2sub(size(Tmp), find(Tmp));
jj = jj - 1;
ii = n_max+1-ii;
assert(length(ii) == m);
clear Tmp;

% define the basis functions and pre-calculate the evaluation of
% each at all points in the domain Omega.
u_all = cell(m,1);
Z = zeros(size(X,1), size(X,2), m);
parfor k = 1:m
    u_all{k} = @(x,y) x.^ii(k) .* y.^jj(k);
    Z(:,:,k) = u_all{k}(X,Y);
end

%% compute the q_i
q = zeros(m,1);

% the first point is selected in a special way (no interpolant yet)
[~,idx] = max(abs(Z(:)))
[r,c,z] = ind2sub(size(Z), idx);

keyboard % TEMP

    
function z = evaluate(f, x1, x2, q, beta)
%  
%   I_M[v(.)] = \sum_{j=1}^M beta_j q_j(.)
%
%     q_j(x) : R^2 -> R
%     v(x)   : R^2 -> R
%
% parfor over the basis function evaluations
    
    x1(:), x2(:)
    z = [] % TODO
    