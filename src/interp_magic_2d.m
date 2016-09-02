function [f_interp, s] = interp_magic_2d(v, X, Y, n_max)
% INTERP_MAGIC_2D   Interpolate a 2d function using "magic points"
%
%   v    : The function to interpolate.  v : R^2 -> R.
%          This code assumes v can handle vector valued inputs (see apply.m)
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
end


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

% define the basis functions 
% Currently uses the polynomial basis:
%
%   W_n(Omega) := { (x^i)*(y^j), (x,y) \in Omega, i+j <= n }, 0 <= n <= n_max
%

% TODO: make U_, Omega parameters so this code can be for any dimension.

U_ = cell(m,1);  % underscore is my notation indicating this is not a matrix
for kk = 1:m
    U_{kk} = @(X2d)  X2d(:,1).^ii(kk) .* X2d(:,2).^jj(kk);
end


%% determine the q_i functions and the "magic points"

% The magic points interpolation data structure:
s.u = {};          % the selected elements of \mathcal{U}
s.q = {};          % the normalized interpolants
s.x = zeros(0,2);  % the "magic points" (interpolation points in \Omega)

% the first point is selected in a special way (no interpolant yet)
Omega = [X(:) Y(:)];
tmp = abs(apply(u, Omega)); 
[~,idx] = max(tmp(:));
[r,c] = ind2sub(size(Omega), idx);
clear tmp;

M(1,:) = Omega(r,:);
q{1} = @(X) u{c}(X) ./ u{c}(M(1,:));
interp{1} = make_interpolant(u{c}, M(1,:), q{1});  % =: \mathcal{I}_i

keyboard % TEMP
for k = 2:m
end




function [u_i, x_i] = choose_next_magic_point(U_, Omega, s)
% CHOOSE_NEXT_MAGIC_POINT  Select the next magic point from Omega.
%
%    U_    := the function space (cell array of function handles)
%    Omega := the evaluation domain, an (m x d) of m points in d dimensions
%    s     := the magic points structure
 
% Special case: if there are not yet any magic points use 
%
%    arg max_{u \in U_}   ||u(.)||_\infty
%
% ie. evaluate all elements of U_ at all points in the domain
%     and pick the largest in magnitude.
if isempty(s.q)
    P = abs(apply(U_, Omega));
    [~,idx] = max(P(:));
    [r,c] = ind2sub(size(Omega), idx);
    x_i = Omega(r,:);
    u_i = U_{c};
    return
end

% Find the worst fitting point across all interpolants and use this as
% the next magic point.
worst = 0;  x_i = NaN;  u_i = NaN;
for k = 1:length(U_)
    interp = make_interpolant(U_{k}, s);
    p = interp(Omega);
    [p_max, idx] = max(p);
    if p_max > worst
        worst = p_max;
        u_i = U_{k};
        x_i = Omega(idx,:);
    end
end


function v_hat = make_interpolant(v, s)
% MAKE_INTERPOLANT  Returns a function v_hat that interpolates v at the points x.
%
%    v :  The function to interpolate;  v : R^d -> R
%    s :  The magic points interpolation structure 
%
    
m = length(s.q);
assert(m > 0);
assert(size(s.x,1) == m);

% Special case: for a single interpolant there is no need 
%               to compute any coefficients.
% See p. 386 in [maday]
if m == 1
    v_hat = @(X) v(X)*s.q{1}(X);
    return
end


% First must obtain the b coefficients (denoted beta in [maday]).
% Technically this does more work than needed since we throw away the
% upper triangular part; revisit this later if it turns out to be a
% performance issue.
Q = apply(s.q, s.x);
Q = tril(Q);

% The diagonal of Q should be all 1s; we force this here just in case there are
% any numerical issues.
Q(1:m+1:m*m) = 1;

% see (3) in [maday]
% This should be well-conditioned due to structure of Q.
b = Q \ v(s.x);

% see (2) in [maday]
v_hat = @(X) sum(bsxfun(@times, b(:)', apply(s.q, X)), 1);
