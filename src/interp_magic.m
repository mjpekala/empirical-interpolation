function [f_interp, s] = interp_magic(v, Omega, U_, m)
% INTERP_MAGIC   Interpolate a scalar-valued function using "magic points"
%
%   v     : The function to interpolate.  v : R^d -> R.
%           This code assumes v can handle vector valued inputs
%           (e.g. see apply.m)
%
%   Omega : a compact subset of R^d over which v will be modeled.
%           An (m x d) matrix of m points in d dimensions.
%   U_    : the function space to use when interpolating 
%           (a cell array containing at least m function handles).
%   m     : the number of magic points to use.
%
% References:
%   Maday et al. "A general multipurpose interpolation procedure:
%                 the magic points," CPAA 2009.


%% Check parameters
assert(0 < m);
assert(m <= length(U_));
assert(m <= size(Omega,1));

a = min(Omega(:));  b = max(Omega(:));
if (a < -1) || (b > 1)
    error('your problem should be rescaled to (a subset of) [-1,1]^2');
end


%% determine the q_i functions and the "magic points"

% The magic points interpolation data structure:
s.u = {};          % the selected elements of \mathcal{U}
s.q = {};          % the normalized interpolants
s.x = zeros(0,2);  % the "magic points" (interpolation points in \Omega)

% the first magic point is kind of a special case.
[s.u{1}, s.x(1,:)] = choose_next_magic_point(U_, Omega, s);
s.q{1} = @(X) s.u{1}(X) / s.u{1}(s.x(1,:));

% compute the remaining magic points and normalized functions
for k = 2:m
    fprintf('[%s]: choosing magic point %d (of %d)\n', mfilename, k, m);
    [u_k, x_k] = choose_next_magic_point(U_, Omega, s);
    I_k = make_interpolant(u_k, s);
    s.u{k} = u_k;
    s.x(k,:) = x_k;
    s.q{k} = @(X) bsxfun(@rdivide, u_k(X) - I_k(X), u_k(x_k) - I_k(x_k));
end


% Now that we know the magic points, can build a dedicated function
% for interpolating v.
f_interp = make_interpolant(v, s);



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
    I_k = make_interpolant(U_{k}, s);
    err_inf = abs(U_{k}(Omega) - I_k(Omega));
    [p_max, idx] = max(err_inf);
    if p_max > worst
        worst = p_max;
        u_i = U_{k};
        x_i = Omega(idx,:);
    end
end

assert(~any(isnan(x_i)));



function v_hat = make_interpolant(v, s)
% MAKE_INTERPOLANT  Returns a *function* v_hat that interpolates v at the points x.
%
%    v :  The function to interpolate;  v : R^d -> R
%    s :  The magic points interpolation structure 
%
    
m = length(s.q);
assert(size(s.x,1) == m);
assert(m > 0);

% Special case: for a single interpolant there is no need 
%               to compute any coefficients.
% See p. 386 in [maday]
if m == 1
    v_hat = @(X) v(s.x(1,:))*s.q{1}(X);
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
%
v_hat = @(X) sum(bsxfun(@times, b(:)', apply(s.q, X)), 2);
