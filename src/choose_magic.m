function [f_interp, s] = interp_magic(Omega, U_, m)
% INTERP_MAGIC  Chooses "magic points" for interpolation [maday]
%
%   Omega : Interpolation domain; a compact subset of R^d.
%           An (n x d) matrix of n points in d dimensions.
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
if (min(Omega(:)) < -1) || (1 < max(Omega(:)))
    error('your problem should be rescaled to (a subset of) [-1,1]^2');
end


%% Variables
[n,d] = size(Omega);

% s := the empirical interpolation data structure:
%
s.u = zeros(m,1);             % active basis functions; an index into U_
s.x = zeros(m,1);             % interpolation points, encoded as indices into Omega
s.sf = zeros(m,1);            % scale factors for the q_i
s.Q = eye(m);                 % Q_{i,j} := q_j(x_i)

% 
U_all = apply(U_, Omega);     % each u_i evaluated at all points in Omega; (n x b)
B_all = zeros(m, length(U_)); % the beta coefficients for all u in U_

% 
Q_all = zeros(n, m);          % the ith column of Q_all is q_i for all points in Omega

keyboard % TEMP

%% determine the q_i functions and the interpolation points

% first point is a kind of special case
[maxval,idx] = max(abs(U_all));
[r,c] = ind2sub(size(Z), idx);
s.u(1) = c;
s.x(1) = r; 
s.sf(1) = maxval;
Q_all(:,1) = U_all(:,s.u(1)) / s.sf(1);

B_all(1,:) = s.v

for jj = 1:m
    
end


return % TEMP




if 0
% the first magic point is kind of a special case.
[s.u{1}, s.x(1,:)] = choose_next_magic_point(U_, Omega, s);
scale = 1.0 / s.u{1}(s.x(1,:));
s.q{1} = @(X) s.u{1}(X) / scale;

% compute the remaining magic points and normalized functions
for k = 2:m
    fprintf('[%s]: choosing magic point %d (of %d)\n', mfilename, k, m);
    [u_k, x_k, worst] = choose_next_magic_point(U_, Omega, s);
    
    assert(~s.basis(worst));
    s.basis(worst) = true;
    
    I_k = make_interpolant(u_k, s);
    s.u{k} = u_k;
    s.x(k,:) = x_k;
    scale = 1.0 / (u_k(x_k) - I_k(x_k));
    s.q{k} = @(X) scale*(u_k(X) - I_k(X));
end


% Now that we know the magic points, can build a dedicated function
% for interpolating v.
f_interp = make_interpolant(v, s);
end


function I_k = calc_I_k_all(U_all, s)
% CALC_I_k_ALL  Computes I_k[u] for all u
%
k_max = find(s.u == 0, 1, 'first') - 1;
assert(k_max > 1);

I_k = zeros(size(U_all));

q_1 = U_all(:,s.u(1)) / s.sf(1);
Q_k




function [u_i, x_i, k_worst] = choose_next_magic_point(U_, Omega, s)
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
    [r,k_worst] = ind2sub(size(Omega), idx);
    x_i = Omega(r,:);
    u_i = U_{k_worst};
    return
end

% Find the worst fitting point across all interpolants and use this as
% the next magic point.
worst = 0;  k_worst = NaN;  x_i = NaN;  u_i = NaN;
for k = 1:length(U_)
    I_k = make_interpolant(U_{k}, s);
    err_inf = abs(U_{k}(Omega) - I_k(Omega));
    [p_max, idx] = max(err_inf);
    if p_max > worst
        worst = p_max;
        k_worst = k;
        u_i = U_{k};
        x_i = Omega(idx,:);
    end
end

assert(~isnan(k_worst));



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

% These next two statements are not necessary; are entailed by algorithm.
%Q = tril(Q);
%Q(1:m+1:m*m) = 1;

% see (3) in [maday]
% This should be well-conditioned due to structure of Q.
b = Q \ v(s.x);

% see (2) in [maday]
%
v_hat = @(X) sum(bsxfun(@times, b(:)', apply(s.q, X)), 2);
