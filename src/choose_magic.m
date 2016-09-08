function [s, lambda_M] = choose_magic(Omega, U_, m)
% CHOOSE_MAGIC  Chooses "magic points" for interpolation [maday]
%
%   Omega : Interpolation domain; a compact subset of R^d.
%           An (n x d) matrix of n points in d dimensions.
%   U_    : the function space to use when interpolating 
%           (a cell array containing at least m function handles).
%   m     : the number of magic points to use (optional).
%
% References:
%   Maday et al. "A general multipurpose interpolation procedure:
%                 the magic points," CPAA 2009.

% mjp, sept 2016

if nargin < 3, m = length(U_); end

if size(Omega,1) == 1, Omega = Omega(:); end

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
Q = zeros(m);                 % Q_{i,j} := q_j(x_i)

% 
U_all = apply(U_, Omega);     % each u_i evaluated at all points in Omega; (n x b)
Q_all = zeros(n, m);          % rescaled active basis functions



%% determine the q_i functions and the interpolation points

% I.  Choose the first interpolation point
[s.sf(1), s.x(1), s.u(1)] = ell_inf_2d(U_all);
Q_all(:,1) = U_all(:,s.u(1)) / s.sf(1);
Q(1,1) = 1;


% I_j := Current interpolant of all basis functions 
%        (using j basis functions).
beta1 = U_all(s.x(1),:) ./ Q_all(s.x(1),1);
I_j = bsxfun(@times, beta1, Q_all(:,1)); 

% the active basis functions should interpolate themselves perfectly
err = abs(U_all(:,s.u(1)) - I_j(:,s.u(1)));
assert(max(err) < 1e-8);

% II.  Choose the remaining m-1 points
for jj = 2:m
    if mod(jj,10) == 0
        fprintf('[%s]: choosing basis element %d (of %d)\n', mfilename, jj, m);
    end
    
    % Choose next basis function and interpolation point
    [s.sf(jj), s.x(jj), s.u(jj)] = ell_inf_2d(U_all - I_j);
    Q_all(:,jj) = (U_all(:,s.u(jj)) - I_j(:,s.u(jj))) / s.sf(jj); % rescale max err to 1

    Q = Q_all(s.x(1:jj), 1:jj);
    
    % sanity checks
    err = max(max(abs(triu(Q,1))));
    if err > 1e-9
        warning(sprintf('[%s]: upper triangular part of Q had nonzero value with mag %0.2e\n', mfilename, err));
        Q = tril(Q);
    end
    
    assert(abs(Q(jj,jj) - 1) < 1e-8);
    assert(length(unique(s.u(1:jj))) == jj);
 
    % update coefficients for each basis function.
    % TODO: rank one update to linear system?
    U_at_magic = U_all(s.x(1:jj), :);
    Q_jj = Q_all(:,1:jj);
    for kk = 1:size(U_all,2)
        beta_j_k = Q \ U_at_magic(:,kk);
        I_j(:,kk) = sum(bsxfun(@times, beta_j_k', Q_jj), 2);
    end
    
    % the active basis functions should interpolate themselves perfectly
    Err = abs(U_all(:,s.u(1:jj)) - I_j(:,s.u(1:jj)));
    assert(max(Err(:)) < 1e-9);
end

% add some data to s before returning
s.U_ = U_;
s.Omega = Omega;
s.Q_all = Q_all;


% III. approximate lebesgue coefficient (optional)
%      see p.387 in [maday] for definitions.
if nargout > 1
    Q_inv = inv(Q); 

    H = zeros(size(s.Q_all,1), length(s.x));
    for ii=1:length(s.x)
        H(:,ii) = s.Q_all * Q_inv(:,ii);
        
        % it should be the case that h_i(x_j) = \delta_ij
        deltaij = zeros(size(s.x));  deltaij(ii) = 1;
        assert(all(abs(deltaij - H(s.x,ii))) < 1e-8);
    end
    lambda_M = max(sum(abs(H), 2)); 
end


function [val, row, col] = ell_inf_2d(A)
% MAX_IN_MATRIX  Finds the maximum magnitude value in a 2d matrix.
%
%  This is just a bit of shorthand for clarity.
[~, idx] = max(abs(A(:)));
[row,col] = ind2sub(size(A), idx);
val = A(row,col);
