function v_hat = interp_magic(v, s)
% INTERP_MAGIC  Interpolate a function v.
%
%   v : the function to interpolate
%   s : a "magic points" structure, as created by choose_magic()
%
%   v_hat : the result of evaluating v at the magic points and then
%           interpolating for all other points in s.Omega.
%
% References:
%   Maday et al. "A general multipurpose interpolation procedure:
%                 the magic points," CPAA 2009.
%


% evaluate v at the magic points
v_magic = v(s.Omega(s.x,:));

% extract the Q matrix
% Note: this is called B^M in [maday]
Q = s.Q_all(s.x, :);

% solve for beta
beta_M = Q \ v_magic;

% interpolate
v_hat = sum(bsxfun(@times, beta_M', s.Q_all), 2);

