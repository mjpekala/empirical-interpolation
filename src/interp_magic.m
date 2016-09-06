function v_hat = interp_magic(v, s)

% evaluate v at the magic points
v_magic = v(s.Omega(s.x,:));

% extract the Q matrix
Q = s.Q_all(s.x, :);

% solve for beta
beta_M = Q \ v_magic;

% interpolate
v_hat = sum(bsxfun(@times, beta_M', s.Q_all), 2);
