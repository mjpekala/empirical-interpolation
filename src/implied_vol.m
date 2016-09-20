function sigma = implied_vol(c, S0, K, r, t)
% IMPLIED_VOL  Computes the implied volatility for a European call.
%
%  Parameters:
%    Note: parameters are either scalars or vectors; in latter
%          case, they must all have the same dimension.
%
%     c  : call option price at t0
%     S0 : stock price at t0
%     K  : strike price
%     r  : risk free interest rate
%     t  : option expiration time (relative to t0)
%
%  References:
%    Hull, "Options, Futures and Other Derivatives," fifth ed.

% mjp, sept 2016

n = max([numel(c) numel(S0) numel(K) numel(r) numel(t)]);
if n == 1
    % Scalar input case.
    f = @(x) f_sigma(x, c, S0, K, r, t);
    sigma = fzero(f, 0.2);
else
    % Vector input case.
    %
    % Note: I don't think fzero can solve multiple problems in parallel.
    %       Could parfor here, but probably not worth the overhead.
    sigma = zeros(n,1);

    for ii = 1:n
        f = @(x) f_sigma(x, ...
                         c(min(numel(c),ii)), ...
                         S0(min(numel(S0),ii)), ...
                         K(min(numel(K),ii)), ...
                         r(min(numel(r),ii)), ...
                         t(min(numel(t),ii)));
        if ii == 1
            sigma(ii) = fzero(f, 0.2);
        else
            sigma(ii) = fzero(f, sigma(ii-1));
        end
    end
end



function rv = f_sigma(sigma, c, S0, K, r, t)
% F_SIGMA  Function whose roots are implied volatilities.

% This works for scalar or vector inputs; however, it will always be
% called with scalar inputs for now.
d1 = (log(S0./K) + (r + sigma.^2/2).*t) ./ (sigma .* sqrt(t));
d2 = d1 - sigma .* sqrt(t);
rv = S0.*normcdf(d1) - K .* exp(-r.*t) .* normcdf(d2) - c;
