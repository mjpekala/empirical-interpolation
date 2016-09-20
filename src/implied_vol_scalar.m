function sigma = implied_vol_scalar(c, S0, K, r, t)
% IMPLIED_VOL_SCALAR  Determines a single implied volatility value.
%
%  Parameters (all are scalars):
%     c  : call option price at t0
%     S0 : stock price at t0
%     K  : strike price
%     r  : risk free interest rate
%     t  : option expiration time (relative to t0)
%
%  References:
%    Hull, "Options, Futures and Other Derivatives," fifth
%    ed. Chapter 12, equations 12.20, 12.21.

% TODO: a more general initial guess

f = @(x) raw_calc(x, c, S0, K, r, t);
sigma = fzero(f, 0.2);


function rv = raw_calc(sigma, c, S0, K, r, t)
  d1 = (log(S0/K) + (r + sigma^2/2)*t) / (sigma * sqrt(t));
  d2 = d1 - sigma * sqrt(t);
  rv = S0*normcdf(d1) - K*exp(-r*t)* normcdf(d2) - c;
