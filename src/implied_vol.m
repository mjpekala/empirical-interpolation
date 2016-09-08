function sigma = implied_vol(c, S0, K, r, t)
% IMPLIED_VOL  Computes the implied volatility for a European call.
%
%  Parameters:
%
%     c  : call option price at t0
%     S0 : stock price at t0
%     K  : strike price
%     r  : risk free interest rate
%     t  : option expiration time (relative to t0)
%
%  References:
%    Hull, "Options, Futures and Other Derivatives," fifth ed.


% see 12.20, 12.21 in [hull]
d1 = @(sigma) (log(S0/K) + (r + sigma^2/2)*t) / (sigma * sqrt(t));
d2 = @(sigma) d1(sigma) - sigma * sqrt(t);

f = @(sigma) S0*normcdf(d1(sigma)) - K*exp(-r*t)* normcdf(d2(sigma)) - c;

sigma = fzero(f, 0.2);


