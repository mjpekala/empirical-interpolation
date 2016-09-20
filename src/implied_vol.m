function sigma = implied_vol(c, S0, K, r, t)
% IMPLIED_VOL  Computes the implied volatility for a European call.
%
%  Parameters (either all scalar or all vectors of same length):
%     c  : call option price at t0
%     S0 : stock price at t0
%     K  : strike price
%     r  : risk free interest rate
%     t  : option expiration time (relative to t0)
%
%  References:
%    Hull, "Options, Futures and Other Derivatives," fifth ed.

% mjp, sept 2016


if length(c) == 1
    sigma = implied_vol_scalar(c, S0, K, r, t);
else
    sigma = zeros(length(c),1);
    for ii = 1:length(c)
        sigma(ii) = implied_vol_scalar(c(ii), S0(ii), K(ii), r(ii), t(ii));
    end
end
