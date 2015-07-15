
from scipy.stats import norm
from scipy.optimize import brentq
from math import sqrt, log, exp

# TODO: Add comments
# TODO: Add dividend yeilds

# An internal method which calculates d_1 which is commonly used in calculations
def _calc_d1(s, tau, k, r, sigma):
    d1  = 1 / (sigma * sqrt(tau)) * (log(s / k) + (r + pow(sigma, 2) / 2) * (tau))
    return d1

# An internal method which calculates d_2 which is commonly used in calculations
def _calc_d2(s, tau, k, r, sigma):
    d2  = _calc_d1(s, tau, k, r, sigma) - sigma * sqrt(tau)
    return d2

# Calculates the price of a standard European call option
# s     The spot price of the underlying instrument.
# tau   The time to expiry of the option
# k     The strike price of the option
# r     The continuously compounded anualized risk free rate
# sigma The anualized standard deviation of returns.
def calc_call_px(s, tau, k, r, sigma):
    d_1     = _calc_d1(s, tau, k, r, sigma)
    d_2     = _calc_d2(s, tau, k, r, sigma)
    px      = s * norm.cdf(d_1) - k * exp(-r * tau) * norm.cdf(d_2)
    return px
# Calculates the price of a standard European put option
# s     The spot price of the underlying instrument.
# tau   The time to expiry of the option
# k     The strike price of the option
# r     The continuously compounded anualized risk free rate
# sigma The anualized standard deviation of returns.
def calc_put_px(s, tau, k, r, sigma):
    px      = k * exp(-r * tau) - s + calc_call_px(s, tau, k, r, sigma)
    return px

# Calculates the implied volatility for a standard European call option.
# s     The spot price of the underlying instrument.
# tau   The time to expiry of the option
# k     The strike price of the option
# r     The continuously compounded anualized risk free rate
# sigma The anualized standard deviation of returns.
def calc_call_sigma(c, s, tau, k, r):
    fn      = lambda x: calc_call_px(s, tau, k, r, x) - c
    # Include all realistic levels for sigma
    a       = 0.0001
    b       = s * 10
    xtol    = 0.00000001 * c
    sigma   = brentq(fn, a, b, xtol = xtol)
    return sigma

# Calculates the implied volatility for a standard European put option.
# s     The spot price of the underlying instrument.
# tau   The time to expiry of the option
# k     The strike price of the option
# r     The continuously compounded anualized risk free rate
# sigma The anualized standard deviation of returns.
def calc_put_sigma(p, s, tau, k, r):
    fn      = lambda x: calc_put_px(s, tau, k, r, x) - p
    # Include all realistic levels for
    a       = 0.0001
    b       = s * 10
    xtol    = 0.00000001 * c
    sigma   = brentq(fn, a, b, xtol = xtol)
    return sigma

# Calculates the implied risk free rate for a standard European call option.
# s     The spot price of the underlying instrument.
# tau   The time to expiry of the option
# k     The strike price of the option
# r     The continuously compounded anualized risk free rate
# sigma The anualized standard deviation of returns.
def calc_call_r(c, s, tau, k, sigma, px):
    fn      = lambda x: calc_call_px(s, tau, k, x, sigma) - p
    # We realisitically expect interest rates to be between -100% and 100%
    a       = -0.99
    b       = 1.0
    xtol    = 0.00000001 * c
    sigma   = brentq(fn, a, b, xtol = xtol)
    return sigma

# Calculates the implied risk free rate for a standard European put option.
# s     The spot price of the underlying instrument.
# tau   The time to expiry of the option
# k     The strike price of the option
# r     The continuously compounded anualized risk free rate
# sigma The anualized standard deviation of returns.
def calc_put_r(p, s, tau, k, sigma, px):
    fn      = lambda x: calc_put_px(s, tau, k, x, sigma) - p
    # We realisitically expect interest rates to be between -100% and 100%
    a       = -0.99
    b       = 1.0
    xtol    = 0.00000001 * c
    sigma   = brentq(fn, a, b, xtol = xtol)
    return sigma

# Calculates the delta of a standard European call option.
# s     The spot price of the underlying instrument.
# tau   The time to expiry of the option
# k     The strike price of the option
# r     The continuously compounded anualized risk free rate
# sigma The anualized standard deviation of returns.
def calc_call_delta(s, tau, k, r, sigma):
    d_1     = _calc_d1(s, tau, k, r, sigma)
    delta   = norm.cdf(d_1)
    return delta

# Calculates the delta of a standard European put option.
# s     The spot price of the underlying instrument.
# tau   The time to expiry of the option
# k     The strike price of the option
# r     The continuously compounded anualized risk free rate
# sigma The anualized standard deviation of returns.
def calc_put_delta(s, tau, k, r, sigma):
    delta   = calc_call_delta(s, tau, k, r, sigma) - 1
    return delta

# Calculates the gamma of a standard European option.
# Note that this is the same for puts and calls.
# s     The spot price of the underlying instrument.
# tau   The time to expiry of the option
# k     The strike price of the option
# r     The continuously compounded anualized risk free rate
# sigma The anualized standard deviation of returns.
def calc_gamma(s, tau, k, r, sigma):
    d_1     = _calc_d1(s, tau, k, r, sigma)
    gamma   = norm.pdf(d_1) / (s * sigma * sqrt(tau))
    return gamma

# Calculates the vega of a standard European option.
# Note that this is the same for puts and calls.
# s     The spot price of the underlying instrument.
# tau   The time to expiry of the option
# k     The strike price of the option
# r     The continuously compounded anualized risk free rate
# sigma The anualized standard deviation of returns.
def calc_vega(s, tau, k, r, sigma):
    d_1     = _calc_d1(s, tau, k, r, sigma)
    vega    = s * norm.pdf(d_1) * sqrt(tau)
    return vega

# Calculates the theta of a standard European call option.
# Note that this is the same for puts and calls.
# s     The spot price of the underlying instrument.
# tau   The time to expiry of the option
# k     The strike price of the option
# r     The continuously compounded anualized risk free rate
# sigma The anualized standard deviation of returns.
def calc_call_theta(s, tau, k, r, sigma):
    d_1     = _calc_d1(s, tau, k, r, sigma)
    d_2     = _calc_d2(s, tau, k, r, sigma)
    theta   = -s * norm.pdf(d_1) * sigma / (2 * sqrt(tau)) - r * k * exp(-r * tau) * norm.cdf(d_2)
    return theta
# Calculates the theta of a standard European put option.
# Note that this is the same for puts and calls.
# s     The spot price of the underlying instrument.
# tau   The time to expiry of the option
# k     The strike price of the option
# r     The continuously compounded anualized risk free rate
# sigma The anualized standard deviation of returns.
def calc_put_theta(s, tau, k, r, sigma):
    d_1     = _calc_d1(s, tau, k, r, sigma)
    d_2     = _calc_d2(s, tau, k, r, sigma)
    theta   = -s * norm.pdf(d_1) * sigma / (2 * sqrt(tau)) + r * k * exp(-r * tau) * norm.cdf(-d_2)
    return theta

# Calculates the rho of a standard European call option.
# Note that this is the same for puts and calls.
# s     The spot price of the underlying instrument.
# tau   The time to expiry of the option
# k     The strike price of the option
# r     The continuously compounded anualized risk free rate
# sigma The anualized standard deviation of returns.
def calc_call_rho(s, tau, k, r, sigma):
    d_2     = _calc_d2(s, tau, k, r, sigma)
    rho     = k * tau * exp(-r * tau) * norm.cdf(d_2)
    return rho

# Calculates the rho of a standard European put option.
# Note that this is the same for puts and calls.
# s     The spot price of the underlying instrument.
# tau   The time to expiry of the option
# k     The strike price of the option
# r     The continuously compounded anualized risk free rate
# sigma The anualized standard deviation of returns.
def calc_put_rho(s, tau, k, r, sigma):
    d_2     = _calc_d2(s, tau, k, r, sigma)
    rho     = -k * tau * exp(-r * tau) * norm.cdf(-d_2)
    return rho

# Calculates the vanna of a standard European option.
# Note that this is the same for puts and calls.
# s     The spot price of the underlying instrument.
# tau   The time to expiry of the option
# k     The strike price of the option
# r     The continuously compounded anualized risk free rate
# sigma The anualized standard deviation of returns.
def calc_vanna(s, tau, k, r, sigma):
    d_1     = _calc_d1(s, tau, k, r, sigma)
    vega    = calc_vega(s, tau, k, r, sigma)
    vanna   = vega / s * (1 - d_1 / (sigma * sqrt(tau)))
    return vanna

# Calculates the charm of a standard European call option.
# Note that this is the same for puts and calls.
# s     The spot price of the underlying instrument.
# tau   The time to expiry of the option
# k     The strike price of the option
# r     The continuously compounded anualized risk free rate
# sigma The anualized standard deviation of returns.
def calc_call_charm(s, tau, k, r, sigma, q = 0):
    d_1     = _calc_d1(s, tau, k, r, sigma)
    d_2     = _calc_d2(s, tau, k, r, sigma)
    a       = q * exp(-q * tau) * norm.cdf(d_1)
    b       = exp(-q * r) * norm.pdf(d_1)
    c       = (2 * (r - q) * tau - d_2 * sigma * sqrt(tau)) / (2 * tau * sigma * sqrt(tau))
    charm   = a - b * c

# TODO: Implement
def calc_put_charm(s, tau, k, r, sigma, q = 0):
    pass
