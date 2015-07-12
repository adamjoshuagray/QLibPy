
from scipy.stats import norm
from scipy.optimize import brentq
from math import sqrt, log, exp

# TODO: Add comments

def _calc_d1(s, tau, k, r, sigma):
    d1  = 1 / (sigma * sqrt(tau)) * (log(s / k) + (r + pow(sigma, 2) / 2) * (tau))
    return d1

def _calc_d2(s, tau, k, r, sigma):
    d2  = _calc_d1(s, tau, k, r, sigma) - sigma * sqrt(tau)
    return d2

def calc_bs_call_px(s, tau, k, r, sigma):
    d_1     = _calc_d1(s, tau, k, r, sigma)
    d_2     = _calc_d2(s, tau, k, r, sigma)
    px      = s * norm.cdf(d_1) - k * exp(-r * tau) * norm.cdf(d_2)
    return px

def calc_bs_put_px(s, tau, k, r, sigma):
    px      = k * exp(-r * tau) - s + calc_bs_call_px(s, tau, k, r, sigma)
    return px

def calc_bs_call_sigma(c, s, tau, k, r):
    fn      = lambda x: calc_bs_call_px(s, tau, k, r, x) - c
    # Include all realistic levels for sigma
    a       = 0.0001
    b       = s * 10
    xtol    = 0.00000001 * c
    sigma   = brentq(fn, a, b, xtol = xtol)
    return sigma

def calc_bs_put_sigma(p, s, tau, k, r):
    fn      = lambda x: calc_bs_put_px(s, tau, k, r, x) - p
    # Include all realistic levels for
    a       = 0.0001
    b       = s * 10
    xtol    = 0.00000001 * c
    sigma   = brentq(fn, a, b, xtol = xtol)
    return sigma

def calc_bs_call_r(c, s, tau, k, sigma, px):
    fn      = lambda x: calc_bs_call_px(s, tau, k, x, )
    # We realisitically expect interest rates to be between -100% and 100%
    a       = -0.99
    b       = 1.0
    xtol    = 0.00000001 * c
    sigma   = brentq(fn, a, b, xtol = xtol)
    return sigma

def calc_bs_put_r(p, s, tau, k, sigma, px):
    fn      = lambda x: calc_bs_put_px(s, tau, k, x, )
    # We realisitically expect interest rates to be between -100% and 100%
    a       = -0.99
    b       = 1.0
    xtol    = 0.00000001 * c
    sigma   = brentq(fn, a, b, xtol = xtol)
    return sigma

def calc_bs_call_delta(s, tau, k, r, sigma):
    d_1     = _calc_d1(s, tau, k, r, sigma)
    delta   = norm.cdf(d_1)
    return delta

def calc_bs_put_delta(s, tau, k, r, sigma):
    delta   = calc_bs_call_delta(s, tau, k, r, sigma) - 1
    return delta

def calc_bs_gamma(s, tau, k, r, sigma):
    d_1     = _calc_d1(s, tau, k, r, sigma)
    gamma   = norm.pdf(d_1) / (s * sigma * sqrt(tau))
    return gamma

def calc_bs_vega(s, tau, k, r, sigma):
    d_1     = _calc_d1(s, tau, k, r, sigma)
    vega    = s * norm.pdf(d_1) * sqrt(tau)
    return vega

def calc_bs_call_theta(s, tau, k, r, sigma):
    d_1     = _calc_d1(s, tau, k, r, sigma)
    d_2     = _calc_d2(s, tau, k, r, sigma)
    theta   = -s * norm.pdf(d_1) * sigma / (2 * sqrt(tau)) - r * k * exp(-r * tau) * norm.cdf(d_2)
    return theta

def calc_bs_put_theta(s, tau, k, r, sigma):
    d_1     = _calc_d1(s, tau, k, r, sigma)
    d_2     = _calc_d2(s, tau, k, r, sigma)
    theta   = -s * norm.pdf(d_1) * sigma / (2 * sqrt(tau)) + r * k * exp(-r * tau) * norm.cdf(-d_2)
    return theta

def calc_bs_call_rho(s, tau, k, r, sigma):
    d_2     = _calc_d2(s, tau, k, r, sigma)
    rho     = k * tau * exp(-r * tau) * norm.cdf(d_2)
    return rho

def calc_bs_put_rho(s, tau, k, r, sigma):
    d_2     = _calc_d2(s, tau, k, r, sigma)
    rho     = -k * tau * exp(-r * tau) * norm.cdf(-d_2)
    return rho
