
from math import *
from scipy.optimize import brentq

# Calculates the present value of a series of amoritzing payments
# and a principle given a certain interest rate.
# fv    - the future value of the instrument
# i     - the interest rate for each period (this may be an array)
# n     - the number of periods
# pmt   - the amortized payment made every period (this may also be an array)
def calc_pv(fv, i, n, pmts):
    if type(pmts) is float:
        pmts    = [pmts] * floor(n)
    if type(i) is float:
        i       = [i] * floor(n)
    exponents   = range(1, int(ceil(n + 1)))
    fact_pairs  = zip(exponents, i)
    factors     = map(lambda (e, j): 1 / pow(1 + j, e), fact_pairs)
    pmt_pairs   = zip(factors, pmts)
    pv_cfs      = map(lambda (f, p): p * f, pmt_pairs)
    pv_pr       = fv * 1 / pow(1 + i[-1], n)
    pv          = sum([pv_pr] + pv_cfs)
    return pv

# Calculates the future value of a series of amorizing payments
# given a certain interest rate
# pv    - the present value of the instrument
# i     - the interest rate for each period (this may be an array)
# n     - the number of periods
# pmt   - the amortized payment made every period (this may also be an array)
def calc_fv(pv, i, n, pmts):
    if type(pmts) is float:
        pmts    = [pmts] * floor(n)
    if type(i) is float:
        i       =[i] * floor(n)
    exponents   = range(1, int(ceil(n + 1)))
    fact_pairs  = zip(exponents, i)
    factors     = map(lambda (e, j): 1 / pow(1 + j, e), fact_pairs)
    pmt_pairs   = zip(factors, pmts)
    pv_cfs      = map(lambda (f, p): -p * f, pmt_pairs)
    fv          = sum([pv] + pv_cfs) * pow(1 + i[-1], n)
    return fv

# Calculates the payment which needs to be made each period
# to give the future value given the present value, interest rate and
# number of periods
# pv    - the present value
# fv    - the future value
# i     - the interest rate for each period
# n     - the number of periods
def calc_pmt(pv, fv, i, n):
    unit_pmts   = calc_pv(0, i, n, 1)
    pmt         = (pv - fv * pow(1 + i, -n)) / unit_pmts
    return pmt

# Calculates the yeild to maturity of a sequence of payments given '
# a present value, future value, number and size of payments
# pv    - the present value
# fv    - the future value
# n     - the number of periods
# pmt   - the amortized payment made each period
def calc_i(pv, fv, n, pmt):
    fn      = lambda x: calc_pv(fv, x, n, pmt) - pv
    # We realistically assume interest rates between - 100% and 100% per period.
    a       = -0.99
    b       = 1.0
    # We want to be within e-3 of pv
    xtol    = 0.001
    i       = brentq(fn, a, b, xtol = xtol)
    return i

# Calculates the number of payments which need to be made in order
# to obtain the future value, give the present value, interest rate and
# payment size per period
# pv    - the present value
# fv    - the future value
# i     - the interest rate per period
# pmt   - the amortized payment made each period
def calc_n(pv, fv, i, pmt):
    fn      = lambda x: calc_pv(fv, i, x, pmt) - pv
    # We can never have a negative number of periods.
    a       = 0
    # After 200 periods additional periods really have no effect.
    b       = 200
    # We want to be within e-3 of pv
    xtol    = 0.003
    n       = brentq(fn, a, b, xtol = xtol)
    return n

def calc_dv01(fv, i, n, pmt):
    pv1     = calc_pv(fv, i + 0.0001, n, pmt)
    pv2     = calc_pv(fv, i - 0.0001, n, pmt)
    dv01    = (pv2 - pv1) / 2.0
    return dv01

def calc_convex(fv, i, n, pmt):
    dv1     = calc_dv01(fv, i + 0.0001, n, pmt)
    dv2     = calc_dv01(fv, i - 0.0001, n, pmt)
    convex  = (dv2 - dv1) / 2.0
    return convex
