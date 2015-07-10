
from math import *

def calc_pv(fv, i, n, pmt):
    exponents   = range(1, int(ceil(n + 1)))
    factors     = map(lambda k: 1 / pow(1 + i, k), exponents)
    pv_cfs      = map(lambda k: pmt * k, factors)
    pv_pr       = fv * 1 / pow(1 + i, n)
    pv          = sum([pv_pr] + pv_cfs)
    return pv


def calc_fv(pv, i, n, pmt):
    exponents   = range(1, n + 1)
    factors     = map(lambda k: pow(1 + i, k), exponents)
    fv_cfs      = map(lambda k: pmt * k, factors)
    fv_pr       = pv * pow(1 + i, n)
    fv          = sum([fv_pr] + fv_cfs)

def calc_pmt(pv, fv, i, n):
    init_pmt    = 50.0
    epsilon     = 0.0001 * pv
    step_size   = 10.0
    cnt         = True
    pmt         = init_pmt
    last_up     = True
    while cnt:
        cpv     = calc_pv(fv, i, n, pmt)
        if abs(cpv - pv) <= epsilon:
            cnt = False
        elif cpv - pv < 0:
            if not last_up:
                step_size       /= 2
                last_up         = True
            pmt                 += step_size
        elif cpv - pv > 0:
            if last_up:
                step_size       /= 2
                last_up         = False
            pmt                 -= step_size
    return pmt

def calc_i(pv, fv, n, pmt):
    init_i      = 0.05
    epsilon     = 0.0001 * pv
    step_size   = 0.01
    cnt         = True
    i           = init_i
    last_up     = True
    while cnt:
        cpv     = calc_pv(fv, i, n, pmt)
        if abs(cpv - pv) <= epsilon:
            cnt = False
        elif cpv - pv > 0:
            if not last_up:
                step_size       /= 2
                last_up         = True
            i   += step_size
        elif cpv - pv < 0:
            if last_up:
                step_size       /= 2
                last_up         = False
            i   -= step_size
    return i

def calc_n(pv, fv, i, pmt):
    init_n      = 100.0
    epsilon     = 0.00001 * pv
    step_size   = 10.0
    cnt         = True
    n           = init_n
    last_up     = True
    while cnt:
        cpv     = calc_pv(fv, i, n, pmt)
        if abs(cpv - pv) <= epsilon:
            cnt = False
        elif cpv - pv > 0:
            if not last_up:
                step_size       = step_size / 2
                last_up         = True
            n   += step_size
        elif cpv - pv < 0:
            if last_up:
                step_size       = step_size / 2
                last_up         = False
            n   -= step_size
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

def  calc_zero():
    pass

def compound(p, i, n):
    c = p * pow(1 + i, n)
    return c
