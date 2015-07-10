
def calc_pv(fv, i, n, pmt):
    exponents   = range(1, n + 1)
    factors     = map(lambda k: 1 / pow(1 + i, k), exponents)
    pv_cfs      = map(lambda k: pmt * k, factors)
    pv_pr       = fv * 1 / pow(1 + i, n)
    pv          = sum([pv_pr] ++ pv_cfs)
    return pv


def calc_fv(pv, i, n, pmt):
    exponents   = range(1, n + 1)
    factors     = map(lambda k: pow(1 + i, k), exponents)
    fv_cfs      = map(lambda k: pmt * k, factors)
    fv_pr       = pv * pow(1 + i, n)
    fv          = sum([fv_pr] ++ fv_cfs)

def calc_pmt(pv, fv, i, n):
    pass

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
                step_size /= 2
            i += step_size

        elif cpv - pv < 0:
            if last_up:
                step_size /= 2
            i -= step_size
    return i

def calc_n(pm, fv, i, pmt):
    pass

def  calc_zero():
    pass

def compound(p, i, n):
    c = p * pow(1 + i, n)
    return c
