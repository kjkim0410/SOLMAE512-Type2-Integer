"""Integer-only discrete Gaussian sampler for strict Type-2 SOLMAE signing.

This is the repository samplerz logic rewritten so the center, sigma-derived
constants, exp polynomial, and Bernoulli rejection test use only integers.
"""
from os import urandom
from fxp import Fxp, rshift_round

RCDT_PREC = 72
RCDT = [
3024686241123004913666,1564742784480091954050,636254429462080897535,
199560484645026482916,47667343854657281903,8595902006365044063,
1163297957344668388,117656387352093658,8867391802663976,
496969357462633,20680885154299,638331848991,14602316184,247426747,
3104126,28824,198,1]
C = [
0x00000004741183A3,0x00000036548CFC06,0x0000024FDCBF140A,
0x0000171D939DE045,0x0000D00CF58F6F84,0x000680681CF796E3,
0x002D82D8305B0FEA,0x011111110E066FD0,0x0555555555070F00,
0x155555555581FF00,0x400000000002B400,0x7FFFFFFFFFFF4800,
0x8000000000000000]
Q = 63
ONE = 1 << Q
LN2 = 6393154322601327616  # ln(2)*2^63
# SOLMAE-512 fixed constants, generated offline from decimal strings.
SIGMA_Q52 = 5944751508129055
SIGMIN_Q52 = 5754851361258101
MAXSIG_Q52 = 8198803121627988


def _div_round(num, den):
    if den == 0: raise ZeroDivisionError
    sign = -1 if (num < 0) ^ (den < 0) else 1
    a,b=abs(num),abs(den); q,r=divmod(a,b)
    if 2*r > b or (2*r == b and (q&1)): q += 1
    return -q if sign < 0 else q


def _ratio_q63(a_raw, a_frac, b_raw, b_frac):
    shift = Q + b_frac - a_frac
    num = a_raw << max(shift,0)
    den = b_raw << max(-shift,0)
    return _div_round(num, den)


def basesampler(randombytes=urandom):
    u = int.from_bytes(randombytes(RCDT_PREC >> 3), 'little')
    z0 = 0
    for elt in RCDT: z0 += int(u < elt)
    return z0


def approxexp_q63(x_raw, ccs_raw):
    y = C[0]
    z = x_raw
    for elt in C[1:]:
        y = elt - ((z * y) >> 63)
    z = ccs_raw << 1
    return (z * y) >> 63


def berexp_q63(x_raw, ccs_raw, randombytes=urandom):
    if x_raw < 0: x_raw = 0
    s = x_raw // LN2
    r = x_raw - s * LN2
    if s > 63: s = 63
    z = (approxexp_q63(r, ccs_raw) - 1) >> s
    w = 0
    for i in range(56, -8, -8):
        p = int.from_bytes(randombytes(1), 'little')
        w = p - ((z >> i) & 0xFF)
        if w: break
    return w < 0


def samplerz_fxp(mu, randombytes=urandom):
    """Sample with SOLMAE-512 sigma=1.32 and sigmin=1.277833..."""
    s = mu.raw >> mu.frac
    rem = mu.raw - (s << mu.frac)
    r_q63 = rshift_round(rem, mu.frac - Q)

    # dss=1/(2*sigma^2), ccs=sigmin/sigma, invmax=1/(2*MAX_SIGMA^2)
    sigma2 = SIGMA_Q52 * SIGMA_Q52
    dss_q63 = _div_round(1 << (Q + 104), 2 * sigma2)
    ccs_q63 = _ratio_q63(SIGMIN_Q52, 52, SIGMA_Q52, 52)
    max2 = MAXSIG_Q52 * MAXSIG_Q52
    invmax_q63 = _div_round(1 << (Q + 104), 2 * max2)

    while True:
        z0 = basesampler(randombytes)
        b = int.from_bytes(randombytes(1), 'little') & 1
        z = b + (2*b - 1) * z0
        a = (z << Q) - r_q63
        x1 = rshift_round(a * a * dss_q63, 126)
        x2 = z0 * z0 * invmax_q63
        x_raw = x1 - x2
        if berexp_q63(x_raw, ccs_q63, randombytes):
            return z + s
