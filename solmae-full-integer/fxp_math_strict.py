"""Strict integer-only elementary functions for SOLMAE Type-2 fixed point.

No Python float, ``math.sin/cos/log`` or NumPy floating-point operation is used.
All constants are pre-quantized integers.  Functions accept/return ``Fxp``.
"""
from math import isqrt
from fxp import Fxp, rshift_round

Q = 60
LN2_Q = 799144290325165979
PI_Q = 3622009729038561421
TWO_PI_Q = 7244019458077122842
HALF_PI_Q = 1811004864519280711
QUARTER_PI_Q = 905502432259640355
ONE_Q = 1 << Q


def _mulq(a, b, q=Q):
    return rshift_round(a * b, q)


def _divq(a, b, q=Q):
    if b == 0:
        raise ZeroDivisionError
    sign = -1 if (a < 0) ^ (b < 0) else 1
    aa, bb = abs(a), abs(b)
    num = aa << q
    qq, rr = divmod(num, bb)
    if 2 * rr > bb or (2 * rr == bb and (qq & 1)):
        qq += 1
    return -qq if sign < 0 else qq


def _to_q60(x):
    return rshift_round(x.raw, x.frac - Q)


def _from_q60(raw, fout, bits=64):
    return Fxp(rshift_round(raw, Q - fout), fout, bits)


def sqrt_fxp(x, fout=None):
    if x.raw < 0:
        raise ValueError("sqrt of negative")
    if fout is None:
        fout = x.frac
    shift = 2 * fout - x.frac
    n = x.raw << shift if shift >= 0 else x.raw >> (-shift)
    return Fxp(isqrt(n), fout, x.bits)


def log_fxp(x, fout):
    """Natural log using normalization and 15-term atanh series.

    x = 2**e * m, m in [1,2), z=(m-1)/(m+1),
    ln(m)=2*(z + z^3/3 + ... + z^29/29).
    """
    if x.raw <= 0:
        raise ValueError("log domain error")

    xr = x.raw
    # Find e such that m in [1,2), without floating point.
    k = xr.bit_length() - 1
    e = k - x.frac
    if e >= 0:
        m_q = rshift_round(xr, e + x.frac - Q)
    else:
        m_q = rshift_round(xr, x.frac + e - Q)
    # Robust normalization in case rounding put m on a boundary.
    while m_q < ONE_Q:
        m_q <<= 1
        e -= 1
    while m_q >= (2 << Q):
        m_q = rshift_round(m_q, 1)
        e += 1

    z = _divq(m_q - ONE_Q, m_q + ONE_Q)
    z2 = _mulq(z, z)
    term = z
    acc = term
    for d in range(3, 30, 2):
        term = _mulq(term, z2)
        acc += term // d
    ln_m = acc << 1
    out = ln_m + e * LN2_Q
    return _from_q60(out, fout, x.bits)


def _reduce_angle_q60(a):
    # Map to [-pi, pi].
    a %= TWO_PI_Q
    if a > PI_Q:
        a -= TWO_PI_Q
    return a


def _sin_poly_q60(x):
    # Taylor through x^13 on |x| <= pi/4.
    x2 = _mulq(x, x)
    # Horner for 1 - x^2/3! + x^4/5! - ... + x^12/13!
    c13 = ONE_Q // 6227020800
    c11 = -(ONE_Q // 39916800)
    c9 = ONE_Q // 362880
    c7 = -(ONE_Q // 5040)
    c5 = ONE_Q // 120
    c3 = -(ONE_Q // 6)
    p = c13
    for c in (c11, c9, c7, c5, c3, ONE_Q):
        p = c + _mulq(p, x2)
    return _mulq(x, p)


def _cos_poly_q60(x):
    # Taylor through x^14 on |x| <= pi/4.
    x2 = _mulq(x, x)
    c14 = -(ONE_Q // 87178291200)
    c12 = ONE_Q // 479001600
    c10 = -(ONE_Q // 3628800)
    c8 = ONE_Q // 40320
    c6 = -(ONE_Q // 720)
    c4 = ONE_Q // 24
    c2 = -(ONE_Q // 2)
    p = c14
    for c in (c12, c10, c8, c6, c4, c2, ONE_Q):
        p = c + _mulq(p, x2)
    return p


def sin_fxp(x, fout):
    a = _reduce_angle_q60(_to_q60(x))
    # Reduce to [-pi/4, pi/4] with exact quadrant identities.
    if a > HALF_PI_Q:
        a = PI_Q - a
    elif a < -HALF_PI_Q:
        a = -PI_Q - a
    if a > QUARTER_PI_Q:
        y = HALF_PI_Q - a
        raw = _cos_poly_q60(y)
    elif a < -QUARTER_PI_Q:
        y = HALF_PI_Q + a
        raw = -_cos_poly_q60(y)
    else:
        raw = _sin_poly_q60(a)
    return _from_q60(raw, fout, x.bits)


def cos_fxp(x, fout):
    a = _reduce_angle_q60(_to_q60(x))
    sign = 1
    if a > HALF_PI_Q:
        a = PI_Q - a
        sign = -1
    elif a < -HALF_PI_Q:
        a = -PI_Q - a
        sign = -1
    if a > QUARTER_PI_Q:
        y = HALF_PI_Q - a
        raw = _sin_poly_q60(y)
    elif a < -QUARTER_PI_Q:
        y = HALF_PI_Q + a
        raw = _sin_poly_q60(y)
    else:
        raw = _cos_poly_q60(a)
    raw *= sign
    return _from_q60(raw, fout, x.bits)
