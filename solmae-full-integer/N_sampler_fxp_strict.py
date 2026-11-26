"""Strict integer-only N-sampler for SOLMAE-512 Type-2 FxP."""

# Shared reference modules; keep the selected implementation first.
import sys as _solmae_sys
from pathlib import Path as _solmae_Path
_solmae_shared = str(_solmae_Path(__file__).resolve().parent.parent / "solmae-ref")
if _solmae_shared not in _solmae_sys.path:
    _solmae_sys.path.append(_solmae_shared)

from os import urandom
from params import SOLMAE_D
from fxp import Fxp, sub, mul
from fxp_math_strict import sqrt_fxp, log_fxp, sin_fxp, cos_fxp
from fxp_formats import F

# Pre-quantized constants. No runtime float construction.
LN2 = Fxp(195103586505167, 48)        # round(ln(2)*2^48)
# TWO_PI = Fxp(452751216685240, 46)     # round(2*pi*2^46) 1st and old
TWO_PI = Fxp(442139859501778, 46)     # round(2 * pi * 2**46) : revised and correct
D = Fxp.from_int(SOLMAE_D, F['RADIUS'])

def ffsll(word):
    num = 0
    if word & 0xffffffff == 0: num += 32; word >>= 32
    if word & 0xffff == 0: num += 16; word >>= 16
    if word & 0xff == 0: num += 8; word >>= 8
    if word & 0xf == 0: num += 4; word >>= 4
    if word & 0x3 == 0: num += 2; word >>= 2
    if word & 0x1 == 0: num += 1
    return num



def n_sampler_fxp_strict(randombytes=urandom):
    coeffs = []
    for _ in range(SOLMAE_D // 2):
        u = int.from_bytes(randombytes(8), 'little')
        u53 = u & 0x1fffffffffffff
        uf = Fxp(u53, 53).rescale(F['ANGLE'])
        uf = mul(TWO_PI, uf, F['ANGLE'])

        v = int.from_bytes(randombytes(8), 'little')
        vf_raw = (1 << 53) + (v & 0x1fffffffffffff)
        vf = Fxp(vf_raw, 54).rescale(F['RADIUS'])

        e1 = int.from_bytes(randombytes(8), 'little')
        e2 = int.from_bytes(randombytes(8), 'little')
        x = 64 + ffsll(e2); y = ffsll(e1); c = ((~e1) & (e1 - 1)) >> 63
        geom = (x & (-(c & 1))) ^ (y & (~(-(c & 1))))

        gl = mul(LN2, Fxp.from_int(geom, F['RADIUS']), F['RADIUS'])
        lv = log_fxp(vf, F['RADIUS'])
        inside = mul(D, sub(gl, lv, F['RADIUS']), F['RADIUS'])
        radius = sqrt_fxp(inside, F['RADIUS'])
        co = cos_fxp(uf, F['RADIUS'])
        si = sin_fxp(uf, F['RADIUS'])
        coeffs.append(mul(radius, co, F['COEFF']))
        coeffs.append(mul(radius, si, F['COEFF']))
    return coeffs
