"""Type-2 quantized N-sampler for SOLMAE-512.
The arithmetic state is fixed-point; log/sin/cos currently use quantized reference
functions from fxp_math for controlled Float-vs-FxP experiments.
"""

# Shared reference modules; keep the selected implementation first.
import sys as _solmae_sys
from pathlib import Path as _solmae_Path
_solmae_shared = str(_solmae_Path(__file__).resolve().parent.parent / "solmae-ref")
if _solmae_shared not in _solmae_sys.path:
    _solmae_sys.path.append(_solmae_shared)

from os import urandom
from params import SOLMAE_D
from N_sampler import ffsll
from fxp import Fxp, add, sub, mul
from fxp_math import sqrt_fxp, log_fxp, sin_fxp, cos_fxp
from fxp_formats import F

LN2=Fxp.from_float(0.693147180559945309417,F['RADIUS'])
TWO_PI=Fxp.from_float(6.28318530717958647692,F['ANGLE'])
D=Fxp.from_int(SOLMAE_D,F['RADIUS'])

def n_sampler_fxp(randombytes=urandom):
    coeffs=[]
    for _ in range(SOLMAE_D//2):
        u=int.from_bytes(randombytes(8),'little')
        u53=u & 0x1fffffffffffff
        uf=Fxp(u53,F['UNIFORM']).rescale(F['ANGLE'])
        uf=mul(TWO_PI,uf,F['ANGLE'])
        v=int.from_bytes(randombytes(8),'little')
        vf_raw=(1<<53)+(v & 0x1fffffffffffff)  # 1/2 + mantissa*2^-54, represented Q54
        vf=Fxp(vf_raw,54).rescale(F['RADIUS'])
        e1=int.from_bytes(randombytes(8),'little'); e2=int.from_bytes(randombytes(8),'little')
        x=64+ffsll(e2); y=ffsll(e1); c=((~e1)&(e1-1))>>63
        geom=(x & (-(c&1))) ^ (y & (~(-(c&1))))
        gl=mul(LN2,Fxp.from_int(geom,F['RADIUS']),F['RADIUS'])
        lv=log_fxp(vf,F['RADIUS'])
        inside=mul(D,sub(gl,lv,F['RADIUS']),F['RADIUS'])
        radius=sqrt_fxp(inside,F['RADIUS'])
        co=cos_fxp(uf,F['RADIUS']); si=sin_fxp(uf,F['RADIUS'])
        coeffs.append(mul(radius,co,F['COEFF']).to_float())
        coeffs.append(mul(radius,si,F['COEFF']).to_float())
    return coeffs
