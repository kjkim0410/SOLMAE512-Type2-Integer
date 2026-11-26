"""Fixed-point elementary functions used by the SOLMAE prototype.
Current implementation uses integer Newton for sqrt. log/sin/cos are quantized
reference functions for validation; replacing them with polynomial/table versions
is the next step toward a strictly FPA-free implementation.
"""
import math
from fxp import Fxp

def sqrt_fxp(x,fout=None):
    if x.raw<0: raise ValueError('sqrt of negative')
    if fout is None: fout=x.frac
    # integer sqrt of scaled value: y_raw ~= sqrt(x)*2^fout
    shift=2*fout-x.frac
    n=x.raw << shift if shift>=0 else x.raw >> (-shift)
    return Fxp(math.isqrt(n),fout,x.bits)

def log_fxp(x,fout):
    if x.raw<=0: raise ValueError
    return Fxp.from_float(math.log(x.to_float()),fout,x.bits)
def sin_fxp(x,fout): return Fxp.from_float(math.sin(x.to_float()),fout,x.bits)
def cos_fxp(x,fout): return Fxp.from_float(math.cos(x.to_float()),fout,x.bits)
