"""Strict integer-only Type-2 fixed-point FFT for SOLMAE-512."""

# Shared reference modules; keep the selected implementation first.
import sys as _solmae_sys
from pathlib import Path as _solmae_Path
_solmae_shared = str(_solmae_Path(__file__).resolve().parent.parent / "solmae-ref")
if _solmae_shared not in _solmae_sys.path:
    _solmae_sys.path.append(_solmae_shared)

from common import split, merge
from fxp import Fxp, _check
from fxp_complex import CFxp, cadd, csub, cmul, cdiv
from fxp_formats import F
from fft_roots_fxp_constants import ROOTS_RAW, TWIDDLE_FRAC

ROOTS = {
    n: [CFxp(Fxp(re, TWIDDLE_FRAC), Fxp(im, TWIDDLE_FRAC)) for re, im in arr]
    for n, arr in ROOTS_RAW.items()
}


def fft_fxp(coeffs, frac=F['FFT_KEY']):
    n = len(coeffs)
    if n > 2:
        a, b = split(coeffs)
        return merge_fft_fxp([fft_fxp(a, frac), fft_fxp(b, frac)], frac)
    if n == 2:
        x0 = Fxp.from_int(coeffs[0], frac)
        x1 = Fxp.from_int(coeffs[1], frac)
        return [CFxp(x0, x1), CFxp(x0, Fxp(-x1.raw, x1.frac, x1.bits))]
    raise ValueError('length must be a power of two >=2')


def split_fft_fxp(v, frac):
    n = len(v); w = ROOTS[n]; a = []; b = []
    for i in range(n // 2):
        s = cadd(v[2*i], v[2*i+1], frac)
        d = csub(v[2*i], v[2*i+1], frac)
        hs = CFxp(Fxp(s.re.raw >> 1, frac), Fxp(s.im.raw >> 1, frac))
        hd = CFxp(Fxp(d.re.raw >> 1, frac), Fxp(d.im.raw >> 1, frac))
        a.append(hs)
        b.append(cmul(hd, w[2*i].conj(), frac))
    return a, b


def merge_fft_fxp(pair, frac):
    a, b = pair; n = 2 * len(a); w = ROOTS[n]; out = [None] * n
    for i in range(n // 2):
        wb = cmul(w[2*i], b[i], frac)
        out[2*i] = cadd(a[i], wb, frac)
        out[2*i+1] = csub(a[i], wb, frac)
    return out


def ifft_fxp(v, frac=F['COEFF']):
    n = len(v)
    if n > 2:
        a, b = split_fft_fxp(v, v[0].re.frac)
        return merge([ifft_fxp(a, frac), ifft_fxp(b, frac)])
    if n == 2:
        return [v[0].re.rescale(frac), v[0].im.rescale(frac)]
    raise ValueError


def solmae_ifft_fxp(v, frac=F['COEFF']):
    """Inverse of SOLMAE half-FFT real storage, all fixed point."""
    n = len(v); h = n // 2
    full = [CFxp(v[i], v[i+h]) for i in range(h)]
    full += [CFxp(v[i], Fxp(-v[i+h].raw, v[i+h].frac, v[i+h].bits)) for i in range(h)]
    return ifft_fxp(full, frac)


def adj_fft_fxp(v): return [z.conj() for z in v]
def add_fft_fxp(a,b,fout): return [cadd(x,y,fout) for x,y in zip(a,b)]
def sub_fft_fxp(a,b,fout): return [csub(x,y,fout) for x,y in zip(a,b)]
def mul_fft_fxp(a,b,fout): return [cmul(x,y,fout) for x,y in zip(a,b)]
def _round_div(num, den):
    if den == 0: raise ZeroDivisionError
    sign=-1 if (num<0) ^ (den<0) else 1
    aa,bb=abs(num),abs(den); q,r=divmod(aa,bb)
    if 2*r > bb or (2*r == bb and (q&1)): q += 1
    return -q if sign<0 else q

def cdiv_wide(a,b,fout):
    # Wide integer intermediate, 64-bit check only on final Type-2 words.
    fa=a.re.frac; fb=b.re.frac
    nr=a.re.raw*b.re.raw + a.im.raw*b.im.raw
    ni=a.im.raw*b.re.raw - a.re.raw*b.im.raw
    den=b.re.raw*b.re.raw + b.im.raw*b.im.raw
    sh=fout + fb - fa
    if sh>=0:
        rr=_round_div(nr<<sh,den); ri=_round_div(ni<<sh,den)
    else:
        rr=_round_div(nr,den<<(-sh)); ri=_round_div(ni,den<<(-sh))
    return CFxp(Fxp(_check(rr),fout),Fxp(_check(ri),fout))

def div_fft_fxp(a,b,fout): return [cdiv_wide(x,y,fout) for x,y in zip(a,b)]
