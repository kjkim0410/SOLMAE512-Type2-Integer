"""Type-2 fixed-point FFT compatible with the repository's recursive FFT layout."""

# Shared reference modules; keep the selected implementation first.
import sys as _solmae_sys
from pathlib import Path as _solmae_Path
_solmae_shared = str(_solmae_Path(__file__).resolve().parent.parent / "solmae-ref")
if _solmae_shared not in _solmae_sys.path:
    _solmae_sys.path.append(_solmae_shared)

from common import split, merge
from fft_constants import roots_dict
from fxp import Fxp
from fxp_complex import CFxp, from_complex, cadd, csub, cmul, cdiv
from fxp_formats import F

# Constants are quantized once when the module loads. For a deployment implementation,
# export these integers to a generated constants file so no float conversion occurs at runtime.
ROOTS={n:[from_complex(z,F['TWIDDLE']) for z in arr] for n,arr in roots_dict.items()}

def fft_fxp(coeffs, frac=F['FFT_KEY']):
    n=len(coeffs)
    if n>2:
        a,b=split(coeffs)
        return merge_fft_fxp([fft_fxp(a,frac),fft_fxp(b,frac)],frac)
    if n==2:
        x0=Fxp.from_int(coeffs[0],frac); x1=Fxp.from_int(coeffs[1],frac)
        z0=CFxp(x0,x1); z1=CFxp(x0,Fxp(-x1.raw,x1.frac,x1.bits))
        return [z0,z1]
    raise ValueError('length must be a power of two >=2')

def split_fft_fxp(v,frac):
    n=len(v); w=ROOTS[n]; a=[]; b=[]
    for i in range(n//2):
        s=cadd(v[2*i],v[2*i+1],frac)
        d=csub(v[2*i],v[2*i+1],frac)
        # exact divide by 2 via rescale/rshift of raw at same frac
        hs=CFxp(Fxp(s.re.raw>>1,frac),Fxp(s.im.raw>>1,frac))
        hd=CFxp(Fxp(d.re.raw>>1,frac),Fxp(d.im.raw>>1,frac))
        a.append(hs); b.append(cmul(hd,w[2*i].conj(),frac))
    return a,b

def merge_fft_fxp(pair,frac):
    a,b=pair; n=2*len(a); w=ROOTS[n]; out=[None]*n
    for i in range(n//2):
        wb=cmul(w[2*i],b[i],frac)
        out[2*i]=cadd(a[i],wb,frac); out[2*i+1]=csub(a[i],wb,frac)
    return out

def ifft_fxp(v,frac=F['COEFF']):
    n=len(v)
    if n>2:
        a,b=split_fft_fxp(v,v[0].re.frac)
        return merge([ifft_fxp(a,frac),ifft_fxp(b,frac)])
    if n==2:
        return [v[0].re.rescale(frac),v[0].im.rescale(frac)]
    raise ValueError

def adj_fft_fxp(v): return [z.conj() for z in v]
def add_fft_fxp(a,b,fout): return [cadd(x,y,fout) for x,y in zip(a,b)]
def sub_fft_fxp(a,b,fout): return [csub(x,y,fout) for x,y in zip(a,b)]
def mul_fft_fxp(a,b,fout): return [cmul(x,y,fout) for x,y in zip(a,b)]
def div_fft_fxp(a,b,fout): return [cdiv(x,y,fout) for x,y in zip(a,b)]
def to_complex_list(v): return [z.to_complex() for z in v]
def from_complex_list(v,frac): return [from_complex(z,frac) for z in v]
