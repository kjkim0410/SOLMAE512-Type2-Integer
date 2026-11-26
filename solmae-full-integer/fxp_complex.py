from dataclasses import dataclass
from fxp import Fxp, add, sub, mul, div, _check, rshift_round
@dataclass(frozen=True, slots=True)
class CFxp:
    re: Fxp
    im: Fxp
    def conj(self): return CFxp(self.re, Fxp(-self.im.raw,self.im.frac,self.im.bits))
    def to_complex(self): return complex(self.re.to_float(),self.im.to_float())

def from_complex(z, frac): return CFxp(Fxp.from_float(z.real,frac),Fxp.from_float(z.imag,frac))
def cadd(a,b,fout): return CFxp(add(a.re,b.re,fout),add(a.im,b.im,fout))
def csub(a,b,fout): return CFxp(sub(a.re,b.re,fout),sub(a.im,b.im,fout))
def _product_raw(a, b, fout):
    # Preserve v4's separate product rounding and intermediate overflow checks.
    return _check(rshift_round(a.raw*b.raw, a.frac+b.frac-fout), a.bits)

def cmul(a,b,fout):
    ac=_product_raw(a.re,b.re,fout); bd=_product_raw(a.im,b.im,fout)
    ad=_product_raw(a.re,b.im,fout); bc=_product_raw(a.im,b.re,fout)
    # v4 sub(ac,bd) uses ac.bits; add(ad,bc) uses ad.bits.
    return CFxp(Fxp(ac-bd,fout,a.re.bits),Fxp(ad+bc,fout,a.re.bits))

def cdiv(a,b,fout,tmp_frac=50):
    # (a * conj(b)) / |b|^2
    bc=b.conj(); num=cmul(a,bc,tmp_frac)
    den=add(mul(b.re,b.re,tmp_frac),mul(b.im,b.im,tmp_frac),tmp_frac)
    return CFxp(div(num.re,den,fout),div(num.im,den,fout))
