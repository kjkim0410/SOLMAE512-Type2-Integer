from dataclasses import dataclass
from fxp import Fxp, add, sub, mul, div
@dataclass(frozen=True)
class CFxp:
    re: Fxp
    im: Fxp
    def conj(self): return CFxp(self.re, Fxp(-self.im.raw,self.im.frac,self.im.bits))
    def to_complex(self): return complex(self.re.to_float(),self.im.to_float())

def from_complex(z, frac): return CFxp(Fxp.from_float(z.real,frac),Fxp.from_float(z.imag,frac))
def cadd(a,b,fout): return CFxp(add(a.re,b.re,fout),add(a.im,b.im,fout))
def csub(a,b,fout): return CFxp(sub(a.re,b.re,fout),sub(a.im,b.im,fout))
def cmul(a,b,fout):
    ac=mul(a.re,b.re,fout); bd=mul(a.im,b.im,fout)
    ad=mul(a.re,b.im,fout); bc=mul(a.im,b.re,fout)
    return CFxp(sub(ac,bd,fout),add(ad,bc,fout))
def cdiv(a,b,fout,tmp_frac=50):
    # (a * conj(b)) / |b|^2
    bc=b.conj(); num=cmul(a,bc,tmp_frac)
    den=add(mul(b.re,b.re,tmp_frac),mul(b.im,b.im,tmp_frac),tmp_frac)
    return CFxp(div(num.re,den,fout),div(num.im,den,fout))
