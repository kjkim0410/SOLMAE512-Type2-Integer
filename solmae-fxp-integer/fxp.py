"""Type-2 fixed-point arithmetic for the SOLMAE-512 research prototype.

A value x is represented by signed integer X and compile-time fraction length f:
    x = X * 2**(-f)
All arithmetic uses Python integers; overflow is checked against a conceptual
signed 64-bit word. This mirrors the Type-2 idea in the fixed-point Falcon work:
uniform word length, per-operand fraction length.
"""
from dataclasses import dataclass
from math import floor

WORD_BITS = 64

class FxpOverflow(OverflowError):
    pass

def _check(raw, bits=WORD_BITS):
    lo = -(1 << (bits-1)); hi = (1 << (bits-1)) - 1
    if raw < lo or raw > hi:
        raise FxpOverflow(f"{raw} does not fit signed {bits}-bit")
    return raw

def rshift_round(x, s):
    if s <= 0:
        return x << (-s)
    a = abs(x)
    q, r = divmod(a, 1 << s)
    half = 1 << (s-1)
    if r > half or (r == half and (q & 1)):
        q += 1
    return -q if x < 0 else q

@dataclass(frozen=True)
class Fxp:
    raw: int
    frac: int
    bits: int = WORD_BITS
    def __post_init__(self):
        _check(self.raw, self.bits)
    @staticmethod
    def from_float(x, frac, bits=WORD_BITS):
        return Fxp(_check(int(round(float(x) * (1 << frac))), bits), frac, bits)
    @staticmethod
    def from_int(x, frac=0, bits=WORD_BITS):
        return Fxp(_check(int(x) << frac, bits), frac, bits)
    def to_float(self):
        return self.raw / float(1 << self.frac)
    def rescale(self, frac):
        raw = rshift_round(self.raw, self.frac-frac)
        return Fxp(_check(raw, self.bits), frac, self.bits)

def add(a, b, fout=None):
    if fout is None: fout=max(a.frac,b.frac)
    aa=a.rescale(fout); bb=b.rescale(fout)
    return Fxp(_check(aa.raw+bb.raw,a.bits),fout,a.bits)

def sub(a,b,fout=None):
    if fout is None: fout=max(a.frac,b.frac)
    aa=a.rescale(fout); bb=b.rescale(fout)
    return Fxp(_check(aa.raw-bb.raw,a.bits),fout,a.bits)

def mul(a,b,fout):
    raw=rshift_round(a.raw*b.raw,a.frac+b.frac-fout)
    return Fxp(_check(raw,a.bits),fout,a.bits)

def div(a,b,fout):
    if b.raw==0: raise ZeroDivisionError
    shift=fout+b.frac-a.frac
    num=a.raw << max(shift,0)
    den=b.raw << max(-shift,0)
    sign=-1 if (num<0) ^ (den<0) else 1
    q,r=divmod(abs(num),abs(den))
    if 2*r > abs(den) or (2*r==abs(den) and (q&1)): q+=1
    raw=sign*q
    return Fxp(_check(raw,a.bits),fout,a.bits)

def floor_int(a):
    return a.raw >> a.frac
