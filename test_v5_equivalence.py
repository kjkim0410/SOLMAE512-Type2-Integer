"""Differential tests including ties, signed limits, and exceptions."""
import random
from benchmark_campaign_v5 import load

def outcome(fn):
    try:
        x=fn()
        if hasattr(x,'re'): return ((x.re.raw,x.re.frac,x.re.bits),(x.im.raw,x.im.frac,x.im.bits))
        if hasattr(x,'raw'): return x.raw,x.frac,x.bits
        return x
    except Exception as e: return type(e).__name__

def main():
    a,b=load('solmae-fxp-integer'),load('solmae-full-integer'); rng=random.Random(20261013); checks=0
    for x in list(range(-128,129))+[-2**63,2**63-1,0,-1,1]+[rng.randrange(-2**63,2**63) for _ in range(1000)]:
        for shift in [-3,0,1,2,16,63,64,100]:
            assert a.fx.rshift_round(x,shift)==b.fx.rshift_round(x,shift); checks+=1
    for _ in range(10000):
        raw1=rng.choice([rng.randrange(-2**63,2**63),0,1,-1]); raw2=rng.choice([rng.randrange(-2**63,2**63),0,1,-1])
        f1,f2,out=[rng.randrange(0,61) for _ in range(3)]; bits1=rng.choice([16,32,64]); bits2=rng.choice([16,32,64])
        raw1=max(-(1<<(bits1-1)),min((1<<(bits1-1))-1,raw1)); raw2=max(-(1<<(bits2-1)),min((1<<(bits2-1))-1,raw2))
        for op in ['add','sub','mul','div']:
            def calc(m): return getattr(m.fx,op)(m.fx.Fxp(raw1,f1,bits1),m.fx.Fxp(raw2,f2,bits2),out)
            assert outcome(lambda:calc(a))==outcome(lambda:calc(b)),(op,raw1,raw2,f1,f2,out); checks+=1
        assert outcome(lambda:a.fx.Fxp(raw1,f1,bits1).rescale(out))==outcome(lambda:b.fx.Fxp(raw1,f1,bits1).rescale(out)); checks+=1
        def complex_calc(m):
            x=m.cx.CFxp(m.fx.Fxp(raw1,f1,bits1),m.fx.Fxp(raw2,f2,bits2))
            y=m.cx.CFxp(m.fx.Fxp(raw2,f2,bits2),m.fx.Fxp(raw1,f1,bits1))
            return m.cx.cmul(x,y,out)
        assert outcome(lambda:complex_calc(a))==outcome(lambda:complex_calc(b)); checks+=1
    for _ in range(1000):
        vals=[rng.randrange(-100000,100001) for _ in range(4)]
        frac=rng.randrange(0,31); out=rng.randrange(0,31)
        for op in ['cadd','csub','cmul','cdiv']:
            def calc(m):
                x=m.cx.CFxp(m.fx.Fxp(vals[0],frac),m.fx.Fxp(vals[1],frac))
                y=m.cx.CFxp(m.fx.Fxp(vals[2],frac),m.fx.Fxp(vals[3],frac))
                return getattr(m.cx,op)(x,y,out)
            assert outcome(lambda:calc(a))==outcome(lambda:calc(b)),(op,vals,frac,out)
            checks+=1
    print(f'PASS: {checks} differential comparisons')
if __name__=='__main__': main()
