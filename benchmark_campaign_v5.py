"""Paired v4/v5 research benchmark; Python 3.10+. Not a production signer."""
import argparse, csv, hashlib, importlib, json, platform, statistics, sys, time
from pathlib import Path
from types import SimpleNamespace
ROOT = Path(__file__).resolve().parent

class Stream:
    """Deterministic SHA-256 counter stream FOR TESTS ONLY; buffered across calls."""
    def __init__(self, seed): self.seed=seed; self.counter=0; self.buf=b''; self.used=0
    def __call__(self,n):
        while len(self.buf)<n:
            self.buf+=hashlib.sha256(self.seed+self.counter.to_bytes(8,'big')).digest()
            self.counter+=1
        out,self.buf=self.buf[:n],self.buf[n:]; self.used+=n
        return out

def load(version):
    # Retained modules keep their own globals; clear only prototype modules.
    names={p.stem for directory in (ROOT/version, ROOT/'solmae-ref', ROOT/'solmae-full-integer') for p in directory.glob('*.py')}
    for name in names: sys.modules.pop(name,None)
    sys.path.insert(0,str(ROOT/'solmae-ref'))
    sys.path.insert(0,str(ROOT/'solmae-full-integer'))
    sys.path.insert(0,str(ROOT/version))
    try:
        fx=importlib.import_module('fxp'); cx=importlib.import_module('fxp_complex')
        sk=importlib.import_module('s_keygen_fxp_strict')
        sig=importlib.import_module('solmae_fxp_strict')
        kg=importlib.import_module('s_keygen')
        return SimpleNamespace(fx=fx,cx=cx,derive=sk.derive_fxp_key_from_basis,
                               sign=sig.sign_fxp_strict,verify=sig.verify,keygen=kg.keygen)
    finally:
        sys.path.pop(0)
        sys.path.pop(0)
        sys.path.pop(0)

def write_csv(path,rows):
    if not rows: return
    with path.open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--trials',nargs='+',type=int,default=[100,1000])
    ap.add_argument('--warmup',type=int,default=3)
    ap.add_argument('--message-bytes',type=int,default=512)
    ap.add_argument('--seed',default='SOLMAE-v5-campaign-2026')
    ap.add_argument('--out',type=Path,default=ROOT/'results')
    args=ap.parse_args()
    if min(args.trials)<1 or args.warmup<0 or args.message_bytes<0: ap.error('invalid count')
    args.out.mkdir(parents=True,exist_ok=True)
    v4,v5=load('solmae-fxp-integer'),load('solmae-full-integer'); impl={'v4':v4,'v5':v5}
    print('Generating shared SOLMAE-512 key...',flush=True)
    t=time.perf_counter(); ref,pk=v4.keygen(); keygen_ms=(time.perf_counter()-t)*1000
    basis=[ref.f,ref.g,ref.F,ref.G]; keys={}; derivation={}
    for name,m in impl.items():
        t=time.perf_counter(); keys[name]=m.derive(*basis); derivation[name]=(time.perf_counter()-t)*1000
    # Explicit numerical-key comparison, independent of signature success.
    def flatten(x):
        if hasattr(x,'raw'): return (x.raw,x.frac,x.bits)
        if hasattr(x,'re'): return (flatten(x.re),flatten(x.im))
        if isinstance(x,list): return [flatten(y) for y in x]
        return x
    if {k:flatten(v) for k,v in vars(keys['v4']).items()} != {k:flatten(v) for k,v in vars(keys['v5']).items()}:
        raise RuntimeError('Derived numerical keys differ')
    meta={'python':sys.version,'platform':platform.platform(),'processor':platform.processor(),
          'trials':args.trials,'warmup':args.warmup,'message_bytes':args.message_bytes,'seed':args.seed,
          'keygen_ms':keygen_ms,'derive_ms':derivation,'numerical_keys_equal':True,
          'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for v in ('solmae-fxp-integer','solmae-full-integer') for p in sorted((ROOT/v).glob('*.py'))},
          'basis_sha256':hashlib.sha256(json.dumps(basis).encode()).hexdigest(),
          'note':'Single shared key; deterministic test randomness; sign timing excludes verification. Keygen uses OS randomness.'}
    (args.out/'environment.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
    for i in range(args.warmup):
        for name,m in impl.items():
            sig=m.sign(keys[name],b'warmup',Stream(hashlib.sha256(f'warmup:{args.seed}:{i}'.encode()).digest()))
            if not m.verify(pk,b'warmup',sig): raise RuntimeError('warmup verification failed')
    summaries=[]; overall_ok=True
    for count in args.trials:
        rows=[]; times={name:[] for name in impl}; verified={name:0 for name in impl}; equal=0
        for i in range(count):
            domain=f'{args.seed}:{count}:{i}'.encode()
            msg=hashlib.shake_256(b'message:'+domain).digest(args.message_bytes)
            seed=hashlib.sha256(b'random:'+domain).digest(); sigs={}; usages={}
            for name in (('v4','v5') if i%2==0 else ('v5','v4')):
                m=impl[name]; rng=Stream(seed); err=''; sig=None; ok=False
                t=time.perf_counter_ns()
                try: sig=m.sign(keys[name],msg,rng)
                except Exception as e: err=f'{type(e).__name__}: {e}'
                elapsed=(time.perf_counter_ns()-t)/1e6
                if sig is not None:
                    try: ok=len(sig)==666 and m.verify(pk,msg,sig)
                    except Exception as e: err=f'verify {type(e).__name__}: {e}'
                if sig is not None: times[name].append(elapsed)
                verified[name]+=int(ok); sigs[name]=sig; usages[name]=rng.used
                rows.append({'implementation':name,'trials':count,'trial':i+1,'sign_ms':elapsed,
                             'verified':ok,'signature_bytes':len(sig) if sig else 0,
                             'random_bytes':rng.used,'error':err})
            match=sigs['v4'] is not None and sigs['v4']==sigs['v5'] and usages['v4']==usages['v5']
            equal+=int(match)
            for row in rows[-2:]: row['paired_signature_and_rng_equal']=match
            # Incremental write preserves evidence if a long campaign is interrupted.
            if (i+1)%10==0 or i+1==count:
                write_csv(args.out/f'campaign_{count}_trials.csv',rows)
                print(f'{count} trials: {i+1}/{count}; verified v4={verified["v4"]}, v5={verified["v5"]}; equal={equal}',flush=True)
        mean4=statistics.mean(times['v4']) if times['v4'] else 0
        mean5=statistics.mean(times['v5']) if times['v5'] else 0
        for name in impl:
            ts=times[name]
            summaries.append({'Implementation':'Strict Type-2 FxP '+name,'Trials':count,'Verified':verified[name],
                'Mean sign time (ms)':statistics.mean(ts) if ts else '',
                'Median sign time (ms)':statistics.median(ts) if ts else '',
                'Std sign time (ms)':statistics.stdev(ts) if len(ts)>1 else 0,
                'Timed completed signatures':len(ts),'Failures':count-verified[name],
                'Paired exact matches':equal,'v4/v5 speedup':mean4/mean5 if mean5 else ''})
        overall_ok &= equal==count and all(v==count for v in verified.values())
        write_csv(args.out/'campaign_summary.csv',summaries)
    print('PASS' if overall_ok else 'FAIL: inspect CSV errors/mismatches',flush=True)
    return 0 if overall_ok else 1
if __name__=='__main__': sys.exit(main())
