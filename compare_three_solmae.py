#!/usr/bin/env python3
"""Compare SOLMAE reference, Type-2 FxP v4 and full-integer v5 (Python 3.10+).

Three independently generated keypairs per trial, one per implementation.  A
shared message and salt are used for all three signers.  Different signatures
are valid: Gaussian sampling is randomized.  Benchmark is research-only.
"""
import argparse
import csv
import hashlib
import importlib
import json
import os
import platform
import statistics
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
VERSIONS = {'ref': 'solmae-ref', 'v4': 'solmae-fxp-integer', 'v5': 'solmae-full-integer'}
MODULE_NAMES = {p.stem for d in VERSIONS.values() for p in (ROOT / d).glob('*.py')}

class SaltFirstRandom:
    """Override only the first SALT_LEN bytes requested by sign().

    Subsequent draws use os.urandom. Not a reproducible sampler stream.
    """
    def __init__(self, salt):
        self.salt = salt
        self.calls = 0
    def __call__(self, n):
        self.calls += 1
        if self.calls == 1:
            if n != len(self.salt):
                raise ValueError('first random request is not SALT_LEN')
            return self.salt
        return os.urandom(n)

class Implementation:
    def __init__(self, name):
        self.name = name
        self.folder = ROOT / VERSIONS[name]
        self.modules = {}
    def __enter__(self):
        # Restore the correct module namespace for each implementation.
        self.old_path = list(sys.path)
        self.saved = {k: v for k, v in sys.modules.items() if k in MODULE_NAMES}
        for k in MODULE_NAMES: sys.modules.pop(k, None)
        sys.path[:0] = [str(self.folder), str(ROOT / 'solmae-ref')]
        try:
            self.kg = importlib.import_module('s_keygen').keygen
            self.verify = importlib.import_module('solmae').verify
            if self.name == 'ref':
                self.sign = importlib.import_module('solmae').sign
            elif self.name == 'v4':
                self.sign = importlib.import_module('solmae_fxp').sign_fxp
                self.convert = importlib.import_module('s_keygen_fxp')
            else:
                self.sign = importlib.import_module('solmae_fxp_strict').sign_fxp_strict
                self.convert = importlib.import_module('s_keygen_fxp_strict').derive_fxp_key_from_basis
        except Exception:
            self.__exit__(None, None, None)
            raise
        return self
    def __exit__(self, *_):
        for k in MODULE_NAMES: sys.modules.pop(k, None)
        sys.modules.update(self.saved)
        sys.path[:] = self.old_path
    def keygen(self):
        sk, pk = self.kg()
        if self.name == 'v4':
            # v4 keygen_fxp() invokes keygen internally. To avoid a second
            # keygen, quantize the *same* reference key using its code path.
            # v4's public API does not expose conversion separately, so use
            # the reference keygen through a temporary module binding.
            old = self.convert.keygen
            self.convert.keygen = lambda *args, **kwargs: (sk, pk)
            try: qsk, _ = self.convert.keygen_fxp()
            finally: self.convert.keygen = old
            sk = qsk
        elif self.name == 'v5':
            sk = self.convert(sk.f, sk.g, sk.F, sk.G)
        return sk, pk

def write_csv(path, rows, columns):
    tmp = path.with_suffix('.tmp')
    with tmp.open('w', newline='', encoding='utf-8-sig') as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    tmp.replace(path)

def measured(call):
    start = time.perf_counter_ns()
    value = call()
    return value, (time.perf_counter_ns() - start) / 1_000_000

def stats(vals):
    if not vals: return {'mean_ms': '', 'median_ms': '', 'min_ms': '', 'max_ms': '', 'stdev_ms': ''}
    return {'mean_ms': statistics.mean(vals), 'median_ms': statistics.median(vals),
            'min_ms': min(vals), 'max_ms': max(vals),
            'stdev_ms': statistics.stdev(vals) if len(vals) > 1 else 0.0}

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--trials', type=int, default=1000)
    ap.add_argument('--message-bytes', type=int, default=512)
    ap.add_argument('--out', type=Path, default=ROOT / 'results_three')
    ap.add_argument('--progress', type=int, default=10)
    ap.add_argument('--seed', default='SOLMAE-three-implementations-2026')
    args = ap.parse_args()
    if args.trials < 1 or args.message_bytes < 0 or args.progress < 1: ap.error('invalid argument')
    args.out.mkdir(parents=True, exist_ok=True)
    from pathlib import Path as _Path
    env = {'python': sys.version, 'python_executable': sys.executable,
           'platform': platform.platform(), 'processor': platform.processor(),
           'machine': platform.machine(), 'cpu_count': os.cpu_count(),
           'trials': args.trials, 'message_bytes': args.message_bytes,
           'seed': args.seed, 'time_unit': 'ms',
           'keygen_definition': 'reference keypair generation plus v4/v5 signing-key conversion',
           'note': 'Each implementation uses its own keypair per trial. The same message and salt are used. '
                   'Sampler randomness is not forced equal; signatures need not match. '
                   'OS randomness is used for key generation and sampler calls.'}
    try:
        import psutil
        env['cpu_model_optional'] = platform.processor()
        env['physical_cores_optional'] = psutil.cpu_count(logical=False)
    except ImportError: pass
    (args.out / 'environment.json').write_text(json.dumps(env, indent=2, ensure_ascii=False), encoding='utf-8')
    columns = ['trial','implementation','keygen_ms','sign_ms','verify_ms','verified',
               'signature_bytes','salt_matches','message_sha256','salt_hex','signature_sha256','error']
    rows = []
    pairs = []
    for i in range(args.trials):
        domain = (args.seed + ':' + str(i)).encode()
        msg = hashlib.shake_256(b'message:' + domain).digest(args.message_bytes)
        salt = hashlib.shake_256(b'salt:' + domain).digest(40)
        sigs = {}
        all_valid = True
        for name in ('ref','v4','v5'):
            row = dict.fromkeys(columns, '')
            row.update(trial=i+1, implementation=name, message_sha256=hashlib.sha256(msg).hexdigest(),
                       salt_hex=salt.hex(), verified=False, salt_matches=False)
            try:
                with Implementation(name) as impl:
                    (sk, pk), row['keygen_ms'] = measured(impl.keygen)
                    rng = SaltFirstRandom(salt)
                    sig, row['sign_ms'] = measured(lambda: impl.sign(sk, msg, randombytes=rng))
                    row['signature_bytes'] = len(sig)
                    row['signature_sha256'] = hashlib.sha256(sig).hexdigest()
                    row['salt_matches'] = sig[1:41] == salt
                    ok, row['verify_ms'] = measured(lambda: impl.verify(pk, msg, sig))
                    row['verified'] = bool(ok) and row['salt_matches'] and len(sig) == 666
                    sigs[name] = sig
            except Exception as exc:
                row['error'] = type(exc).__name__ + ': ' + str(exc)
                sigs[name] = None
            all_valid &= bool(row['verified'])
            rows.append(row)
        exact = sigs['ref'] is not None and sigs['ref'] == sigs['v4'] == sigs['v5']
        pairs.append({'trial': i+1, 'all_three_valid': all_valid, 'all_three_signatures_identical': exact,
                      'same_message': True, 'same_salt': all(bool(r['salt_matches']) for r in rows[-3:])})
        if (i+1) % args.progress == 0 or i+1 == args.trials:
            write_csv(args.out / 'compare_three_trials.csv', rows, columns)
            write_csv(args.out / 'compare_three_pairs.csv', pairs, list(pairs[0]))
            good = sum(p['all_three_valid'] for p in pairs)
            print(f'{i+1}/{args.trials}: all-three verified {good}/{i+1}, exact signatures '
                  f'{sum(p["all_three_signatures_identical"] for p in pairs)}/{i+1}', flush=True)
    summary = []
    for name in ('ref','v4','v5'):
        rs = [r for r in rows if r['implementation'] == name]
        verified = sum(bool(r['verified']) for r in rs)
        for metric in ('keygen_ms','sign_ms','verify_ms'):
            values = [r[metric] for r in rs if isinstance(r[metric], (int,float))]
            summary.append({'implementation': name, 'operation': metric.removesuffix('_ms'),
                            'trials': args.trials, 'completed': len(values),
                            'verified_signatures': verified, 'failures': args.trials-verified,
                            **stats(values)})
    write_csv(args.out / 'compare_three_summary.csv', summary, list(summary[0]))
    print('\nHardware/OS:', platform.platform(), '/', platform.processor())
    print('Python:', sys.version.split()[0], sys.executable)
    print('All three verified:', sum(p['all_three_valid'] for p in pairs), '/', args.trials)
    print('All three byte-identical:', sum(p['all_three_signatures_identical'] for p in pairs), '/', args.trials)
    print('\nMean times (ms):')
    print(f'{"Implementation":<16} {"Keygen":>14} {"Sign":>14} {"Verify":>14} {"Valid":>12}')
    for name in ('ref','v4','v5'):
        sr = [r for r in summary if r['implementation'] == name]
        means = [r['mean_ms'] for r in sr]
        fmt = lambda x: f'{x:.4f}' if isinstance(x, (int,float)) else 'N/A'
        valid = sr[0]['verified_signatures']
        print(f'{name:<16} {fmt(means[0]):>14} {fmt(means[1]):>14} {fmt(means[2]):>14} {valid:>5}/{args.trials}')
    print('Results:', args.out.resolve())
    return 0 if all(p['all_three_valid'] for p in pairs) else 1

if __name__ == '__main__':
    sys.exit(main())
