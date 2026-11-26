# SOLMAE-512 Strict Type-2 Fixed-Point Prototype

This directory contains the v5 strict integer/fixed-point signing core.

## What changed

- `fxp_math_strict.py`: integer-only `sqrt`, `log`, `sin`, and `cos`.
  - `log`: binary normalization plus an atanh odd-power polynomial.
  - `sin/cos`: integer 2*pi range reduction and polynomial approximation on a reduced interval.
- `fft_roots_fxp_constants.py`: pre-quantized Q61 integer FFT twiddle table for n <= 512.
- `fft_fxp_strict.py`: no runtime conversion from Python complex/float twiddles.
- `N_sampler_fxp_strict.py`: N-sampler stays in `Fxp` objects; no `.to_float()` conversion.
- `samplerz_fxp_strict.py`: integer-only discrete Gaussian rejection sampler, including FACCT exp polynomial and Bernoulli test.
- `PeikertSampler_fxp_strict.py` and `h_sampler_fxp_strict.py`: fixed-point sampler pipeline with no reference FFT conversion.
- `s_keygen_fxp_strict.py`: derives the numerical Type-2 signing-key view from the exact integer NTRU basis `(f,g,F,G)` without floating-point arithmetic.
- `solmae_fxp_strict.py`: signature coefficients are rounded directly from fixed-point integers; no `.to_float()` call.

The SOLMAE salt length remains **40 bytes**, and the encoded SOLMAE-512 signature remains **666 bytes**.

## Important scope statement

The **strict Type-2 numerical signing core** and the derivation of its numerical signing-key view from an existing integer NTRU basis are floating-point-free. The original reference files are retained in the package for comparison, KAT/regression testing, public-key verification, and generation of the integer NTRU basis. Thus, importing or executing the original reference key-generation path may still use NumPy/Python floating-point code. For an experiment that isolates the fixed-point signing core, generate/load `(f,g,F,G)` first and then call `derive_fxp_key_from_basis()` followed by `sign_fxp_strict()`.

## Test

```bash
python test_fxp_strict_solmae512.py
```

Expected final line:

```text
STRICT TYPE-2 SOLMAE-512: PASS
```

## Directory boundary

The v4 signing module is `../solmae-fxp-integer/solmae_fxp.py`; the v5 signing module is `../solmae-full-integer/solmae_fxp_strict.py`. Shared reference modules reside in `../solmae-ref/`. Do not copy strict modules into v4.

From the package root: `py -3.10 run_module.py v5`.
