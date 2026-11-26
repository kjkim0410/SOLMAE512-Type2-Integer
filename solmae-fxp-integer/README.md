# SOLMAE-512 Type-2 Fixed-Point Extension

This directory contains the v4 Type-2 fixed-point implementation and preserves the reference
SOLMAE parameters, including the **40-byte salt**.

Added experimental Type-2 fixed-point modules:

- `fxp.py`, `fxp_complex.py`, `fxp_formats.py`, `fxp_math.py`
- `fft_fxp.py`
- `N_sampler_fxp.py`
- `PeikertSampler_fxp.py`
- `h_sampler_fxp.py`
- `s_keygen_fxp.py`
- `solmae_fxp.py`
- `test_fxp_core.py`, `test_fxp_solmae512.py`

The exact NTRU key generation, hashing, encoding, public-key arithmetic and
verification remain on the original/reference path. Numerical signing-key data,
FFT-domain state and the hybrid/Peikert sampler path are quantized with a
64-bit Type-2 representation (uniform word length with per-variable fractional
lengths).

`fxp_math.py` still uses quantized reference implementations of log/sin/cos for
controlled Float-vs-FxP validation. Therefore this package is a **research
prototype**, not yet a strictly floating-point-free production implementation.

Run:

```bash
python test_fxp_core.py
python test_fxp_solmae512.py
```

## Directory boundary

The v4 signing module is `../solmae-fxp-integer/solmae_fxp.py`; the v5 signing module is `../solmae-full-integer/solmae_fxp_strict.py`. Shared reference modules reside in `../solmae-ref/`. Do not copy strict modules into v4.

From the package root: `py -3.10 run_module.py v4`.
