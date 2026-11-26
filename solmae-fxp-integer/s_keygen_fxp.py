"""Type-2 fixed-point view of a SOLMAE secret key.

The exact NTRU key generation remains the reference implementation in
``s_keygen.py``.  This module quantizes the numerical FFT-domain signing-key
objects into 64-bit Type-2 fixed-point formats for the experimental FxP sampler.
"""

# Shared reference modules; keep the selected implementation first.
import sys as _solmae_sys
from pathlib import Path as _solmae_Path
_solmae_shared = str(_solmae_Path(__file__).resolve().parent.parent / "solmae-ref")
if _solmae_shared not in _solmae_sys.path:
    _solmae_sys.path.append(_solmae_shared)

from s_keygen import keygen
from fxp_complex import from_complex
from fxp import Fxp
from fxp_formats import F


class secret_key_fxp:
    """Container for the quantized numerical signing-key representation."""
    pass


def keygen_fxp(randombytes=None):
    """Generate a normal SOLMAE key pair and quantize the signing-key view."""
    sk, pk = keygen() if randombytes is None else keygen(randombytes)
    qsk = secret_key_fxp()

    complex_fields = (
        ("b10_fft", F["FFT_KEY"]),
        ("b11_fft", F["FFT_KEY"]),
        ("b20_fft", F["FFT_KEY"]),
        ("b21_fft", F["FFT_KEY"]),
        ("beta10_fft", F["BETA"]),
        ("beta11_fft", F["BETA"]),
        ("beta21_fft", F["BETA"]),
    )
    for name, frac in complex_fields:
        setattr(qsk, name, [from_complex(z, frac) for z in getattr(sk, name)])

    qsk.Sigma1 = [Fxp.from_float(x, F["SIGMA"]) for x in sk.Sigma1]
    qsk.Sigma2 = [Fxp.from_float(x, F["SIGMA"]) for x in sk.Sigma2]
    qsk.float_sk = sk
    return qsk, pk
