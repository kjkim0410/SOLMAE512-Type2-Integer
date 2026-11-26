"""Type-2 fixed-point Peikert sampler adapter for SOLMAE-512.

FFT-domain targets and covariance values are stored in Type-2 fixed point.
The N-sampler geometry is quantized in ``N_sampler_fxp``.  The final scalar
``samplerz`` remains the repository's validated integer-output Gaussian sampler;
its center is converted from the quantized state at that interface.
"""

# Shared reference modules; keep the selected implementation first.
import sys as _solmae_sys
from pathlib import Path as _solmae_Path
_solmae_shared = str(_solmae_Path(__file__).resolve().parent.parent / "solmae-ref")
if _solmae_shared not in _solmae_sys.path:
    _solmae_sys.path.append(_solmae_shared)

from os import urandom
from params import Params, SOLMAE_D
from N_sampler_fxp import n_sampler_fxp
from samplerz import samplerz
from fft_fxp import to_complex_list, from_complex_list
from fxp_formats import F
from fxp import Fxp, mul
from fft import Solmae_ifft, ifft, fft


def PeikertSampler_fxp(t, Sigma, randombytes=urandom,
                       nu=Params[SOLMAE_D]["smoothing"]):
    rand = n_sampler_fxp(randombytes)

    # Sigma is the SOLMAE half-FFT real storage. Quantize each product before
    # entering the reference Solmae_ifft layout conversion.
    p = []
    for s, r in zip(Sigma, rand):
        rf = Fxp.from_float(r, F["COEFF"])
        p.append(mul(s, rf, F["COEFF"]).to_float())
    p = Solmae_ifft(p)

    t_time = ifft(to_complex_list(t))
    sigma = Params[SOLMAE_D]["sigma"]
    x = [
        samplerz(t_time[i] - p[i], sigma, nu, randombytes)
        for i in range(SOLMAE_D)
    ]
    return from_complex_list(fft(x), F["SAMPLE_FFT"])
