"""Strict integer-only Peikert sampler for SOLMAE-512 Type-2 FxP."""

# Shared reference modules; keep the selected implementation first.
import sys as _solmae_sys
from pathlib import Path as _solmae_Path
_solmae_shared = str(_solmae_Path(__file__).resolve().parent.parent / "solmae-ref")
if _solmae_shared not in _solmae_sys.path:
    _solmae_sys.path.append(_solmae_shared)

from os import urandom
from params import SOLMAE_D
from N_sampler_fxp_strict import n_sampler_fxp_strict
from samplerz_fxp_strict import samplerz_fxp
from fft_fxp_strict import fft_fxp, ifft_fxp, solmae_ifft_fxp
from fxp import mul, sub
from fxp_formats import F


def PeikertSampler_fxp_strict(t, Sigma, randombytes=urandom):
    rand = n_sampler_fxp_strict(randombytes)
    p_half = [mul(s, r, F['COEFF']) for s, r in zip(Sigma, rand)]
    p_time = solmae_ifft_fxp(p_half, F['COEFF'])
    t_time = ifft_fxp(t, F['COEFF'])
    x = [samplerz_fxp(sub(t_time[i], p_time[i], F['COEFF']), randombytes)
         for i in range(SOLMAE_D)]
    return fft_fxp(x, F['SAMPLE_FFT'])
