"""Experimental Type-2 fixed-point signing front-end for SOLMAE-512.

Key generation/public-key arithmetic and verification stay exact/reference.
The numerical signing path (FFT state, key view, hybrid/Peikert sampling) uses
64-bit Type-2 fixed-point objects with per-variable fractional lengths.

This module deliberately imports HEAD_LEN and SALT_LEN from ``params.py`` so the
correct SOLMAE salt size (40 bytes in this package) is used consistently.
"""

# Shared reference modules; keep the selected implementation first.
import sys as _solmae_sys
from pathlib import Path as _solmae_Path
_solmae_shared = str(_solmae_Path(__file__).resolve().parent.parent / "solmae-ref")
if _solmae_shared not in _solmae_sys.path:
    _solmae_sys.path.append(_solmae_shared)

from os import urandom
from params import SOLMAE_D, Params, HEAD_LEN, SALT_LEN
from common import q
from encoding import compress
from solmae import hash_to_point, verify, logn
from fft_fxp import fft_fxp, ifft_fxp, sub_fft_fxp
from fxp_formats import F
from h_sampler_fxp import h_sampler_fxp


def sign_fxp(sk, message=None, randombytes=urandom, max_attempts=1000):
    if message is None:
        message = urandom(SOLMAE_D)

    header = (0x30 + logn[SOLMAE_D]).to_bytes(1, "little")
    salt = randombytes(SALT_LEN)
    hashed = hash_to_point(message, salt)

    c1_fft = fft_fxp([0] * SOLMAE_D, F["TARGET"])
    c2_fft = fft_fxp(hashed, F["TARGET"])

    for _ in range(max_attempts):
        v0, v1 = h_sampler_fxp(c2_fft, sk, randombytes)
        s1f = sub_fft_fxp(c1_fft, v0, F["TARGET"])
        s2f = sub_fft_fxp(c2_fft, v1, F["TARGET"])

        s1 = [x.to_float() for x in ifft_fxp(s1f, F["COEFF"])]
        s2 = [x.to_float() for x in ifft_fxp(s2f, F["COEFF"])]
        s1 = [(round(x) + (q >> 1)) % q - (q >> 1) for x in s1]
        s2 = [(round(x) + (q >> 1)) % q - (q >> 1) for x in s2]

        norm_sign = sum(x * x for x in s1) + sum(x * x for x in s2)
        if norm_sign <= Params[SOLMAE_D]["sig_bound"]:
            enc = compress(
                s1,
                Params[SOLMAE_D]["sig_bytelen"] - HEAD_LEN - SALT_LEN,
            )
            if enc is not False:
                return header + salt + enc

    raise RuntimeError("FxP signer exceeded max_attempts")


__all__ = ["sign_fxp", "verify"]
