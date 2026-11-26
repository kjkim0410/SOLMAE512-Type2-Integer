"""Hybrid SOLMAE sampler using the Type-2 fixed-point FFT path."""
from os import urandom
from fxp_formats import F
from fft_fxp import adj_fft_fxp, mul_fft_fxp, add_fft_fxp, sub_fft_fxp
from fxp_complex import CFxp
from fxp import Fxp
from PeikertSampler_fxp import PeikertSampler_fxp


def _neg_fft(v):
    return [CFxp(Fxp(-z.re.raw, z.re.frac), Fxp(-z.im.raw, z.im.frac)) for z in v]


def h_sampler_fxp(c2_fft, sk, randombytes=urandom):
    z2 = PeikertSampler_fxp(
        mul_fft_fxp(adj_fft_fxp(sk.beta21_fft), c2_fft, F["TARGET"]),
        sk.Sigma2, randombytes
    )
    v0 = mul_fft_fxp(z2, sk.b20_fft, F["TARGET"])
    v1 = mul_fft_fxp(z2, sk.b21_fft, F["TARGET"])

    temp_c1 = _neg_fft(v0)
    temp_c2 = sub_fft_fxp(c2_fft, v1, F["TARGET"])
    mu1 = add_fft_fxp(
        mul_fft_fxp(adj_fft_fxp(sk.beta10_fft), temp_c1, F["TARGET"]),
        mul_fft_fxp(adj_fft_fxp(sk.beta11_fft), temp_c2, F["TARGET"]),
        F["TARGET"]
    )

    z1 = PeikertSampler_fxp(mu1, sk.Sigma1, randombytes)
    v0 = add_fft_fxp(v0, mul_fft_fxp(z1, sk.b10_fft, F["TARGET"]), F["TARGET"])
    v1 = add_fft_fxp(v1, mul_fft_fxp(z1, sk.b11_fft, F["TARGET"]), F["TARGET"])
    return v0, v1
