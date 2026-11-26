############################################################################
# This is testing the hybrid sampler for fixed SOLMAE_D defined in params.py 
# (e.g. SOLMAE_D = 512) 
# SAFE version with defensive fixes (no dependency changes)
#############################################################################

from fft import mul_fft, adj_fft, neg, sub_fft, add_fft, ifft
from PeikertSampler_GPT import PeikertSampler
from params import SOLMAE_D, Params
from os import urandom

def ensure_vector(x, n):
    """
    Ensure x is a vector of length n.
    If scalar → expand to constant vector.
    """
    if hasattr(x, "__len__"):
        return x
    else:
        return [x] * n


def h_sampler(c2_fft, sk, randombytes=urandom):

    nu = Params[SOLMAE_D]["smoothing"]
    n = len(c2_fft)

    # ---------- First loop ----------
    mu2 = mul_fft(adj_fft(sk.beta21_fft), c2_fft)
    mu2 = ensure_vector(mu2, n)

    z2 = PeikertSampler(mu2, sk.Sigma2, randombytes, nu)
    z2 = ensure_vector(z2, n)

    v0 = mul_fft(z2, sk.b20_fft)
    v1 = mul_fft(z2, sk.b21_fft)

    v0 = ensure_vector(v0, n)
    v1 = ensure_vector(v1, n)

    temp_c1 = neg(v0)                 # -v0
    temp_c2 = sub_fft(c2_fft, v1)     # c2 - v1

    temp_c1 = ensure_vector(temp_c1, n)
    temp_c2 = ensure_vector(temp_c2, n)

    # ---------- Second loop ----------
    mu1 = add_fft(
        mul_fft(adj_fft(sk.beta10_fft), temp_c1),
        mul_fft(adj_fft(sk.beta11_fft), temp_c2)
    )
    mu1 = ensure_vector(mu1, n)

    z1 = PeikertSampler(mu1, sk.Sigma1, randombytes, nu)
    z1 = ensure_vector(z1, n)

    v0 = add_fft(v0, mul_fft(z1, sk.b10_fft))
    v1 = add_fft(v1, mul_fft(z1, sk.b11_fft))

    v0 = ensure_vector(v0, n)
    v1 = ensure_vector(v1, n)

    return v0, v1


############################################################################
# ✅ TEST HARNESS (SAFE VERSION)
############################################################################

class DummySK:
    """
    Minimal dummy secret key for testing only.
    Replace with real SOLMAE secret key.
    """
    def __init__(self, n):
        self.beta21_fft = [0]*n
        self.beta10_fft = [0]*n
        self.beta11_fft = [0]*n

        self.b20_fft = [0]*n
        self.b21_fft = [0]*n
        self.b10_fft = [0]*n
        self.b11_fft = [0]*n

        self.Sigma1 = [1.0]*n   # ensure vector
        self.Sigma2 = [1.0]*n   # ensure vector


def safe_len(x):
    return len(x) if hasattr(x, "__len__") else "scalar"


if __name__ == "__main__":

    print("\n=== h_sampler self-test for SOLMAE-", SOLMAE_D)

    n = SOLMAE_D

    # dummy inputs
    c2_fft = [0]*n
    sk = DummySK(n)

    try:
        v0, v1 = h_sampler(c2_fft, sk)

        print("h_sampler executed successfully")
        print("len(v0) =", safe_len(v0))
        print("len(v1) =", safe_len(v1))

        print("type(v0[0]) =", type(v0[0]) if hasattr(v0, "__len__") else type(v0))
        print("type(v1[0]) =", type(v1[0]) if hasattr(v1, "__len__") else type(v1))

    except Exception as e:
        print("Error:", e)

