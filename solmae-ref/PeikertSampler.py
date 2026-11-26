#####################################################
# This is Peikert sampler for Solmae only.
#####################################################
from params import Params, SOLMAE_D #same as falcon python code
from fft import mul_fft, Solmae_ifft, fft, ifft #same as falcon python code
from N_sampler import n_sampler #coded by JH
from samplerz import samplerz #same as falcon python code
from os import urandom

def PeikertSampler(t, Sigma, randombytes=urandom,
                   nu=Params[SOLMAE_D]["smoothing"]):
    """
    Peikert Gaussian sampler (Algorithm 5 style)

    Args:
        t      : target vector (FFT domain)
        Sigma  : covariance (FFT domain)
        randombytes : randomness source
        nu     : smoothing parameter

    Returns:
        Sampled vector in FFT domain
    """

    # --- parameters (local binding for speed) ---
    sigma = Params[SOLMAE_D]["sigma"]
    n = SOLMAE_D

    # --- Step 1: sample random vector ---
    rand_vec = n_sampler(randombytes)

    # --- Step 2: compute perturbation p ---
    p = mul_fft(Sigma, rand_vec)
    p = Solmae_ifft(p)

    # --- Step 3: transform t to time domain ---
    t_time = ifft(t)

    # --- Step 4: coordinate-wise Gaussian sampling ---
    x = [
        samplerz(t_time[i] - p[i], sigma, nu)
        for i in range(n)
    ]

    # --- Step 5: return in FFT domain ---
    return fft(x)

if __name__ == '__main__':
    t = [i for i in range(SOLMAE_D)]
    Sigma = [i for i in range(SOLMAE_D)]
    print(PeikertSampler(t, Sigma))


# def PeikertSampler(t, Sigma, randombytes=urandom, nu = Params[SOLMAE_D]["smoothing"]):
#     rand_vec = n_sampler(randombytes)
#     p = mul_fft(Sigma, rand_vec)
#     p = Solmae_ifft(p)
#     t = ifft(t)
#     # x = []
#     # for i in range(SOLMAE_D):
#     #     x.append(samplerz(t[i]-p[i], Params[SOLMAE_D]["sigma"], nu))
#     x = list(map(lambda i: samplerz(t[i]-p[i], Params[SOLMAE_D]["sigma"], nu), list(range(SOLMAE_D))))
#     return fft(x)

# if __name__ == '__main__':
#     t = [i for i in range(SOLMAE_D)]
#     Sigma = [i for i in range(SOLMAE_D)]
#     print(PeikertSampler(t, Sigma))