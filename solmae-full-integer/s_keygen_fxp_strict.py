"""Integer-only derivation of the Type-2 numerical signing-key view.

Input f,g,F,G are exact integer NTRU basis polynomials.  All FFT-domain basis,
Gram-Schmidt beta values and Sigma arrays are derived without Python floats.
"""
from fxp import Fxp, sub
from fxp_complex import CFxp, cadd, csub, cmul, cdiv
from fxp_math_strict import sqrt_fxp
from fft_fxp_strict import fft_fxp, adj_fft_fxp, mul_fft_fxp, add_fft_fxp, sub_fft_fxp, div_fft_fxp
from fxp_formats import F

# Decimal-source constants prequantized offline.
SIGWIDTH_Q34 = 2981394498191       # 173.54 * 2^34
SMOOTH_Q48 = 376613518838858       # 1.338 * 2^48

class secret_key_fxp_strict:
    pass


def _real_const_fft(raw, frac, n):
    z = CFxp(Fxp(raw, frac), Fxp(0, frac))
    return [z for _ in range(n)]


def _cut_half_realimag(v, frac):
    h=len(v)//2
    return [v[i].re.rescale(frac) for i in range(h)] + [v[i].im.rescale(frac) for i in range(h)]


def derive_fxp_key_from_basis(f,g,Fpoly,Gpoly):
    n=len(f); sk=secret_key_fxp_strict()
    sk.b10_fft=fft_fxp(f,F['FFT_KEY']); sk.b11_fft=fft_fxp(g,F['FFT_KEY'])
    b1norm=add_fft_fxp(
        mul_fft_fxp(adj_fft_fxp(sk.b10_fft),sk.b10_fft,F['FFT_NORM']),
        mul_fft_fxp(adj_fft_fxp(sk.b11_fft),sk.b11_fft,F['FFT_NORM']), F['FFT_NORM'])
    sk.beta10_fft=div_fft_fxp(sk.b10_fft,b1norm,F['BETA'])
    sk.beta11_fft=div_fft_fxp(sk.b11_fft,b1norm,F['BETA'])

    # signature_width^2 and smoothing^2 in Q34.
    sw2_raw=(SIGWIDTH_Q34*SIGWIDTH_Q34) >> 34
    eta2_raw=(SMOOTH_Q48*SMOOTH_Q48) >> (96-F['FFT_NORM'])
    sw2=_real_const_fft(sw2_raw,F['FFT_NORM'],n)
    eta2=_real_const_fft(eta2_raw,F['FFT_NORM'],n)
    sig1_fft=sub_fft_fxp(div_fft_fxp(sw2,b1norm,F['FFT_NORM']),eta2,F['FFT_NORM'])
    half1=_cut_half_realimag(sig1_fft,F['SIGMA'])
    # Imaginary half is mathematically zero. Clamp tiny fixed rounding residues to zero.
    sk.Sigma1=[sqrt_fxp(x if x.raw>=0 else Fxp(0,x.frac),F['SIGMA']) for x in half1]

    sk.b20_fft=fft_fxp(Fpoly,F['FFT_KEY']); sk.b21_fft=fft_fxp(Gpoly,F['FFT_KEY'])
    temp=add_fft_fxp(
        mul_fft_fxp(adj_fft_fxp(sk.beta10_fft),sk.b20_fft,F['FFT_KEY']),
        mul_fft_fxp(adj_fft_fxp(sk.beta11_fft),sk.b21_fft,F['FFT_KEY']),F['FFT_KEY'])
    b20t=sub_fft_fxp(sk.b20_fft,mul_fft_fxp(temp,sk.b10_fft,F['FFT_KEY']),F['FFT_KEY'])
    b21t=sub_fft_fxp(sk.b21_fft,mul_fft_fxp(temp,sk.b11_fft,F['FFT_KEY']),F['FFT_KEY'])
    b2norm=add_fft_fxp(
        mul_fft_fxp(adj_fft_fxp(b20t),b20t,F['FFT_NORM']),
        mul_fft_fxp(adj_fft_fxp(b21t),b21t,F['FFT_NORM']),F['FFT_NORM'])
    sk.beta21_fft=div_fft_fxp(b21t,b2norm,F['BETA'])
    sig2_fft=sub_fft_fxp(div_fft_fxp(sw2,b2norm,F['FFT_NORM']),eta2,F['FFT_NORM'])
    half2=_cut_half_realimag(sig2_fft,F['SIGMA'])
    sk.Sigma2=[sqrt_fxp(x if x.raw>=0 else Fxp(0,x.frac),F['SIGMA']) for x in half2]
    sk.f=list(f); sk.g=list(g); sk.F=list(Fpoly); sk.G=list(Gpoly)
    return sk
