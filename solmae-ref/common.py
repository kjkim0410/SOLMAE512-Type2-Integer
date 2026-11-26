##########################################################
# (Note) Same as common.py in falcon (23.8.10, kkj)
##########################################################

"""This file contains methods and objects which are reused through multiple files such as
FFT, NTT, NTRUgen and Solmae etc. (23.8.10, kkj)"""

"""q is the integer modulus, 12289  used in Solmae-512 and Solmae-1024."""

q = 12 * 1024 + 1


def split(f):
    """Split a polynomial f in two polynomials.

    Args:
        f: a polynomial

    Format: coefficient
    """
    n = len(f)
    f0 = [f[2 * i + 0] for i in range(n // 2)]
    f1 = [f[2 * i + 1] for i in range(n // 2)]
    return [f0, f1]


def merge(f_list):
    """Merge two polynomials into a single polynomial f.

    Args:
        f_list: a list of polynomials

    Format: coefficient
    """
    f0, f1 = f_list
    n = 2 * len(f0)
    f = [0] * n
    for i in range(n // 2):
        f[2 * i + 0] = f0[i]
        f[2 * i + 1] = f1[i]
    return f


def sqnorm(v):
    """Compute the square euclidean norm of the vector v."""
    res = 0
    for elt in v:
        for coef in elt:
            res += coef ** 2
    return res
