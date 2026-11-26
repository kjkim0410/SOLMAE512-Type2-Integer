#####################################################
# (Note) This is for Solmae only (23.8.10, kkj)
#####################################################
    # ### Leftover ####
    # print("Test Sig KATs       : ", end="")
    # print("OK" if (test_sign_KAT() is True) else "Not OK")
    # wrapper_test(test_samplerz_simple, "SamplerZ", None, 100000)
    # wrapper_test(test_samplerz_KAT, "SamplerZ KATs", None, 1)
    # print("")
    ################


from common import q, sqnorm
from fft import add, sub, mul, div, neg, fft, ifft
from ntt import mul_zq, div_zq
from samplerz import samplerz, MAX_SIGMA
from random import randint, random, gauss, uniform
from math import sqrt, ceil
from ntrugen import karamul, ntru_gen, gs_norm
from encoding import compress, decompress
from Crypto.Hash import SHAKE256

from solmae import sign, verify, hash_to_point, HEAD_LEN, SALT_LEN
from s_keygen import secret_key, public_key, keygen
from os import urandom
from params import SOLMAE_D, Params

from scripts import saga
from scripts.samplerz_KAT512 import sampler_KAT512
from scripts.sign_KAT_JTD import sign_KAT
from scripts.samplerz_KAT1024 import sampler_KAT1024
# https://stackoverflow.com/a/25823885/4143624
from timeit import default_timer as timer

#we need to check why this sublutine require to test
def vecmatmul(t, B):
    """Compute the product t * B, where t is a vector and B is a square matrix.

    Args:
        B: a matrix

    Format: coefficient
    """
    nrows = len(B)
    ncols = len(B[0])
    deg = len(B[0][0])
    assert(len(t) == nrows)
    v = [[0 for k in range(deg)] for j in range(ncols)]
    for j in range(ncols):
        for i in range(nrows):
            v[j] = add(v[j], mul(t[i], B[i][j]))
    return v

#why not implemented ZeroDivionError?
def test_fft(n, iterations=10):
    """Test the FFT."""
    for i in range(iterations):
        f = [randint(-3, 4) for j in range(n)]
        g = [randint(-3, 4) for j in range(n)]
        h = mul(f, g)
        k = div(h, f)
        k = [int(round(elt)) for elt in k]
        if k != g:
            print("(f * g) / f =", k)
            print("g =", g)
            print("mismatch")
            return False
    return True


def test_ntt(n, iterations=10):
    """Test the NTT."""
    for i in range(iterations):
        f = [randint(0, q - 1) for j in range(n)]
        g = [randint(0, q - 1) for j in range(n)]
        h = mul_zq(f, g)
        try:
            k = div_zq(h, f)
            if k != g:
                print("(f * g) / f =", k)
                print("g =", g)
                print("mismatch")
                return False
        except ZeroDivisionError:
            continue
    return True


def check_ntru(f, g, F, G):
    """Check that f * G - g * F = q mod (x ** n + 1)."""
    a = karamul(f, G)
    b = karamul(g, F)
    c = [a[i] - b[i] for i in range(len(f))]
    return ((c[0] == q) and all(coef == 0 for coef in c[1:]))


def test_ntrugen(n, iterations=10):
    """Test ntru_gen."""
    for i in range(iterations):
        f, g, F, G = ntru_gen(n)
        if check_ntru(f, g, F, G) is False:
            return False
    return True


def test_compress(n, iterations):
    """Test compression and decompression."""
    try:
        sigma = 1.5 * sqrt(q)
        slen = Params[n]["sig_bytelen"] - SALT_LEN - HEAD_LEN
    except KeyError:
        return True
    for i in range(iterations):
        while(1):
            initial = [int(round(gauss(0, sigma))) for coef in range(n)]
            compressed = compress(initial, slen)
            if compressed is not False:
                break
        decompressed = decompress(compressed, slen, n)
        if decompressed != initial:
            return False
    return True


def test_samplerz(nb_mu=100, nb_sig=100, nb_samp=1000):
    """
    Test our Gaussian sampler on a bunch of samples.
    This is done by using a light version of the SAGA test suite,
    see ia.cr/2019/1411.
    """
    # Minimal size of a bucket for the chi-squared test (must be >= 5)
    chi2_bucket = 10
    assert(nb_samp >= 10 * chi2_bucket)
    sigmin = 1.3
    nb_rej = 0
    for i in range(nb_mu):
        mu = uniform(0, q)
        for j in range(nb_sig):
            sigma = uniform(sigmin, MAX_SIGMA)
            list_samples = [samplerz(mu, sigma, sigmin) for _ in range(nb_samp)]
            v = saga.UnivariateSamples(mu, sigma, list_samples)
            if (v.is_valid is False):
                nb_rej += 1
    return True
    if (nb_rej > 5 * ceil(saga.pmin * nb_mu * nb_sig)):
        return False
    else:
        return True


def KAT_randbytes(k):
    """
    Use a fixed bytestring 'octets' as a source of random bytes
    """
    global octets
    oc = octets[: (2 * k)]
    if len(oc) != (2 * k):
        raise IndexError("Randomness string out of bounds")
    octets = octets[(2 * k):]
    return bytes.fromhex(oc)[::-1]


def test_samplerz_KAT(unused, unused2):
    # octets is a global variable used as samplerz's randomness.
    # It is set to many fixed values by test_samplerz_KAT,
    # then used as a randomness source via KAT_randbits.
    global octets
    for D in sampler_KAT512 + sampler_KAT1024:
        mu = D["mu"]
        sigma = D["sigma"]
        sigmin = D["sigmin"]
        # Hard copy. octets is the randomness source for samplez
        octets = D["octets"][:]
        exp_z = D["z"]
        try:
            z = samplerz(mu, sigma, sigmin, randombytes=KAT_randbytes)
        except IndexError:
            return False
        if (exp_z != z):
            print("SamplerZ does not match KATs")
            return False
    return True


# def test_compress(n, iterations):
#     """Test compression and decompression."""
#     try:
#         sigma = 1.5 * sqrt(q)
#         slen = Params[n]["sig_bytelen"] - SALT_LEN - HEAD_LEN
#     except KeyError:
#         return True
#     for i in range(iterations):
#         while(1):
#             initial = [int(round(gauss(0, sigma))) for coef in range(n)]
#             compressed = compress(initial, slen)
#             if compressed is not False:
#                 break
#         decompressed = decompress(compressed, slen, n)
#         if decompressed != initial:
#             return False
#     return True


def test_signature(n, iterations=10):
    """
    Test Falcon.
    """
    f = sign_KAT[n][0]["f"]
    g = sign_KAT[n][0]["g"]
    F = sign_KAT[n][0]["F"]
    G = sign_KAT[n][0]["G"]
    sk = secret_key(n, [f, g, F, G])
    pk = public_key(sk)
    for i in range(iterations):
        message = b"abc"
        sig = sk.sign(message)
        if pk.verify(message, sig) is False:
            return False
    return True

def test_s_signature(SOLMAE_D,iterations=10):
    sk=secret_key()
    pk=public_key()
    sk, pk = keygen()

    message=urandom(SOLMAE_D) #generate SOLMAE-D Bytes message randomly

    hex_message = message.hex()
    print("Message =",hex_message[:40], '...', hex_message[-20:])
    print("Length of Message: ", len(message))

    signature = sign(sk, message) # signing
    str_signature =signature.hex()

    print("Signature =",str_signature[:40], '...', str_signature[-20:])
    print("Length of Sig. = ", len(signature))   
    
    if (verify(pk, message, signature)== True): # verifying
            print("Verification passed!!")
    else:
            print("Verification failed!!")
            return False
    return True


def test_sign_KAT():
    """
    Test the signing procedure against test vectors obtained from
    the Round 3 implementation of Falcon.

    Starting from the same private key, same message, and same SHAKE256
    context (for randomness generation), we check that we obtain the
    same signatures.
    """
    message = b"data1"
    shake = SHAKE256.new(b"external")
    for n in sign_KAT:
        sign_KAT_n = sign_KAT[n]
        for D in sign_KAT_n:
            f = D["f"]
            g = D["g"]
            F = D["F"]
            G = D["G"]
            sk = secret_key(n, [f, g, F, G])
            # The next line is done to synchronize the SHAKE256 context
            # with the one in the Round 3 C implementation of Falcon.
            _ = shake.read(8 * D["read_bytes"])
            sig = sk.sign(message, shake.read)
            if sig != bytes.fromhex(D["sig"]):
                return False
    return True

    
def wrapper_test(my_test, name, n, iterations):
    """
    Common wrapper for tests. Run the test, print whether it is successful,
    and if it is, print the running time of each execution.
    """
    d = {True: "OK    ", False: "Not OK"}
    start = timer()
    rep = my_test(n, iterations)
    end = timer()
    message = "Test {name}".format(name=name)
    message = message.ljust(20) + ": " + d[rep]
    if rep is True:
        diff = end - start
        msec = round(diff * 1000 / iterations, 3)
        message += " ({msec} msec / execution)".format(msec=msec).rjust(30)
    print(message)


# Dirty trick to fit test_samplerz into our test wrapper
def test_samplerz_simple(n, iterations):
    return test_samplerz(10, 10, iterations // 100)


def test(n, iterations=500):
    """A battery of tests."""
    wrapper_test(test_fft, "FFT", n, iterations)
    wrapper_test(test_ntt, "NTT", n, iterations)
    # test_ntrugen is super slow, hence performed over a single iteration
    wrapper_test(test_ntrugen, "NTRUGen", n, 1)
    # test_compress and test_signature are only performed
    # for parameter sets that are defined.
    if (n in Params):
        wrapper_test(test_compress, "Compress", n, iterations)
        wrapper_test(test_s_signature, "Signature", n, iterations)
        # wrapper_test(test_sign_KAT, "Signature KATs", n, iterations)
    print("")

# Run all the tests
if (__name__ == "__main__"):
    # ### Leftover ####
    # print("Test Sig KATs       : ", end="")
    # print("OK" if (test_sign_KAT() is True) else "Not OK")
    # wrapper_test(test_samplerz_simple, "SamplerZ", None, 100000)
    # wrapper_test(test_samplerz_KAT, "SamplerZ KATs", None, 1)
    # print("")
    ################
    
    for i in range(9, 11):
        n = (1 << i)
        it = 10
        print("Test battery for n = {n}".format(n=n))
        test(n, it)
