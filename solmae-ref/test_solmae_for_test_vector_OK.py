"""
Test keygen, sign, verify for SOLMAE-512 and SOLMAE-1024
"""

import params
from solmae import sign, verify
from s_keygen import secret_key, public_key, keygen
from os import urandom


# -----------------------------
# SINGLE TEST
# -----------------------------
def test_s_signature():
    sk, pk = keygen()

    print("== Leading values of keygen for SOLMAE-", params.SOLMAE_D)

    print("f =", sk.f[:5], "...")
    print("g =", sk.g[:5], "...")
    print("F =", sk.F[:5], "...")
    print("G =", sk.G[:5], "...")
    print("h =", pk.h[:5], "...")

    print("b10_fft =", sk.b10_fft[:3], "...")
    print("b11_fft =", sk.b11_fft[:3], "...")
    print("b20_fft =", sk.b20_fft[:3], "...")
    print("b21_fft =", sk.b21_fft[:3], "...")

    print("Sigma1 =", sk.Sigma1[:3], "...")
    print("Sigma2 =", sk.Sigma2[:3], "...")

    # -----------------------------
    # MESSAGE
    # -----------------------------
    message = urandom(params.SOLMAE_D)

    hex_message = message.hex()
    print("Message =", hex_message[:40], "...", hex_message[-20:])
    print("Length =", len(message))

    # -----------------------------
    # SIGN
    # -----------------------------
    signature = sign(sk, message)

    print("Signature =", signature.hex()[:40], "...")

    # -----------------------------
    # VERIFY
    # -----------------------------
    if verify(pk, message, signature):
        print("Verification passed!!")
        return True
    else:
        print("Verification failed!!")
        return False


# -----------------------------
# MAIN LOOP
# -----------------------------
if __name__ == "__main__":

    cases = 2

    print("** SOLMAE test for 512 and 1024 **")

    for D in [512, 1024]:
        print("\n==============================")
        print(f" Testing SOLMAE-{D}")
        print("==============================")

        params.SOLMAE_D = D  # 🔥 dynamic switch

        for i in range(cases):
            print(f"\n<< Test Case {i+1} >>")

            ok = test_s_signature()

            if not ok:
                print("\n❌ Test FAILED at SOLMAE-", D)
                exit(1)

    print("\n✅ ALL TESTS PASSED")