"""
Average timing test for SOLMAE-512 and SOLMAE-1024
(KKJ updated version)
"""

import time
import params
from solmae import sign, verify
from s_keygen import keygen
from os import urandom

ITER = 1000
DEBUG = False


def run_benchmark(D):
    # 🔥 set parameter dynamically
    params.SOLMAE_D = D

    step1_time = 0.0  # keygen
    step2_time = 0.0  # sign
    step3_time = 0.0  # verify

    print("\n====================================")
    print(f" SOLMAE-{D} Benchmark")
    print("====================================")
    print(f"Iterations: {ITER}")

    for i in range(ITER):

        # message generation
        message = urandom(params.SOLMAE_D)

        if DEBUG:
            print("Message length:", len(message))

        if i % 10 == 0:
            print("i =", i)

        # -----------------
        # Keygen
        # -----------------
        t0 = time.perf_counter()
        sk, pk = keygen()
        t1 = time.perf_counter()

        # -----------------
        # Sign
        # -----------------
        signature = sign(sk, message)
        t2 = time.perf_counter()

        # -----------------
        # Verify
        # -----------------
        if not verify(pk, message, signature):
            print(f"\n❌ Verification failed at iteration {i} (D={D})")
            return

        t3 = time.perf_counter()

        # accumulate
        step1_time += (t1 - t0)
        step2_time += (t2 - t1)
        step3_time += (t3 - t2)

    # -----------------
    # Print results
    # -----------------
    print("\n--- Average Times ---")
    print(f"KeyGen : {step1_time / ITER * 1000:.4f} ms")
    print(f"Sign   : {step2_time / ITER * 1000:.4f} ms")
    print(f"Verify : {step3_time / ITER * 1000:.4f} ms")

    print(f"Total  : {(step1_time + step2_time + step3_time) / ITER * 1000:.4f} ms")

    print("✔ All verifications passed")


# -----------------------------
# MAIN
# -----------------------------
if __name__ == "__main__":

    print("**Avg time of SOLMAE-512 and SOLMAE-1024**")

    for D in [512, 1024]:
        run_benchmark(D)

    print("\nDone.")