# Shared reference modules; keep the selected implementation first.
import sys as _solmae_sys
from pathlib import Path as _solmae_Path
_solmae_shared = str(_solmae_Path(__file__).resolve().parent.parent / "solmae-ref")
if _solmae_shared not in _solmae_sys.path:
    _solmae_sys.path.append(_solmae_shared)

from os import urandom
from params import SOLMAE_D, SALT_LEN
from s_keygen_fxp import keygen_fxp
from solmae_fxp import sign_fxp
from solmae import verify


def main():
    assert SOLMAE_D == 512, "This prototype is tuned for SOLMAE-512"
    assert SALT_LEN == 40, "Expected corrected SOLMAE salt length of 40 bytes"
    print("Generating SOLMAE-512 key...")
    sk, pk = keygen_fxp()
    message = urandom(512)
    print("Signing with Type-2 fixed-point path...")
    sig = sign_fxp(sk, message)
    ok = verify(pk, message, sig)
    print("verify =", ok, "signature bytes =", len(sig), "salt bytes =", SALT_LEN)
    if not ok:
        raise SystemExit("FxP signature verification failed")


if __name__ == "__main__":
    main()
