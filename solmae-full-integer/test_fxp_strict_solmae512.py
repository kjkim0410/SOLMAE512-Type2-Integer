"""End-to-end validation of strict Type-2 SOLMAE-512 signing core."""

# Shared reference modules; keep the selected implementation first.
import sys as _solmae_sys
from pathlib import Path as _solmae_Path
_solmae_shared = str(_solmae_Path(__file__).resolve().parent.parent / "solmae-ref")
if _solmae_shared not in _solmae_sys.path:
    _solmae_sys.path.append(_solmae_shared)

## from s_keygen import keygen # Old and no module 
from s_keygen import keygen # new ####
from s_keygen_fxp_strict import derive_fxp_key_from_basis
from solmae_fxp_strict import sign_fxp_strict, verify

if __name__ == '__main__':
    ref_sk, pk = keygen()
    sk = derive_fxp_key_from_basis(ref_sk.f, ref_sk.g, ref_sk.F, ref_sk.G)
    for i in range(3):
        msg = (b'SOLMAE-512 strict Type-2 FxP test #' + str(i).encode())
        sig = sign_fxp_strict(sk, msg)
        ok = verify(pk, msg, sig)
        print('trial', i+1, 'sig_bytes=', len(sig), 'verify=', ok)
        assert len(sig) == 666 and ok
    print('STRICT TYPE-2 SOLMAE-512: PASS')
