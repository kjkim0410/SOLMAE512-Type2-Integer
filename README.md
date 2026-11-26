# SOLMAE512-Type2-Integer

Python implementations of SOLMAE-512 using original floating-point, Type-2 fixed-point (v4), and full-integer numerical signing (v5).

Research prototype; not production cryptography. Python 3.10+.

## Layout
- `solmae-ref/`: original/reference SOLMAE modules, parameter scripts and regression tests (no FxP modules).
- `solmae-fxp-integer/`: v4 Type-2 fixed-point modules, helpers, and tests.
- `solmae-full-integer/`: v5 strict integer-path modules, helpers, and tests.
- `benchmark_campaign_v5.py`: v4/v5 paired comparison using the v5 strict algorithm with each version's arithmetic backend.

The benchmark deliberately imports the v5 strict signing/key-derivation modules for both backends to make the backend comparison meaningful. It is **not** a comparison of independent v4 and v5 end-to-end algorithm implementations.

## Signing module mapping

| Version | Module | Function |
|---|---|---|
| FP | `solmae-ref/solmae.py` | `sign` |
| v4 | `solmae-fxp-integer/solmae_fxp.py` | `sign_fxp` |
| v5 | `solmae-full-integer/solmae_fxp_strict.py` | `sign_fxp_strict` |

Strict modules belong only in v5. The v4 prototype retains floating-point helper operations; v5 targets integer-only numerical signing, not all key generation.

## Install
```powershell
py -3.10 -m pip install -r requirements.txt
```

## Run
```powershell
py -3.10 run_module.py v4
py -3.10 run_module.py v5
py -3.10 test_v5_equivalence.py
py -3.10 benchmark_campaign_v5.py --trials 100 1000 --out results_full
```

Use `run_module.py` for explicit version selection. The two end-to-end test scripts also include relative reference-path setup.

## Validation status
All 58 Python sources were syntax-checked. Packaging smoke tests on Linux / Python 3.13.5 passed: 74,096 differential comparisons, v4 signing (one trial), v5 signing (three trials), one three-way trial and one paired-backend trial. These are not new 1,000-trial campaigns. See `RELEASE_NOTES.md` for limitations and provenance.

## Three-implementation comparison (1,000 trials)

Run from this repository root (Python 3.10):

```powershell
py -3.10 -m pip install -r requirements.txt
py -3.10 compare_three_solmae.py --trials 1000 --out results_full/three
```

A quick smoke check is `py -3.10 compare_three_solmae.py --trials 1`.
The script compares the original reference, v4 fixed-point, and v5 full-integer
implementations using the same message and explicitly injected 40-byte salt per
trial. Each implementation independently generates a keypair. Keygen timing
includes the conversion/derivation of the version-specific signing-key view;
sign and verify timings are separate. The signatures are expected to verify,
not necessarily be byte-identical (Gaussian sampling is randomized). Results
are written to the selected output directory as `compare_three_trials.csv`,
`compare_three_pairs.csv`, `compare_three_summary.csv`, and `environment.json`.
All timings and hardware metadata are recorded on the computer that runs the
script. **No 1,000-trial performance numbers are asserted without execution.**

Caveat: The original reference sampler uses its own randomness in some paths;
only the first salt draw is controlled. This comparison tests common message,
salt and verification, not identical internal sampler coins. Results for key
creation include separate random keypairs, not the same shared secret key.

## Included experimental data
`results/` contains the supplied paired-backend 100/1,000-trial data. `results_three/` contains the supplied independent three-implementation 1,000-trial data. Preserve these directories; use `--out results_full/...` for new runs. The supplied three-way results report 1,000 valid signatures per implementation. No statistical equivalence of Gaussian distributions or side-channel resistance is established by these tests.

## License
No license was supplied with the source archive. Confirm the rights and intended license before public distribution; no license has been assigned by this packaging process.
