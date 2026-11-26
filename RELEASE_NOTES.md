# SOLMAE512-Type2-Integer: packaging and validation notes

Final archive naming and documentation cleanup: removed obsolete development-package names. Source code and experimental data remain unchanged. The execution checks below refer to the preceding packaging pass.

All 58 Python files are preserved byte-for-byte, including legacy reference KATs and parameter scripts. Version-local arithmetic modules are intentionally separate. Removed superseded documentation and historical smoke/validation outputs; retained both current experimental datasets. Removed the local Python executable path from results_three/environment.json. No numerical changes were made.

## New checks
- Python syntax: 58/58.
- Differential arithmetic: 74,096 comparisons passed.
- v4 end-to-end: one valid 666-byte signature, 40-byte salt.
- v5 end-to-end: three valid 666-byte signatures.
- Independent three-way smoke test: 1/1 all valid, 0/1 byte-identical.
- Paired backend smoke test: 1/1 both valid and equal.

Environment: Linux, Python 3.13.5. Windows/Python 3.10, full campaigns and legacy reference KAT scripts were not rerun. Historical campaign source hashes match the packaged files except test_fxp_solmae512.py; original metadata is retained. The three-way data does not contain source hashes. No license was supplied. This is not a security audit.

## Removed original entries
- `PATH_FIX_REPORT.md`
- `UPDATE_REPORT.md`
- `VALIDATION.md`
- `STRUCTURE_MANIFEST.md`
- `README_v5_KO.md`
- `solmae-full-integer/README_Strict_Type2_FxP.md`
- `solmae-fxp-integer/README_Type2_FxP.md`
- `results_smoke`
- `validation_100`
