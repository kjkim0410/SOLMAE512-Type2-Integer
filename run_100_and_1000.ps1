Set-Location $PSScriptRoot
py -3.10 test_v5_equivalence.py
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
py -3.10 benchmark_campaign_v5.py --trials 100 1000 --out results_full
exit $LASTEXITCODE
