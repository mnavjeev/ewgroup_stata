# Validation of ewgroup 0.2.0

Validated on 11 September 2026 using macOS 14.6.1 (Intel), R 4.1.3,
Rcpp 1.0.8.2, testthat 3.1.4, and licensed Stata/MP 17.0.
The mathematical reference is `overleaf/draft/sections/estimation.tex`.

## Correctness and installation

- R source installation, package build, code and documentation checks, examples,
  and all tests passed. `R CMD check --no-manual --no-build-vignettes` reports one
  warning from the installed Rcpp dependency's use of deprecated `sprintf` with
  the local Xcode SDK; no package test or example failed.
- Stata original certification, the expanded regression suite, isolated local
  `net install`, and the bundled examples passed.
- Fixtures agree across installed R/C++, the pure-R fallback, an independent
  implementation of the draft equations, and Mata. Tests include finite-difference
  derivatives, scalar/diagonal/full covariances, SURE endpoints, unit scaling,
  output safety, missing data, raw QR OLS/HC0 comparisons, and Mata upgrades.
- Changes of units by factors from 1e-110 through 1e110 preserve the SURE decision
  in the regression tests. Distances or raw squared residuals outside finite
  arithmetic range now produce errors instead of silently omitted terms.
- Both source trees pass `git diff --check`. Existing edits were retained.

## Exact-estimator benchmarks

These use identical datasets in both languages. Every Stata generated estimate
matches its R counterpart to absolute tolerance 1e-8. Estimator time excludes
CSV input/output; peak resident memory includes the entire process and inputs.
Stata uses `nomatrices` and both languages omit full weights. The Stata times
below were checked using the final numerical code, including the extra checks
for overflow and extreme measurement units.

| Cells | Coefficients | R seconds | Stata seconds | R peak MiB | Stata peak MiB |
|---:|---:|---:|---:|---:|---:|
| 10,000 | 1 | 4.20 | 7.21 | 94.6 | 139.2 |
| 50,000 | 1 | 100.77 | 172.74 | 135.3 | 229.8 |
| 10,000 | 3 | 6.18 | 45.65 | 107.1 | 167.4 |

At fixed coefficient dimension, working memory grows linearly in cell count
unless full weights are requested. R tests additionally profile the fallback
with derivatives enabled and check that it does not allocate a full pairwise
matrix. Exact computation still has quadratic runtime; these timings are
machine-specific, not runtime guarantees. A 50,000 by 50,000 double weights
matrix alone requires about 20 GB.

## Compatibility and reproducibility

The declared minimum versions remain R 4.1 and Stata 16. This run directly
verified R 4.1.3 and Stata 17; Stata 16 and Windows/Linux were not available for
runtime testing. No new external Stata packages or compiled Stata plugin are
required. The R and Stata estimators retain the draft's HC0 covariance formula,
SURE recombination, and default gamma rule.

Run `python3 tests/run_checks.py --benchmarks` from `ewgroup_stata/` with the two
package repositories side by side to reproduce the checks. See [README.md](README.md)
for individual tests and fixture regeneration. The recorded numeric benchmark
results are in [benchmark_results.csv](benchmark_results.csv).
