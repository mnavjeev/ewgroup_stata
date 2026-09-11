# Checking both packages

Keep the `ewgroup/` and `ewgroup_stata/` repositories beside each other. You need
R with Rcpp and testthat, Python 3, a C++ toolchain, and a licensed Stata 16 or
newer. The checks do not install into your normal R library or Stata PLUS path.

From `ewgroup_stata/`, run:

```sh
python3 tests/run_checks.py
```

The runner copies both repositories into a temporary directory, installs and
checks the R source package, regenerates reference fixtures, and runs Stata's
certification, regression, and installation tests. It prints the location of
all logs and retains it for inspection. Specify a Stata executable using
`--stata /path/to/stata-mp` or the `STATA_BIN` environment variable. Use
`--work-dir /path/to/check-output` to choose a retained output directory.

`certify.do` covers the original scalar, vector, and raw-data examples.
`regressions.do` covers the new public interface, data preservation on errors,
integer output promotion, sample selection, invalid covariance/tuning inputs,
unit changes, overflow rejection, reloads, and independently estimated raw-data
OLS/HC0 quantities. `check_net_install.do` installs the manifest into a unique
temporary PLUS directory and excludes the source checkout from command lookup.

`make_reference.R` requires an installed R package and compares every saved
fixture against both the pure-R fallback and `paper_reference.R`, a small direct
implementation of the draft's equations. The R package's own tests additionally
compare analytical derivatives with numerical finite differences. To regenerate
fixtures without touching this directory:

```sh
Rscript tests/make_reference.R --library=/path/to/R/library --output=/tmp/fixtures
```

## Benchmarks

```sh
python3 tests/run_checks.py --benchmarks
```

This also runs 10,000- and 50,000-cell scalar fits and a 10,000-cell fit with
three coefficients and full covariance matrices. Each Stata fit reads exactly
the input used by R and must match its estimates to absolute tolerance `1e-8`.
The reported estimator time excludes CSV generation/import. Peak resident
memory includes process startup and inputs; it is measured in a separate
Python process using `getrusage`, not inferred from object sizes. Results are
saved as `benchmarks.csv` in the output directory.

For benchmarks alone using an already installed R build:

```sh
python3 tests/run_checks.py --benchmark-only --library=/path/to/R/library
```

Full weights are intentionally omitted. The exact computation remains
quadratic in the number of cells, while working memory is linear at fixed
coefficient dimension. Full weights alone need about 20 GB for 50,000 cells.

See [VALIDATION.md](VALIDATION.md) for the recorded validation environment,
results, and remaining compatibility limitations.
