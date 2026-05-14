# ewgroup for Stata

`ewgroup` implements an exponentially weighted grouped-heterogeneity estimator
for Stata. The estimator starts from noisy cell-specific estimates, pools nearby
cells using covariance-adjusted exponential weights, applies the debiasing
correction from the paper, and combines the pooled estimate with the naive
cell-by-cell estimate using a feasible SURE rule.

The package has two user-facing commands:

- `ewgroup`: starts from raw data, estimates cell-by-cell least-squares
  coefficients, computes scaled HC0 covariance estimates, and applies the
  estimator.
- `ewgroup_core`: starts from one row per cell with precomputed naive estimates
  and standard errors, variances, or scaled covariance estimates. For scalar
  estimates this is close to the workflow used by empirical Bayes commands:
  give the command a point estimate and its standard error, and ask for an
  adjusted estimate.

## Repository Layout

```text
ewgroup.ado              raw-data command
ewgroup_core.ado         cell-summary command
ewgroup_mata.mata        shared Mata implementation
ewgroup.sthlp            Stata help for ewgroup
ewgroup_core.sthlp       Stata help for ewgroup_core
ewgroup.pkg              net install package manifest
stata.toc                Stata package table of contents
examples/                runnable example do-files
tests/                   certification tests and R reference fixtures
```

## Requirements

- Stata 16 or newer.
- No external Stata packages are required.
- R is only needed if you want to regenerate the reference fixtures in
  `tests/` from the companion R implementation.

## Installation

### Local development checkout

From Stata:

```stata
adopath ++ "/path/to/ewgroup_stata"
help ewgroup
```

This is the easiest setup while editing the package.

### Local net install

From Stata:

```stata
net install ewgroup, from("/path/to/ewgroup_stata") replace
help ewgroup
```

### Public GitHub install

After pushing this directory to a public repository, users can install from the
raw GitHub URL. For example, if the repository is
`https://github.com/USER/ewgroup_stata`, users can run:

```stata
net install ewgroup, from("https://raw.githubusercontent.com/USER/ewgroup_stata/main") replace
```

If the package is kept in a subdirectory of a larger repository, the `from()`
URL should point to the directory containing `stata.toc` and `ewgroup.pkg`.

## Quick Start

### Empirical-Bayes-style use with estimates and standard errors

Many applied users will already have one estimate per teacher, firm, plan,
judge, region, or other unit, together with a standard error. In that case use
`ewgroup_core` with `se()`:

```stata
clear
input double beta_hat se_hat
-1.10 .200
-0.95 .235
 0.20 .212
 0.27 .245
 1.30 .224
end

ewgroup_core beta_hat, se(se_hat) generate(theta) replace
list beta_hat se_hat theta
```

Here `beta_hat` is the cell-by-cell estimate and `theta` is the final adjusted
estimate. This is the simplest workflow when the first-stage estimates and
standard errors have already been computed.

If you have variances instead of standard errors, use `variance()`:

```stata
generate double variance_hat = se_hat^2
ewgroup_core beta_hat, variance(variance_hat) generate(theta2) replace
```

The `se()` and `variance()` options are for scalar estimates. For vector-valued
coefficients, use the `sigma()` covariance input described below.

### Cell means from raw data

Use this when the model is

```text
y_i = beta(w_i) + error_i
```

and you want one final estimate per cell.

```stata
clear
input byte w double y
1 -1.25
1 -1.10
1 -1.05
1 -0.95
1 -0.90
2 -0.85
2 -0.78
2 -0.72
2 -0.70
2 -0.65
3  0.35
3  0.42
3  0.55
3  0.60
3  0.63
4  1.05
4  1.10
4  1.22
4  1.28
4  1.35
end

ewgroup y, group(w) generate(theta) tildeprefix(tilde_) replace
list w y theta tilde_1, sepby(w)
matrix list e(theta)
```

The generated variable `theta` repeats the final cell-level estimate for each
observation in the same cell. The matrix `e(theta)` contains one row per cell.

### Cell-specific regressions from raw data

Use this when the model is

```text
y_i = beta(w_i)' x_i + error_i
```

The command includes an intercept by default, following Stata's regression
convention.

```stata
clear
set obs 80
generate byte w = ceil(_n / 20)
generate double x = mod(_n - 1, 20) / 10
generate double y = 1 + .2*w + (.1*w)*x + sin(_n)/20

ewgroup y x, group(w) prefix(theta_) tildeprefix(tilde_) replace

describe theta_1 theta_2
matrix list e(beta_hat)
matrix list e(theta)
```

With one regressor and an intercept, `theta_1` corresponds to the slope on `x`
and `theta_2` corresponds to the intercept. The column names of `e(theta)` show
the same coefficient order.

To omit the intercept:

```stata
ewgroup y x, group(w) noconstant prefix(theta_) replace
```

## Low-Level Interface: `ewgroup_core`

Use `ewgroup_core` when you already have one row per cell:

- `beta_vars`: naive estimates, one column per coefficient.
- `sigma()`: scaled covariance estimates.
- `sigma2()`: scale parameter such that
  `Var(beta_hat_j) = sigma2 * Sigma_hat_j`.

### Scalar case

```stata
clear
input double beta_hat se_hat
-1.10 .200
-0.95 .235
 0.20 .212
 0.27 .245
 1.30 .224
end

ewgroup_core beta_hat, se(se_hat) generate(theta) tildeprefix(tilde_) ///
    returnweights replace

list beta_hat se_hat theta tilde_1
matrix list r(weights)
```

If you want to work with the paper's scaled covariance notation directly, pass
`sigma()` and `sigma2()`. In that case `sigma()` contains `Sigma_hat_j` and the
variance of `beta_hat_j` is `sigma2 * Sigma_hat_j`.

### Vector case with diagonal covariance estimates

If there are `d` coefficient columns and `sigma()` contains `d` variables,
`ewgroup_core` treats those variables as the diagonal entries of each
cell-specific covariance matrix.

```stata
ewgroup_core b1 b2, sigma(s11 s22) sigma2(.1) gamma(.05) prefix(theta_) replace
```

### Vector case with full covariance estimates

If `sigma()` contains `d*d` variables, the variables are interpreted in
row-major order:

```text
S11 S12 ... S1d S21 S22 ... Sdd
```

For `d = 2`:

```stata
ewgroup_core b1 b2, sigma(s11 s12 s21 s22) sigma2(.1) gamma(.05) ///
    prefix(theta_) tildeprefix(tilde_) replace
```

## Tuning Parameters

### `sigma2()`

For `ewgroup`, the default is `J/N`, where `J` is the number of observed cells
and `N` is the estimation sample size. This matches the normalization used in
the paper.

For `ewgroup_core`, `sigma2()` is required because the command does not know how
the preliminary estimates were constructed.

### `gamma()`

`gamma()` controls the exponential weighting. If omitted, the package uses

```text
gamma = factor / max_j lambda_max(Sigma_hat_j)
```

where `factor` defaults to `0.2`.

The command requires

```text
gamma * max_j lambda_max(Sigma_hat_j) < 1
```

so that the covariance-adjustment matrices are invertible.

## Stored Results

### `ewgroup`

`ewgroup` is an e-class command. The main stored results are:

- `e(theta)`: final cell-level estimates.
- `e(beta_hat)`: naive cell-level estimates.
- `e(tilde)`: exponentially weighted estimates before SURE recombination.
- `e(cells)`: cell ids corresponding to rows of the matrices.
- `e(ncell)`: cell sample sizes.
- `e(alpha)`: SURE mixing weight.
- `e(gamma)`: tuning parameter used.
- `e(sigma2)`: scale parameter used.

### `ewgroup_core`

`ewgroup_core` is an r-class command. The main stored results are:

- `r(theta)`: final estimates.
- `r(tilde)`: exponentially weighted estimates before SURE recombination.
- `r(weights)`: exponential weights, if `returnweights` is specified.
- `r(alpha)`: SURE mixing weight.
- `r(gamma)`: tuning parameter used.
- `r(sigma2)`: scale parameter used.

See `help ewgroup` and `help ewgroup_core` for the full list.

## Running Examples

The example script can be run from Stata:

```stata
cd "/path/to/ewgroup_stata"
do examples/basic_usage.do
```

## Testing

The certification script compares Stata results against fixtures generated from
the companion R package and also runs a raw-data vector regression smoke test.

```stata
cd "/path/to/ewgroup_stata"
do tests/certify.do
```

To regenerate the R fixtures:

```bash
Rscript tests/make_reference.R
```

Then rerun the Stata certification script.

## Notes for Public Release

Before tagging a public release:

- Run `do tests/certify.do` in Stata.
- Check that `ewgroup.pkg` lists every file that should be installed.
- Replace the example GitHub URL in this README with the real repository URL.
- Add a license file if the package is released independently of the companion
  R package.

## Authors

Manu Navjeevan, Denis Chetverikov, and Andrei Voronin.
