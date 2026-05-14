# ewgroup for Stata

`ewgroup` is a Stata package for estimating many noisy group- or cell-level
parameters. It is useful when you have one estimate per unit, such as a teacher,
firm, plan, judge, occupation, county, or other discrete group, and you think
some units may have similar underlying effects.

The package starts from the usual cell-by-cell estimates. It then averages
nearby estimates using weights based on both the estimates and their sampling
uncertainty. A final SURE step decides how much of that averaging to use. If the
data do not support much pooling, the final estimate can stay close to the
cell-by-cell estimate.

## Install

From Stata:

```stata
net install ewgroup, from("https://raw.githubusercontent.com/mnavjeev/ewgroup_stata/main") replace
```

Then check that Stata can find the command:

```stata
help ewgroup
help ewgroup_core
```

If you are working from a local checkout of the repository instead of GitHub,
use:

```stata
adopath ++ "/path/to/ewgroup_stata"
```

## Which Command Should I Use?

Use `ewgroup_core` if you already have one row per unit with an estimate and a
standard error. This is the closest workflow to empirical Bayes commands.

```stata
ewgroup_core beta_hat, se(se_hat) generate(theta)
```

Use `ewgroup` if you have the raw data and want Stata to compute the first-stage
cell estimates for you.

```stata
ewgroup y, group(unit) generate(theta)
```

Use `ewgroup` with regressors if each unit has its own regression coefficient:

```stata
ewgroup y x1 x2, group(unit) prefix(theta_)
```

Use the `sigma()` option only when you are working with covariance matrices
directly, especially for vector-valued estimates.

## What Should the Data Look Like?

There are two common data layouts.

### Raw data layout

Use this layout with `ewgroup`. Each row is one observation. The group variable
tells Stata which cell the observation belongs to.

```text
w       y
1      -1.25
1      -1.10
1      -1.05
2      -0.85
2      -0.78
3       0.35
3       0.42
```

For a regression with cell-specific slopes, the raw data also include the
regressors:

```text
w       y        x
1      1.21     0.0
1      1.28     0.1
1      1.31     0.2
2      1.40     0.0
2      1.51     0.1
2      1.59     0.2
```

In this layout, use `ewgroup`.

### Estimate-and-standard-error layout

Use this layout with `ewgroup_core`. Each row is one unit or cell, not one raw
observation.

```text
unit    beta_hat    se_hat
1       -1.10       0.200
2       -0.95       0.235
3        0.20       0.212
```

In this layout, use `ewgroup_core`.

## Step-by-Step: Estimates and Standard Errors

This is the simplest use case. Suppose you have one estimate per unit and a
standard error for each estimate.

Your data should look like this:

```text
unit    beta_hat    se_hat
1       -1.10       0.200
2       -0.95       0.235
3        0.20       0.212
4        0.27       0.245
5        1.30       0.224
```

If those variables are already in memory, run:

```stata
ewgroup_core beta_hat, se(se_hat) generate(theta) replace
list unit beta_hat se_hat theta
```

Here:

- `beta_hat` is the original cell-by-cell estimate.
- `se_hat` is its standard error.
- `theta` is the adjusted estimate produced by `ewgroup_core`.

If you have variances rather than standard errors:

```stata
generate double variance_hat = se_hat^2
ewgroup_core beta_hat, variance(variance_hat) generate(theta2) replace
```

## Step-by-Step: Cell Means From Raw Data

Use this section when your dataset has one row per observation, not one row per
cell. The outcome is `y`, and `w` says which cell the observation belongs to.

The model is:

```text
y_i = beta(w_i) + error_i
```

In this example there are four cells. Cell 1 has five observations, cell 2 has
five observations, and so on. Your data should be in long form:

```text
observation    w       y
1              1      -1.25
2              1      -1.10
3              1      -1.05
4              1      -0.95
5              1      -0.90
6              2      -0.85
7              2      -0.78
...            ...     ...
18             4       1.22
19             4       1.28
20             4       1.35
```

Enter the example data:

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
```

Estimate the adjusted cell means:

```stata
ewgroup y, group(w) generate(theta) replace
```

This command does three things:

1. It computes the usual mean of `y` in each cell `w`.
2. It estimates the sampling variance of each cell mean from the observations
   inside the cell.
3. It applies the `ewgroup` adjustment and writes the adjusted cell mean to
   `theta`.

For example, before the adjustment, the command treats the five observations
with `w == 1` as the data for cell 1, the five observations with `w == 2` as the
data for cell 2, and so on. You do not need to collapse the data first.

Now look at the result:

```stata
list w y theta, sepby(w)
matrix list e(theta)
```

The variable `theta` is constant within each value of `w`, because the estimate
is a cell-level estimate. The matrix `e(theta)` gives the same estimates with
one row per cell. The row order follows the cell ids in `e(cells)`.

## Step-by-Step: Cell-Specific Regressions

Use this when each cell has its own slope, intercept, or vector of regression
coefficients:

```text
y_i = beta(w_i)' x_i + error_i
```

Example with one regressor and an intercept:

```text
observation    w       y        x
1              1      1.24     0.0
2              1      1.20     0.1
3              1      1.27     0.2
...            ...    ...      ...
21             2      1.41     0.0
22             2      1.44     0.1
```

Each cell must have enough observations to estimate its own regression. With
one regressor and an intercept, each cell needs at least two observations, and
the regressor must vary within the cell.

If those variables are already in memory, run:

```stata
ewgroup y x, group(w) prefix(theta_) replace
describe theta_1 theta_2
matrix list e(beta_hat)
matrix list e(theta)
```

By default Stata includes an intercept. With one regressor, `theta_1` is the
slope on `x` and `theta_2` is the intercept. The column names of `e(theta)` show
the same order.

To omit the intercept:

```stata
ewgroup y x, group(w) noconstant prefix(theta_) replace
```

## Reading the Main Results

For `ewgroup_core`, the adjusted estimates are stored in `r(theta)`.

```stata
matrix list r(theta)
return list
```

For `ewgroup`, the adjusted estimates are stored in `e(theta)`.

```stata
matrix list e(theta)
ereturn list
```

The most useful scalars are:

- `alpha`: the SURE weight on the averaged estimate.
- `gamma`: the tuning parameter used in the exponential weights.
- `sigma2`: the scale parameter used for the covariance estimates.

When `alpha` is close to zero, the final estimate is close to the original
cell-by-cell estimate. When `alpha` is close to one, the final estimate puts
more weight on the averaged estimate.

## Vector-Valued Estimates With Covariance Matrices

For scalar estimates, use `se()` or `variance()` unless you have a reason to use
the paper's scaled covariance notation.

For vector-valued estimates, each unit has more than one coefficient. For
example, an occupation might have a baseline wage premium and a wage-growth
coefficient, or a region might have an intercept and a slope.

If you are starting from raw data, use `ewgroup`:

```stata
ewgroup y x, group(w) prefix(theta_) replace
```

Stata will compute the cell-specific regression estimates and their covariance
matrices for you.

If you already have the first-stage estimates and covariance matrices, use
`ewgroup_core`. The data should have one row per unit:

```text
unit      b1       b2       s11      s12      s21      s22
1       -1.00     0.20      1.00     0.10     0.10     1.30
2       -0.80     0.10      1.10     0.00     0.00     0.90
3        0.60    -0.30      0.80     0.20     0.20     1.40
4        0.70    -0.40      1.20     0.10     0.10     1.00
```

Here `b1` and `b2` are the two preliminary coefficients. The variables `s11`,
`s12`, `s21`, and `s22` are the entries of the covariance matrix for those two
coefficients.

### Diagonal covariance matrices

If you only want to use the variances of each coefficient and ignore the
covariance between coefficients, pass the diagonal entries to `sigma()`:

```text
unit      b1       b2       s11      s22
1       -1.00     0.20      1.00     1.30
2       -0.80     0.10      1.10     0.90
3        0.60    -0.30      0.80     1.40
```

```stata
ewgroup_core b1 b2, sigma(s11 s22) sigma2(.1) gamma(.05) prefix(theta_) replace
```

This writes `theta_1` and `theta_2`, the adjusted versions of `b1` and `b2`.

### Full covariance matrices

If you have the full covariance matrix, list the entries in row order. For a
two-coefficient estimate, the order is:

```text
s11 s12 s21 s22
```

```stata
ewgroup_core b1 b2, sigma(s11 s12 s21 s22) sigma2(.1) gamma(.05) ///
    prefix(theta_) replace
```

For three coefficients, the order would be:

```text
s11 s12 s13 s21 s22 s23 s31 s32 s33
```

The covariance matrix should be symmetric. In practice, `s12` and `s21` should
usually be the same up to rounding.

With `sigma()`, the variance of the preliminary estimate is interpreted as:

```text
Var(beta_hat_j) = sigma2 * Sigma_hat_j
```

This is the notation used in the paper and in the companion R package.

If your covariance variables are ordinary variances and covariances, a simple
choice is `sigma2(1)`. If you are using the normalization from the paper, pass
the corresponding `sigma2()` value.

## Tuning Parameters

Most users can start with the defaults.

For `ewgroup`, the default `sigma2` is `J/N`, where `J` is the number of cells
and `N` is the number of observations.

For `ewgroup_core` with `se()` or `variance()`, the default is `sigma2(1)`,
because the standard errors or variances are already on the usual scale.

If `gamma()` is not supplied, the package uses:

```text
gamma = factor / max_j lambda_max(Sigma_hat_j)
```

with `factor(0.2)` by default. A larger `gamma` makes the weights more sensitive
to differences between preliminary estimates. The command checks that `gamma`
is small enough for the calculation to be well defined.

## Common Mistakes

### The output variable already exists

Use `replace` if you want to overwrite generated output variables:

```stata
ewgroup_core beta_hat, se(se_hat) generate(theta) replace
```

### I have raw data but used `ewgroup_core`

`ewgroup_core` expects one row per cell. If your data have one row per person,
case, observation, or transaction, use `ewgroup`.

### I have estimates and standard errors but used `ewgroup`

Use `ewgroup_core`:

```stata
ewgroup_core beta_hat, se(se_hat) generate(theta)
```

### My group variable is a string

`ewgroup` accepts string group variables. It encodes them internally before
estimation. The stored `e(cells)` matrix contains the encoded numeric ids.

## Running the Example File

The repository includes a short example script:

```stata
cd "/path/to/ewgroup_stata"
do examples/basic_usage.do
```

## Testing

The certification script compares the Stata output against fixtures generated
from the companion R package:

```stata
cd "/path/to/ewgroup_stata"
do tests/certify.do
```

The package manifest can also be checked with a local `net install`:

```stata
do tests/check_net_install.do
```

To regenerate the R fixtures:

```bash
Rscript tests/make_reference.R
```

Then rerun the Stata certification script.

## Files in This Repository

```text
ewgroup.ado              raw-data command
ewgroup_core.ado         cell-summary command
ewgroup_mata.mata        shared Mata code
ewgroup.sthlp            Stata help for ewgroup
ewgroup_core.sthlp       Stata help for ewgroup_core
ewgroup.pkg              net install package manifest
stata.toc                Stata package table of contents
examples/                runnable example do-files
tests/                   certification tests and reference fixtures
```

## Requirements

- Stata 16 or newer.
- No external Stata packages are required.
- R is only needed if you want to regenerate the reference fixtures in
  `tests/`.

## Authors

Manu Navjeevan, Denis Chetverikov, and Andrei Voronin.
