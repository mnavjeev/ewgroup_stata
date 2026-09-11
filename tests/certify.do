version 16.0

// Start with no existing data or Mata objects from a previous run.
clear all

// Let the script run without pausing after each screen of output.
set more off

// Keep numerical comparisons precise.
set type double

// The script can be run from ewgroup_stata/ or ewgroup_stata/tests/. This block
// finds the directory that contains the reference CSV files.
local testdir "."
capture confirm file "`testdir'/reference_scalar.csv"
if (_rc) {
    local testdir "tests"
}
capture confirm file "`testdir'/reference_scalar.csv"
if (_rc) {
    di as err "run this certification script from ewgroup_stata/ or ewgroup_stata/tests/"
    exit 601
}

// Add the package root to Stata's ado search path so ewgroup and ewgroup_core
// can be found without installing them.
adopath ++ "`testdir'/.."

// Tolerance used when comparing Stata results with reference values generated
// by the R implementation.
local tol = 1e-8

// First fixture: scalar cell estimates supplied directly to ewgroup_core.
di as text "checking scalar ewgroup_core fixture"
import delimited using "`testdir'/reference_scalar.csv", clear

// `import delimited` reads the CSV fixture. ewgroup_core writes Stata estimates
// into new variables and, with returnweights, leaves the weight matrix in r().
ewgroup_core beta_hat, sigma(sigma_hat) sigma2(.05) gamma(.2) ///
    generate(theta_st) tildeprefix(tilde_st) weightsprefix(w_st) ///
    returnweights replace

// Compare generated variables and returned scalars to the reference columns.
assert abs(theta_st - theta) < `tol'
assert abs(tilde_st1 - tilde) < `tol'
assert abs(r(alpha) - 1) < `tol'
assert abs(r(gamma) - .2) < `tol'

forvalues k = 1/5 {
    assert abs(w_st`k' - w`k') < `tol'
}

// With d = 1, the default remains factor divided by the largest eigenvalue.
ewgroup_core beta_hat, sigma(sigma_hat) sigma2(.05)
assert abs(r(gamma) * r(max_lambda) - .2) < `tol'

// Check the standard-error input path. Standard errors should produce the same
// result as the scaled covariance input above.
generate double se_hat = sqrt(.05 * sigma_hat)
ewgroup_core beta_hat, se(se_hat) sigma2(.05) gamma(.2) ///
    generate(theta_se) replace
assert abs(theta_se - theta) < `tol'

// Check the variance input path as well.
generate double variance_hat = se_hat^2
ewgroup_core beta_hat, variance(variance_hat) sigma2(.05) gamma(.2) ///
    generate(theta_var) replace
assert abs(theta_var - theta) < `tol'

scalar got_alpha = r(alpha)
scalar got_alpha_unconstrained = r(alpha_unconstrained)
scalar got_gamma = r(gamma)
scalar got_sure_A = r(sure_A)
scalar got_sure_D = r(sure_D)

// Load expected scalar return values from a small CSV. `preserve` and `restore`
// temporarily swap in that CSV without losing the main fixture data.
preserve
import delimited using "`testdir'/reference_scalar_scalars.csv", clear
foreach s in alpha alpha_unconstrained gamma sure_A sure_D {
    quietly summarize value if name == "`s'", meanonly
    scalar exp_`s' = r(mean)
}
restore

foreach s in alpha alpha_unconstrained gamma sure_A sure_D {
    assert abs(got_`s' - exp_`s') < `tol'
}

// Second fixture: vector-valued cell estimates supplied directly to ewgroup_core.
di as text "checking vector ewgroup_core fixture"
import delimited using "`testdir'/reference_vector.csv", clear

// sigma() receives all entries of the two-by-two covariance matrix for each
// cell: s11, s12, s21, and s22.
ewgroup_core b1 b2, sigma(s11 s12 s21 s22) sigma2(.1) gamma(.05) ///
    prefix(theta_st_) tildeprefix(tilde_st_) replace

// Compare both final coefficients and both intermediate smoothed coefficients.
assert abs(theta_st_1 - theta1) < `tol'
assert abs(theta_st_2 - theta2) < `tol'
assert abs(tilde_st_1 - tilde1) < `tol'
assert abs(tilde_st_2 - tilde2) < `tol'

scalar got_alpha = r(alpha)
scalar got_alpha_unconstrained = r(alpha_unconstrained)
scalar got_gamma = r(gamma)
scalar got_sure_A = r(sure_A)
scalar got_sure_D = r(sure_D)

// Compare scalar diagnostics for the vector fixture.
preserve
import delimited using "`testdir'/reference_vector_scalars.csv", clear
foreach s in alpha alpha_unconstrained gamma sure_A sure_D {
    quietly summarize value if name == "`s'", meanonly
    scalar exp_`s' = r(mean)
}
restore

foreach s in alpha alpha_unconstrained gamma sure_A sure_D {
    assert abs(got_`s' - exp_`s') < `tol'
}

// For vector estimates, the default also divides by their dimension d.
ewgroup_core b1 b2, sigma(s11 s12 s21 s22) sigma2(.1)
assert abs(r(gamma) * 2 * r(max_lambda) - .2) < `tol'

// Third fixture: raw observations, where ewgroup first estimates one mean per
// group before running the smoothing step.
di as text "checking raw-data ewgroup fixture"
import delimited using "`testdir'/reference_raw.csv", clear
ewgroup y, group(group) generate(theta_st) tildeprefix(tilde_st) replace

// Generated observation-level variables should match the reference values.
assert abs(theta_st - theta) < `tol'
assert abs(tilde_st1 - tilde) < `tol'

scalar got_alpha = e(alpha)
scalar got_alpha_unconstrained = e(alpha_unconstrained)
scalar got_gamma = e(gamma)
scalar got_sigma2 = e(sigma2)
scalar got_sure_A = e(sure_A)
scalar got_sure_D = e(sure_D)

// Compare e() scalar diagnostics for the raw-data fixture.
preserve
import delimited using "`testdir'/reference_raw_scalars.csv", clear
foreach s in alpha alpha_unconstrained gamma sigma2 sure_A sure_D {
    quietly summarize value if name == "`s'", meanonly
    scalar exp_`s' = r(mean)
}
restore

foreach s in alpha alpha_unconstrained gamma sigma2 sure_A sure_D {
    assert abs(got_`s' - exp_`s') < `tol'
}

// Final smoke test: raw observations with one covariate. This checks that the
// command can estimate and store both a slope and an intercept.
di as text "checking raw-data ewgroup vector smoke test"
clear
set obs 24
generate int group = ceil(_n / 6)
generate double x = mod(_n - 1, 6) - 2.5
generate double y = 1 + .4 * group + (.1 * group) * x + ///
    cond(mod(_n, 2), .03, -.02)
ewgroup y x, group(group) prefix(theta_) tildeprefix(tilde_) replace

// Confirm expected result sizes and verify generated variables are not missing.
assert e(N) == 24
assert e(J) == 4
assert e(d) == 2
assert !missing(e(alpha))
assert !missing(theta_1)
assert !missing(theta_2)
assert !missing(tilde_1)
assert !missing(tilde_2)

di as result "ewgroup certification checks passed"
