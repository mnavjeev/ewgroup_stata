version 16.0
clear all
set more off
set type double

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

adopath ++ "`testdir'/.."
local tol = 1e-8

di as text "checking scalar ewgroup_core fixture"
import delimited using "`testdir'/reference_scalar.csv", clear
ewgroup_core beta_hat, sigma(sigma_hat) sigma2(.05) gamma(.2) ///
    generate(theta_st) tildeprefix(tilde_st) weightsprefix(w_st) ///
    returnweights replace

assert abs(theta_st - theta) < `tol'
assert abs(tilde_st1 - tilde) < `tol'
assert abs(r(alpha) - 1) < `tol'
assert abs(r(gamma) - .2) < `tol'

forvalues k = 1/5 {
    assert abs(w_st`k' - w`k') < `tol'
}

generate double se_hat = sqrt(.05 * sigma_hat)
ewgroup_core beta_hat, se(se_hat) sigma2(.05) gamma(.2) ///
    generate(theta_se) replace
assert abs(theta_se - theta) < `tol'

generate double variance_hat = se_hat^2
ewgroup_core beta_hat, variance(variance_hat) sigma2(.05) gamma(.2) ///
    generate(theta_var) replace
assert abs(theta_var - theta) < `tol'

scalar got_alpha = r(alpha)
scalar got_alpha_unconstrained = r(alpha_unconstrained)
scalar got_gamma = r(gamma)
scalar got_sure_A = r(sure_A)
scalar got_sure_D = r(sure_D)

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

di as text "checking vector ewgroup_core fixture"
import delimited using "`testdir'/reference_vector.csv", clear
ewgroup_core b1 b2, sigma(s11 s12 s21 s22) sigma2(.1) gamma(.05) ///
    prefix(theta_st_) tildeprefix(tilde_st_) replace

assert abs(theta_st_1 - theta1) < `tol'
assert abs(theta_st_2 - theta2) < `tol'
assert abs(tilde_st_1 - tilde1) < `tol'
assert abs(tilde_st_2 - tilde2) < `tol'

scalar got_alpha = r(alpha)
scalar got_alpha_unconstrained = r(alpha_unconstrained)
scalar got_gamma = r(gamma)
scalar got_sure_A = r(sure_A)
scalar got_sure_D = r(sure_D)

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

di as text "checking raw-data ewgroup fixture"
import delimited using "`testdir'/reference_raw.csv", clear
ewgroup y, group(group) generate(theta_st) tildeprefix(tilde_st) replace

assert abs(theta_st - theta) < `tol'
assert abs(tilde_st1 - tilde) < `tol'

scalar got_alpha = e(alpha)
scalar got_alpha_unconstrained = e(alpha_unconstrained)
scalar got_gamma = e(gamma)
scalar got_sigma2 = e(sigma2)
scalar got_sure_A = e(sure_A)
scalar got_sure_D = e(sure_D)

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

di as text "checking raw-data ewgroup vector smoke test"
clear
set obs 24
generate int group = ceil(_n / 6)
generate double x = mod(_n - 1, 6) - 2.5
generate double y = 1 + .4 * group + (.1 * group) * x + ///
    cond(mod(_n, 2), .03, -.02)
ewgroup y x, group(group) prefix(theta_) tildeprefix(tilde_) replace

assert e(N) == 24
assert e(J) == 4
assert e(d) == 2
assert !missing(e(alpha))
assert !missing(theta_1)
assert !missing(theta_2)
assert !missing(tilde_1)
assert !missing(tilde_2)

di as result "ewgroup certification checks passed"
