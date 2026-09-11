version 16.0
clear all
set more off
set type double
local testdir "."
capture confirm file "`testdir'/reference_cases.csv"
if (_rc) local testdir "tests"
confirm file "`testdir'/reference_cases.csv"
adopath ++ "`testdir'/.."

// Cross-language fixtures: compare scaled errors to allow different units.
import delimited using "`testdir'/reference_cases.csv", clear
levelsof case, local(cases)
foreach case of local cases {
    preserve
    keep if case == "`case'"
    local d = d[1]
    local sig2 = sigma2[1]
    local gam = gamma[1]
    local bvars
    local svars
    forvalues k = 1/`d' {
        local bvars `bvars' b`k'
        forvalues ell = 1/`d' {
            local svars `svars' s`k'`ell'
        }
    }
    quietly ewgroup `bvars', sigma(`svars') sigma2(`sig2') gamma(`gam') ///
        prefix(got_) tildeprefix(pool_)
    assert abs(e(alpha) - alpha) < 1e-8
    assert abs(e(sure_A) - sure_a) <= 1e-8 * max(abs(sure_a), 1e-300)
    assert abs(e(sure_D) - sure_d) <= 1e-8 * max(abs(sure_d), 1e-300)
    forvalues k = 1/`d' {
        quietly summarize b`k', meanonly
        local scale = max(abs(r(min)), abs(r(max)), 1e-300)
        assert abs(got_`k' - theta`k') < 1e-8 * `scale'
        assert abs(pool_`k' - tilde`k') < 1e-8 * `scale'
    }
    matrix tr = e(derivative_trace)
    forvalues j = 1/`=_N' {
        assert abs(tr[`j',1] - derivative_trace[`j']) <= ///
            1e-8 * max(abs(derivative_trace[`j']), 1e-300)
    }
    assert e(sample) == 1
    restore
}

// Summary interface, legacy compatibility, sample, and nomatrices.
clear
input double(b se)
-1.1 .200
-.95 .235
.20 .212
.27 .245
1.30 .224
end
generate long rowid = _n
generate double original_se = se
quietly ewgroup_core b, se(se) generate(legacy)
scalar original_alpha = r(alpha)
quietly ewgroup b, se(se) generate(theta) returnweights
assert abs(theta - legacy) < 1e-12
assert abs(e(alpha) - original_alpha) < 1e-12
assert e(N) == 5 & e(J) == 5 & e(d) == 1
assert e(sample) == 1
matrix summary_weights = e(weights)
assert rowsof(summary_weights) == 5
quietly ewgroup b, se(se) generate(small) nomatrices
assert abs(theta - small) < 1e-12
local matrices : e(matrices)
assert "`matrices'" == ""
assert e(sample) == 1
capture ewgroup b, se(se) nomatrices
assert _rc == 198
capture ewgroup b, se(se) generate(new) nomatrices returnweights
assert _rc == 198
capture confirm variable new
assert _rc != 0
generate double variance = se^2
quietly ewgroup b, variance(variance) generate(from_variance)
assert abs(theta - from_variance) < 1e-12

// Never overwrite any original uncertainty input, even with replace.
foreach cmd in ewgroup ewgroup_core {
    capture `cmd' b, se(se) generate(se) replace
    assert _rc == 110
    assert se == original_se
    capture `cmd' b, variance(variance) generate(variance) replace
    assert _rc == 110
    assert variance == original_se^2
    capture `cmd' b, se(se) generate(shared1) tildeprefix(shared) replace
    assert _rc == 198
    capture confirm variable shared1
    assert _rc != 0
    foreach opt in gamma sigma2 factor {
        capture `cmd' b, se(se) `opt'(.) generate(new)
        assert _rc == 198
        capture confirm variable new
        assert _rc != 0
    }
}

// Failed computation preserves outputs and previous estimation results/sample.
quietly regress b rowid
matrix before_b = e(b)
generate byte before_sample = e(sample)
generate byte existing = 99
foreach cmd in ewgroup ewgroup_core {
    capture `cmd' b, se(se) gamma(100) generate(existing) replace
    assert _rc == 198
    assert existing == 99
    local oldtype : type existing
    assert "`oldtype'" == "byte"
    assert "`e(cmd)'" == "regress"
    assert e(sample) == before_sample
    matrix after_b = e(b)
    assert mreldif(before_b, after_b) == 0
}

// Promote integer output to double and leave rows outside if/in unchanged.
quietly ewgroup b in 1/3, se(se) generate(existing) replace
local newtype : type existing
assert "`newtype'" == "double"
assert existing == 99 in 4/5
assert e(sample) == (_n <= 3)
matrix partial = e(theta)
forvalues j = 1/3 {
    assert abs(existing[`j'] - partial[`j',1]) < 1e-12
}
assert rowid == _n

// A stale or absent Mata version marker must trigger a safe reload.
capture mata: mata drop _ewgroup_version()
mata:
real scalar _ewgroup_version()
{
    return(100)
}
end
quietly ewgroup b, se(se) generate(reloaded)
mata: assert(_ewgroup_version() == 200)
assert abs(reloaded - theta) < 1e-12
mata: mata drop _ewgroup_version()
quietly ewgroup_core b, se(se)
mata: assert(_ewgroup_version() == 200)

// Missing inputs are excluded; explicit missing options are errors.
replace se = . in 2
quietly ewgroup b if rowid <= 4, se(se) generate(selected)
assert e(sample) == (rowid <= 4 & rowid != 2)
assert missing(selected) if !e(sample)
capture ewgroup b, se(se) group(rowid)
assert _rc == 198
replace se = -1 in 1
capture ewgroup b, se(se) generate(invalid)
assert _rc != 0
capture confirm variable invalid
assert _rc != 0

// Invalid covariance inputs cannot be silently symmetrized or rescaled away.
clear
input double(b1 b2 s11 s12 s21 s22)
-1 .2 1 .6 0 1.3
-.8 .1 1.1 .6 0 .9
.6 -.3 .8 .6 0 1.4
.7 -.4 1.2 .6 0 1
end
capture ewgroup b1 b2, sigma(s11 s12 s21 s22) sigma2(.1) prefix(bad_)
assert _rc != 0
capture confirm variable bad_1
assert _rc != 0
replace s12 = 0
replace s11 = -1
capture ewgroup b1 b2, sigma(s11 s12 s21 s22) sigma2(.1)
assert _rc != 0
foreach v of varlist s11 s12 s21 s22 {
    replace `v' = `v' * 1e-20
}
capture ewgroup b1 b2, sigma(s11 s12 s21 s22) sigma2(.1)
assert _rc != 0

// Finite but extreme changes of units preserve the dimensionless SURE rule.
clear
input double(b s)
-1.2 1.6
.4 1.9
-.3 .6
-.5 1.8
1 .9
-.2 1.2
end
quietly ewgroup b, sigma(s) sigma2(.1) gamma(.1) generate(base)
scalar base_alpha = e(alpha)
foreach unit in 1e-110 1e110 {
    generate double scaled_b = b * `unit'
    generate double scaled_s = s * (`unit'^2)
    local scaled_gamma = .1 / (`unit'^2)
    quietly ewgroup scaled_b, sigma(scaled_s) sigma2(.1) gamma(`scaled_gamma') generate(scaled)
    assert abs(e(alpha) - base_alpha) < 1e-10
    assert abs(scaled / `unit' - base) < 1e-10
    drop scaled_b scaled_s scaled
}

// Mata reductions must not silently omit arithmetic overflow as missing data.
clear
input double(b s)
-1e200 1
1e200 1
end
generate double untouched = 42
capture ewgroup b, sigma(s) sigma2(1) gamma(.1) generate(untouched) replace returnweights
assert _rc != 0
assert untouched == 42
clear
input double(y group)
-1e200 1
0 1
1e200 1
end
generate double untouched = 42
capture ewgroup y, group(group) gamma(.1) generate(untouched) replace
assert _rc != 0
assert untouched == 42

// Raw regression fixtures verify QR OLS and HC0 against the draft separately.
foreach mode in regression noconstant {
    import delimited using "`testdir'/reference_raw_`mode'.csv", clear
    local nc
    local d 2
    if ("`mode'" == "noconstant") {
        local nc noconstant
        local d 1
    }
    quietly ewgroup y x, group(group) `nc' prefix(got_) tildeprefix(pool_)
    assert e(sample) == 1
    assert e(N) == _N & e(J) == 4 & e(d) == `d'
    assert abs(e(alpha) - alpha) < 1e-8
    assert abs(e(sigma2) - sigma2) < 1e-12
    matrix bhat = e(beta_hat)
    levelsof group, local(groups)
    local j 0
    foreach g of local groups {
        local ++j
        forvalues k = 1/`d' {
            assert abs(bhat[`j',`k'] - beta`k') < 1e-10 if group == `g'
        }
    }
    forvalues k = 1/`d' {
        assert abs(got_`k' - theta`k') < 1e-8
        assert abs(pool_`k' - tilde`k') < 1e-8
    }
    assert row_id == _n
    generate str12 string_group = "cell_" + string(group)
    quietly ewgroup y x, group(string_group) `nc' prefix(strings_) nomatrices
    forvalues k = 1/`d' {
        assert abs(strings_`k' - theta`k') < 1e-8
    }
}

// Raw group errors must not damage data or existing estimation results.
clear
set obs 8
generate double y = sin(_n)
generate double x = _n
generate byte group = ceil(_n / 2)
generate double untouched = 12
quietly regress y x
capture ewgroup y x, group(group) prefix(untouched) replace
assert _rc != 0
assert "`e(cmd)'" == "regress"
assert untouched == 12
capture confirm variable untouched1
assert _rc != 0
replace group = ceil(_n / 4)
replace x = 1
generate double untouched1 = 17
capture ewgroup y x, group(group) prefix(untouched) replace
assert _rc == 506
assert untouched == 12
assert untouched1 == 17

display as result "ewgroup regression checks passed"
