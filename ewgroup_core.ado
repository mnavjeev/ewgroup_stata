*! version 0.1.0 14may2026
program define ewgroup_core, rclass
    version 16.0
    syntax varlist(numeric min=1) [if] [in], ///
        [SIGma(varlist numeric min=1) SE(varname numeric) ///
        VARiance(varname numeric) SIGMA2(real -1) GAMMA(real -1) ///
        FACTOR(real 0.2) GENerate(name) PREFix(name) ///
        TILDEPrefix(name) WEIGHTSPrefix(name) REPLACE RETURNWeights]

    marksample touse
    markout `touse' `varlist' `sigma' `se' `variance'
    quietly count if `touse'
    if (r(N) == 0) {
        error 2000
    }
    local Jall = r(N)

    local d : word count `varlist'
    local d2 = `d' * `d'

    local input_modes = ("`sigma'" != "") + ("`se'" != "") + ("`variance'" != "")
    if (`input_modes' != 1) {
        di as err "specify exactly one of sigma(), se(), or variance()"
        exit 198
    }

    if ("`se'" != "" | "`variance'" != "") {
        if (`d' != 1) {
            di as err "se() and variance() are allowed only for scalar estimates"
            di as err "use sigma() for vector-valued estimates"
            exit 198
        }
        if (`sigma2' == -1) {
            local sigma2 1
        }
        if (`sigma2' <= 0) {
            di as err "sigma2() must be positive"
            exit 198
        }
        tempvar sigma_from_se
        if ("`se'" != "") {
            capture assert `se' >= 0 if `touse'
            if (_rc) {
                di as err "se() must contain nonnegative standard errors"
                exit 459
            }
            quietly generate double `sigma_from_se' = (`se'^2) / `sigma2' if `touse'
        }
        else {
            capture assert `variance' >= 0 if `touse'
            if (_rc) {
                di as err "variance() must contain nonnegative variances"
                exit 459
            }
            quietly generate double `sigma_from_se' = `variance' / `sigma2' if `touse'
        }
        local sigma `sigma_from_se'
    }
    else if (`sigma2' == -1) {
        di as err "sigma2() is required when using sigma()"
        exit 198
    }

    local q : word count `sigma'

    if (`d' == 1) {
        if (`q' != 1) {
            di as err "sigma() must contain one variable when beta is scalar"
            exit 198
        }
        local covtype 1
    }
    else if (`q' == `d') {
        local covtype 2
    }
    else if (`q' == `d2') {
        local covtype 3
    }
    else {
        di as err "sigma() must contain 1, d, or d*d variables"
        di as err "for d = `d', sigma() has `q' variables"
        exit 198
    }

    if (`sigma2' <= 0) {
        di as err "sigma2() must be positive"
        exit 198
    }
    if (`factor' <= 0) {
        di as err "factor() must be positive"
        exit 198
    }
    if (`gamma' != -1 & `gamma' <= 0) {
        di as err "gamma() must be positive"
        exit 198
    }

    if ("`generate'" != "" & "`prefix'" != "") {
        di as err "only one of generate() and prefix() may be specified"
        exit 198
    }
    if ("`generate'" != "" & `d' > 1) {
        di as err "generate() is allowed only for scalar beta; use prefix() for d > 1"
        exit 198
    }

    local theta_vars
    if ("`generate'" != "") {
        local theta_vars `generate'
    }
    else if ("`prefix'" != "") {
        forvalues k = 1/`d' {
            local theta_vars `theta_vars' `prefix'`k'
        }
    }

    local tilde_vars
    if ("`tildeprefix'" != "") {
        forvalues k = 1/`d' {
            local tilde_vars `tilde_vars' `tildeprefix'`k'
        }
    }

    local weight_vars
    if ("`weightsprefix'" != "") {
        quietly count if `touse'
        local J = r(N)
        forvalues k = 1/`J' {
            local weight_vars `weight_vars' `weightsprefix'`k'
        }
    }

    local protected `varlist' `sigma'
    foreach v of local theta_vars {
        if strpos(" `protected' ", " `v' ") {
            di as err "output variable `v' conflicts with an input variable"
            exit 110
        }
    }
    foreach v of local tilde_vars {
        if strpos(" `protected' ", " `v' ") {
            di as err "output variable `v' conflicts with an input variable"
            exit 110
        }
    }
    foreach v of local weight_vars {
        if strpos(" `protected' ", " `v' ") {
            di as err "output variable `v' conflicts with an input variable"
            exit 110
        }
    }

    local all_output `theta_vars' `tilde_vars' `weight_vars'
    foreach v of local all_output {
        capture confirm variable `v'
        if (_rc == 0) {
            if ("`replace'" == "") {
                di as err "`v' already exists; specify replace to overwrite"
                exit 110
            }
            capture confirm numeric variable `v'
            if (_rc) {
                di as err "`v' exists and is not numeric"
                exit 109
            }
            quietly replace `v' = .
        }
        else {
            quietly generate double `v' = .
        }
    }

    capture mata: _ewgroup_store_scalar("", 0)
    if (_rc) {
        findfile ewgroup_mata.mata
        quietly do "`r(fn)'"
    }

    tempname Theta Tilde Deriv Weights
    tempname alpha alpha_u gamma_used sure_A sure_D max_lambda
    local weights_mat
    if ("`returnweights'" != "") {
        local weights_mat `Weights'
    }

    mata: _ewgroup_core_stata("`varlist'", "`sigma'", "`touse'", ///
        `sigma2', `gamma', `factor', `covtype', "`theta_vars'", ///
        "`tilde_vars'", "`weight_vars'", "`Theta'", "`Tilde'", ///
        "`weights_mat'", "`Deriv'", "`alpha'", "`alpha_u'", ///
        "`gamma_used'", "`sure_A'", "`sure_D'", "`max_lambda'")

    return scalar alpha = `alpha'
    return scalar alpha_unconstrained = `alpha_u'
    return scalar gamma = `gamma_used'
    return scalar sigma2 = `sigma2'
    return scalar sure_A = `sure_A'
    return scalar sure_D = `sure_D'
    return scalar max_lambda = `max_lambda'
    return scalar J = `Jall'
    return scalar d = `d'
    return matrix theta = `Theta'
    return matrix tilde = `Tilde'
    return matrix derivative_trace = `Deriv'
    if ("`returnweights'" != "") {
        return matrix weights = `Weights'
    }
end
