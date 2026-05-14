*! version 0.1.0 14may2026
program define ewgroup, eclass
    version 16.0
    syntax varlist(numeric min=1) [if] [in], GROUP(varname) ///
        [NOCONstant SIGMA2(real -1) GAMMA(real -1) FACTOR(real 0.2) ///
        GENerate(name) PREFix(name) TILDEPrefix(name) REPLACE]

    marksample touse
    gettoken depvar xvars : varlist

    capture confirm numeric variable `group'
    if (_rc == 0) {
        markout `touse' `group'
        local group_use `group'
    }
    else {
        quietly replace `touse' = 0 if `group' == ""
        tempvar group_num
        quietly egen long `group_num' = group(`group') if `touse'
        local group_use `group_num'
    }

    quietly count if `touse'
    local N = r(N)
    if (`N' == 0) {
        error 2000
    }

    if ("`xvars'" == "" & "`noconstant'" != "") {
        di as err "noconstant is not allowed when no covariates are specified"
        exit 198
    }

    local d_x : word count `xvars'
    if ("`xvars'" == "") {
        local d 1
        local coefnames _cons
    }
    else if ("`noconstant'" != "") {
        local d `d_x'
        local coefnames `xvars'
    }
    else {
        local d = `d_x' + 1
        local coefnames `xvars' _cons
    }

    if (`sigma2' != -1 & `sigma2' <= 0) {
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

    local protected `varlist' `group'
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

    local all_output `theta_vars' `tilde_vars'
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

    tempname Theta Beta Tilde Cells Ncell Deriv
    tempname alpha alpha_u gamma_used sigma2_used sure_A sure_D max_lambda
    local nocons = ("`noconstant'" != "")

    mata: _ewgroup_raw_stata("`depvar'", "`xvars'", "`group_use'", ///
        "`touse'", `nocons', `sigma2', `gamma', `factor', ///
        "`theta_vars'", "`tilde_vars'", "`Theta'", "`Beta'", "`Tilde'", ///
        "`Cells'", "`Ncell'", "`Deriv'", "`alpha'", "`alpha_u'", ///
        "`gamma_used'", "`sigma2_used'", "`sure_A'", "`sure_D'", ///
        "`max_lambda'")

    matrix colnames `Theta' = `coefnames'
    matrix colnames `Beta' = `coefnames'
    matrix colnames `Tilde' = `coefnames'
    matrix colnames `Deriv' = derivative_trace
    matrix colnames `Cells' = cell
    matrix colnames `Ncell' = N

    ereturn clear
    ereturn scalar N = `N'
    ereturn scalar J = rowsof(`Theta')
    ereturn scalar d = `d'
    ereturn scalar alpha = `alpha'
    ereturn scalar alpha_unconstrained = `alpha_u'
    ereturn scalar gamma = `gamma_used'
    ereturn scalar sigma2 = `sigma2_used'
    ereturn scalar sure_A = `sure_A'
    ereturn scalar sure_D = `sure_D'
    ereturn scalar max_lambda = `max_lambda'
    ereturn matrix theta = `Theta'
    ereturn matrix beta_hat = `Beta'
    ereturn matrix tilde = `Tilde'
    ereturn matrix cells = `Cells'
    ereturn matrix ncell = `Ncell'
    ereturn matrix derivative_trace = `Deriv'
    ereturn local cmd "ewgroup"
    ereturn local depvar "`depvar'"
    ereturn local xvars "`xvars'"
    ereturn local group "`group'"
    ereturn local coefnames "`coefnames'"
    if ("`noconstant'" != "") {
        ereturn local noconstant "noconstant"
    }

    di as text "Exponentially weighted grouped estimator"
    di as text "  observations: " as result %9.0g e(N)
    di as text "  cells:        " as result %9.0g e(J)
    di as text "  dimension:    " as result %9.0g e(d)
    di as text "  gamma:        " as result %9.6g e(gamma)
    di as text "  alpha:        " as result %9.6g e(alpha)
end
