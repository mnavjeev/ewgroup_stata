*! version 0.2.0 11sep2026
// Retain the low-level r() interface, with rollback on errors and interrupts.
program define ewgroup_core, rclass
    version 16.0
    preserve
    _ewgroup_core_work `0'
    nobreak {
        restore, not
        return add
    }
end

program define _ewgroup_core_work, rclass
    version 16.0
    syntax varlist(numeric min=1) [if] [in], ///
        [SIGma(varlist numeric min=1) SE(varname numeric) ///
        VARiance(varname numeric) SIGMA2(string) GAMMA(string) ///
        FACTOR(real 0.2) GENerate(name) PREFix(name) ///
        TILDEPrefix(name) WEIGHTSPrefix(name) REPLACE RETURNWeights NOMATRices]

    local protected `varlist' `sigma' `se' `variance'
    local input_modes = ("`sigma'" != "") + ("`se'" != "") + ("`variance'" != "")
    if `input_modes' != 1 {
        di as err "specify exactly one of sigma(), se(), or variance()"
        exit 198
    }
    local d : word count `varlist'
    local d2 = `d' * `d'
    if ("`se'" != "" | "`variance'" != "") & `d' != 1 {
        di as err "se() and variance() require scalar estimates; use sigma() for vectors"
        exit 198
    }
    if "`sigma2'" == "" {
        if "`sigma'" != "" {
            di as err "sigma2() is required when using sigma()"
            exit 198
        }
        local sigma2 1
    }
    capture confirm number `sigma2'
    if _rc {
        di as err "sigma2() must be finite and positive"
        exit 198
    }
    if missing(`sigma2') | `sigma2' <= 0 {
        di as err "sigma2() must be finite and positive"
        exit 198
    }
    if "`gamma'" == "" {
        local gamma .
    }
    else {
        capture confirm number `gamma'
        if _rc {
            di as err "gamma() must be finite and positive"
            exit 198
        }
        if missing(`gamma') | `gamma' <= 0 {
            di as err "gamma() must be finite and positive"
            exit 198
        }
    }
    if missing(`factor') | `factor' <= 0 {
        di as err "factor() must be finite and positive"
        exit 198
    }

    marksample touse
    markout `touse' `sigma' `se' `variance'
    quietly count if `touse'
    local Jall = r(N)
    if `Jall' == 0 error 2000

    if "`se'" != "" | "`variance'" != "" {
        tempvar sigma_from_se
        if "`se'" != "" {
            capture assert `se' >= 0 if `touse'
            if _rc {
                di as err "se() must contain nonnegative standard errors"
                exit 459
            }
            quietly generate double `sigma_from_se' = (`se' / sqrt(`sigma2'))^2 if `touse'
        }
        else {
            capture assert `variance' >= 0 if `touse'
            if _rc {
                di as err "variance() must contain nonnegative variances"
                exit 459
            }
            quietly generate double `sigma_from_se' = `variance' / `sigma2' if `touse'
        }
        capture assert !missing(`sigma_from_se') if `touse'
        if _rc {
            di as err "scaled covariance overflowed; rescale the estimates or sigma2()"
            exit 459
        }
        local sigma `sigma_from_se'
    }
    local q : word count `sigma'
    if `d' == 1 & `q' == 1 local covtype 1
    else if `d' > 1 & `q' == `d' local covtype 2
    else if `d' > 1 & `q' == `d2' local covtype 3
    else {
        di as err "sigma() must contain d diagonal entries or d*d full covariance entries"
        exit 198
    }

    if "`generate'" != "" & "`prefix'" != "" {
        di as err "only one of generate() and prefix() may be specified"
        exit 198
    }
    if "`generate'" != "" & `d' > 1 {
        di as err "generate() requires scalar estimates; use prefix() for vectors"
        exit 198
    }
    if "`nomatrices'" != "" {
        if "`generate'`prefix'" == "" {
            di as err "nomatrices requires generate() or prefix() for the estimates"
            exit 198
        }
        if "`returnweights'" != "" {
            di as err "nomatrices may not be combined with returnweights"
            exit 198
        }
    }
    else if max(`Jall', `d') > c(max_matsize) {
        di as err "returned matrices exceed this Stata edition's matrix limit; use nomatrices"
        exit 908
    }

    // Reject impossible wide outputs before building a potentially huge list
    // of weight variable names.
    local nout = ("`generate'" != "") + `d' * ("`prefix'" != "") + ///
        `d' * ("`tildeprefix'" != "") + `Jall' * ("`weightsprefix'" != "")
    if c(k) + `nout' > c(maxvar) {
        di as err "insufficient variable capacity to stage the requested outputs"
        di as err "increase maxvar or request fewer generated variables"
        exit 900
    }

    local theta_vars `generate'
    if "`prefix'" != "" {
        forvalues k = 1/`d' {
            local theta_vars `theta_vars' `prefix'`k'
        }
    }
    local tilde_vars
    if "`tildeprefix'" != "" {
        forvalues k = 1/`d' {
            local tilde_vars `tilde_vars' `tildeprefix'`k'
        }
    }
    local weight_vars
    if "`weightsprefix'" != "" {
        forvalues k = 1/`Jall' {
            local weight_vars `weight_vars' `weightsprefix'`k'
        }
    }
    local all_output `theta_vars' `tilde_vars' `weight_vars'
    _ewgroup_validate_outputs, outputs("`all_output'") protected("`protected'") `replace'

    // Mata writes temporary doubles; originals remain untouched during computation.
    local theta_stage
    local tilde_stage
    local weight_stage
    foreach kind in theta tilde weight {
        foreach v of local `kind'_vars {
            tempvar stage
            quietly generate double `stage' = .
            local `kind'_stage ``kind'_stage' `stage'
        }
    }
    local all_stage `theta_stage' `tilde_stage' `weight_stage'

    capture mata: assert(_ewgroup_version() == 200)
    if _rc {
        findfile ewgroup_mata.mata
        quietly do "`r(fn)'"
    }
    tempname Theta Tilde Deriv Weights
    tempname alpha alpha_u gamma_used sure_A sure_D max_lambda
    local theta_mat
    local tilde_mat
    local deriv_mat
    local weights_mat
    if "`nomatrices'" == "" {
        local theta_mat `Theta'
        local tilde_mat `Tilde'
        local deriv_mat `Deriv'
    }
    if "`returnweights'" != "" local weights_mat `Weights'

    mata: _ewgroup_core_stata("`varlist'", "`sigma'", "`touse'", ///
        `sigma2', `gamma', `factor', `covtype', "`theta_stage'", ///
        "`tilde_stage'", "`weight_stage'", "`theta_mat'", "`tilde_mat'", ///
        "`weights_mat'", "`deriv_mat'", "`alpha'", "`alpha_u'", ///
        "`gamma_used'", "`sure_A'", "`sure_D'", "`max_lambda'")

    if "`nomatrices'" == "" {
        matrix colnames `Theta' = `varlist'
        matrix colnames `Tilde' = `varlist'
        matrix colnames `Deriv' = derivative_trace
    }
    local k 0
    foreach v of local all_output {
        local ++k
        local stage : word `k' of `all_stage'
        capture confirm variable `v', exact
        if !_rc {
            quietly recast double `v'
            quietly replace `v' = `stage' if `touse'
        }
        else {
            quietly rename `stage' `v'
        }
    }

    return scalar alpha = `alpha'
    return scalar alpha_unconstrained = `alpha_u'
    return scalar gamma = `gamma_used'
    return scalar sigma2 = `sigma2'
    return scalar factor = `factor'
    return scalar sure_A = `sure_A'
    return scalar sure_D = `sure_D'
    return scalar max_lambda = `max_lambda'
    return scalar J = `Jall'
    return scalar d = `d'
    return local coefnames "`varlist'"
    if "`nomatrices'" == "" {
        return matrix theta = `Theta'
        return matrix tilde = `Tilde'
        return matrix derivative_trace = `Deriv'
    }
    if "`returnweights'" != "" return matrix weights = `Weights'
end
