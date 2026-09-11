*! version 0.2.0 11sep2026
// Summary estimates plus uncertainty, or raw observations plus group().
program define ewgroup, eclass
    version 16.0
    // Keep both data and prior estimation results intact on errors or Break.
    tempname previous
    _estimates hold `previous', restore nullok copy
    preserve
    _ewgroup_work `0'
    nobreak {
        _estimates unhold `previous', not
        restore, not
    }
end

program define _ewgroup_work, eclass
    version 16.0
    syntax varlist(numeric min=1) [if] [in], ///
        [GROUP(varname) SIGma(varlist numeric min=1) SE(varname numeric) ///
        VARiance(varname numeric) NOCONstant SIGMA2(string) GAMMA(string) ///
        FACTOR(real 0.2) GENerate(name) PREFix(name) TILDEPrefix(name) ///
        WEIGHTSPrefix(name) RETURNWeights NOMATRices REPLACE]

    local covmodes = ("`sigma'" != "") + ("`se'" != "") + ("`variance'" != "")
    if ("`group'" != "" & `covmodes' != 0) | ("`group'" == "" & `covmodes' != 1) {
        di as err "specify group() for raw data, or exactly one of se(), variance(), sigma() for estimates"
        exit 198
    }
    if "`group'" == "" & "`noconstant'" != "" {
        di as err "noconstant is allowed only with group()"
        exit 198
    }
    if "`group'" != "" & "`weightsprefix'`returnweights'" != "" {
        di as err "weightsprefix() and returnweights are supported with summary estimates only"
        exit 198
    }
    foreach opt in sigma2 gamma {
        if "``opt''" != "" {
            capture confirm number ``opt''
            if _rc {
                di as err "`opt'() must be finite and positive"
                exit 198
            }
            if missing(``opt'') | ``opt'' <= 0 {
                di as err "`opt'() must be finite and positive"
                exit 198
            }
        }
    }
    if missing(`factor') | `factor' <= 0 {
        di as err "factor() must be finite and positive"
        exit 198
    }

    marksample touse
    tempname Theta Beta Tilde Cells Ncell Deriv Weights
    tempname alpha alpha_u gamma_used sigma2_used sure_A sure_D max_lambda Jcells
    local depvar
    local xvars

    if "`group'" == "" {
        local mode summary
        markout `touse' `sigma' `se' `variance'
        quietly count if `touse'
        local N = r(N)
        if `N' == 0 error 2000
        local d : word count `varlist'
        local coefnames `varlist'
        local options `replace' `returnweights' `nomatrices' factor(`factor')
        foreach opt in sigma se variance sigma2 gamma generate prefix tildeprefix weightsprefix {
            if "``opt''" != "" local options `options' `opt'(``opt'')
        }
        quietly ewgroup_core `varlist' if `touse', `options'
        scalar `alpha' = r(alpha)
        scalar `alpha_u' = r(alpha_unconstrained)
        scalar `gamma_used' = r(gamma)
        scalar `sigma2_used' = r(sigma2)
        scalar `sure_A' = r(sure_A)
        scalar `sure_D' = r(sure_D)
        scalar `max_lambda' = r(max_lambda)
        scalar `Jcells' = r(J)
        if "`nomatrices'" == "" {
            matrix `Theta' = r(theta)
            matrix `Tilde' = r(tilde)
            matrix `Deriv' = r(derivative_trace)
            if "`returnweights'" != "" matrix `Weights' = r(weights)
            mata: st_matrix("`Beta'", st_data(selectindex(st_data(., "`touse'")), tokens("`varlist'")))
            mata: st_matrix("`Cells'", selectindex(st_data(., "`touse'")))
            mata: st_matrix("`Ncell'", J(`N', 1, 1))
        }
    }
    else {
        local mode raw
        gettoken depvar xvars : varlist
        capture confirm numeric variable `group'
        if !_rc {
            markout `touse' `group'
            local group_use `group'
        }
        else {
            quietly replace `touse' = 0 if `group' == ""
            tempvar group_num
            quietly egen long `group_num' = group(`group') if `touse', label
            local group_use `group_num'
        }
        quietly count if `touse'
        local N = r(N)
        if `N' == 0 error 2000
        if "`xvars'" == "" & "`noconstant'" != "" {
            di as err "noconstant is not allowed when no covariates are specified"
            exit 198
        }
        local d : word count `xvars'
        local coefnames `xvars'
        if "`noconstant'" == "" {
            local ++d
            local coefnames `coefnames' _cons
        }
        if "`sigma2'" == "" local sigma2 -1
        if "`gamma'" == "" local gamma .
        if "`generate'" != "" & "`prefix'" != "" {
            di as err "only one of generate() and prefix() may be specified"
            exit 198
        }
        if "`generate'" != "" & `d' > 1 {
            di as err "generate() requires scalar estimates; use prefix() for vectors"
            exit 198
        }
        if "`nomatrices'" != "" & "`generate'`prefix'" == "" {
            di as err "nomatrices requires generate() or prefix() for the estimates"
            exit 198
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
        // Count groups once for a useful preflight matrix-size error.
        tempvar tag
        quietly egen byte `tag' = tag(`group_use') if `touse'
        quietly count if `tag' & `touse'
        scalar `Jcells' = r(N)
        if "`nomatrices'" == "" & max(`Jcells', `d') > c(max_matsize) {
            di as err "returned matrices exceed this Stata edition's matrix limit; use nomatrices"
            exit 908
        }
        local all_output `theta_vars' `tilde_vars'
        _ewgroup_validate_outputs, outputs("`all_output'") protected("`varlist' `group'") `replace'
        local theta_stage
        local tilde_stage
        foreach kind in theta tilde {
            foreach v of local `kind'_vars {
                tempvar stage
                quietly generate double `stage' = .
                local `kind'_stage ``kind'_stage' `stage'
            }
        }
        local all_stage `theta_stage' `tilde_stage'
        capture mata: assert(_ewgroup_version() == 200)
        if _rc {
            findfile ewgroup_mata.mata
            quietly do "`r(fn)'"
        }
        local theta_mat
        local beta_mat
        local tilde_mat
        local cells_mat
        local ncell_mat
        local deriv_mat
        if "`nomatrices'" == "" {
            local theta_mat `Theta'
            local beta_mat `Beta'
            local tilde_mat `Tilde'
            local cells_mat `Cells'
            local ncell_mat `Ncell'
            local deriv_mat `Deriv'
        }
        local nocons = ("`noconstant'" != "")
        mata: _ewgroup_raw_stata("`depvar'", "`xvars'", "`group_use'", ///
            "`touse'", `nocons', `sigma2', `gamma', `factor', ///
            "`theta_stage'", "`tilde_stage'", "`theta_mat'", "`beta_mat'", "`tilde_mat'", ///
            "`cells_mat'", "`ncell_mat'", "`deriv_mat'", "`alpha'", "`alpha_u'", ///
            "`gamma_used'", "`sigma2_used'", "`sure_A'", "`sure_D'", ///
            "`max_lambda'", "`Jcells'")

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
    }

    if "`nomatrices'" == "" {
        matrix colnames `Theta' = `coefnames'
        matrix colnames `Beta' = `coefnames'
        matrix colnames `Tilde' = `coefnames'
        matrix colnames `Deriv' = derivative_trace
        matrix colnames `Cells' = cell
        matrix colnames `Ncell' = N
    }
    ereturn post, esample(`touse') obs(`N')
    ereturn scalar J = `Jcells'
    ereturn scalar d = `d'
    ereturn scalar alpha = `alpha'
    ereturn scalar alpha_unconstrained = `alpha_u'
    ereturn scalar gamma = `gamma_used'
    ereturn scalar sigma2 = `sigma2_used'
    ereturn scalar factor = `factor'
    ereturn scalar sure_A = `sure_A'
    ereturn scalar sure_D = `sure_D'
    ereturn scalar max_lambda = `max_lambda'
    if "`nomatrices'" == "" {
        ereturn matrix theta = `Theta'
        ereturn matrix beta_hat = `Beta'
        ereturn matrix tilde = `Tilde'
        ereturn matrix cells = `Cells'
        ereturn matrix ncell = `Ncell'
        ereturn matrix derivative_trace = `Deriv'
        if "`returnweights'" != "" ereturn matrix weights = `Weights'
    }
    ereturn local cmd "ewgroup"
    ereturn local cmdline `"ewgroup `0'"'
    ereturn local mode "`mode'"
    ereturn local depvar "`depvar'"
    ereturn local xvars "`xvars'"
    ereturn local group "`group'"
    ereturn local coefnames "`coefnames'"
    ereturn local noconstant "`noconstant'"
    ereturn local nomatrices "`nomatrices'"

    di as text "Exponentially weighted grouped estimator (`mode' data)"
    di as text "  observations: " as result %9.0g e(N)
    di as text "  cells:        " as result %9.0g e(J)
    di as text "  dimension:    " as result %9.0g e(d)
    di as text "  gamma:        " as result %9.6g e(gamma)
    di as text "  alpha:        " as result %9.6g e(alpha)
end
