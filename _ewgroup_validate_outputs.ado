*! version 0.2.0 11sep2026
// Validate every output before changing any original variable.
program define _ewgroup_validate_outputs
    version 16.0
    syntax, [OUTPUTS(string) PROTECTED(string) REPLACE]
    local seen
    foreach v of local outputs {
        confirm name `v'
        if substr("`v'", 1, 2) == "__" {
            di as err "output variable `v' uses Stata's reserved temporary-name prefix"
            exit 198
        }
        if strpos(" `seen' ", " `v' ") {
            di as err "output variable `v' is requested more than once"
            exit 198
        }
        local seen `seen' `v'
        if strpos(" `protected' ", " `v' ") {
            di as err "output variable `v' conflicts with an input variable"
            exit 110
        }
        capture confirm variable `v', exact
        if !_rc {
            if "`replace'" == "" {
                di as err "`v' already exists; specify replace to overwrite"
                exit 110
            }
            confirm numeric variable `v', exact
        }
        else {
            confirm new variable `v'
        }
    }
    // One temporary double variable is needed for every requested output.
    local nout : word count `outputs'
    if c(k) + `nout' > c(maxvar) {
        di as err "insufficient variable capacity to stage the requested outputs"
        di as err "increase maxvar or request fewer generated variables"
        exit 900
    }
end
