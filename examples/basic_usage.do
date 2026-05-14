version 16.0
clear all
set more off
set type double

capture confirm file "ewgroup.ado"
if (_rc == 0) {
    adopath ++ "."
}
else {
    capture confirm file "../ewgroup.ado"
    if (_rc == 0) {
        adopath ++ ".."
    }
}

di as text "Example 1: estimating grouped cell means"
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

di as text "Example 2: estimating cell-specific slopes and intercepts"
clear
set obs 80
generate byte w = ceil(_n / 20)
generate double x = mod(_n - 1, 20) / 10
generate double y = 1 + .2*w + (.1*w)*x + sin(_n)/20

ewgroup y x, group(w) prefix(theta_) tildeprefix(tilde_) replace
summarize theta_1 theta_2 tilde_1 tilde_2
matrix list e(beta_hat)
matrix list e(theta)

di as text "Example 3: using precomputed cell estimates"
clear
input double beta_hat se_hat
-1.10 .200
-0.95 .235
 0.20 .212
 0.27 .245
 1.30 .224
end

ewgroup_core beta_hat, se(se_hat) ///
    generate(theta) tildeprefix(tilde_) returnweights replace
list beta_hat se_hat theta tilde_1
matrix list r(weights)
