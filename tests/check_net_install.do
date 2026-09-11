version 16.0
// Run this do-file in a fresh batch session. Isolate installation from the
// user's PLUS directory and restore paths and data if a check fails.
set more off
local oldpwd "`c(pwd)'"
local pkgdir "`c(pwd)'"
capture confirm file "`pkgdir'/ewgroup.pkg"
if (_rc) local pkgdir "`pkgdir'/.."
confirm file "`pkgdir'/ewgroup.pkg"
local oldplus "`c(sysdir_plus)'"
local oldado : copy global S_ADO
tempfile location
local installroot "`location'_install"
mkdir "`installroot'"
mkdir "`installroot'/plus"
preserve
capture noisily {
    sysdir set PLUS "`installroot'/plus"
    // Exclude the checkout and any user-installed copy from lookup.
    global S_ADO "BASE;PLUS"
    net install ewgroup, from("`pkgdir'") replace
    cd "`installroot'"
    capture program drop ewgroup
    capture program drop ewgroup_core
    capture program drop _ewgroup_validate_outputs
    mata: mata clear
    which ewgroup
    which ewgroup_core
    which _ewgroup_validate_outputs
    findfile ewgroup_mata.mata
    clear
    input double(beta_hat se_hat)
    -1.10 .200
    -0.95 .235
     0.20 .212
     0.27 .245
     1.30 .224
    end
    ewgroup beta_hat, se(se_hat) generate(theta)
    assert !missing(theta)
    assert e(sample) == 1
    ewgroup_core beta_hat, se(se_hat) generate(legacy)
    assert abs(theta - legacy) < 1e-12
    ewgroup beta_hat, se(se_hat) generate(compact) nomatrices
    assert abs(theta - compact) < 1e-12
}
local rc = _rc
restore
sysdir set PLUS "`oldplus'"
global S_ADO `"`oldado'"'
cd "`oldpwd'"
if (`rc') exit `rc'
display as result "local net install check passed"
