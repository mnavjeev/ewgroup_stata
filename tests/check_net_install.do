version 16.0
clear all
set more off
set type double

local pkgdir = c(pwd)
capture confirm file "`pkgdir'/ewgroup.pkg"
if (_rc) {
    local pkgdir "`pkgdir'/.."
}
capture confirm file "`pkgdir'/ewgroup.pkg"
if (_rc) {
    di as err "run this script from ewgroup_stata/ or ewgroup_stata/tests/"
    exit 601
}

local oldplus "`c(sysdir_plus)'"
sysdir set PLUS "`c(tmpdir)'"
capture noisily net uninstall ewgroup
net install ewgroup, from("`pkgdir'") replace

cd "`c(tmpdir)'"
clear
input double beta_hat se_hat
-1.10 .200
-0.95 .235
 0.20 .212
 0.27 .245
 1.30 .224
end

ewgroup_core beta_hat, se(se_hat) generate(theta) replace
assert !missing(theta)

sysdir set PLUS "`oldplus'"
di as result "local net install check passed"
