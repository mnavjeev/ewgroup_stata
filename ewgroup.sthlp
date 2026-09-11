{smcl}
{* *! version 0.2.0 11sep2026}{...}
{vieweralsosee "ewgroup_core" "help ewgroup_core"}{...}
{title:Title}

{p2colset 5 18 20 2}{...}
{p2col:{hi:ewgroup} {hline 2}}Exponentially weighted grouped-heterogeneity estimator for summary or raw data{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 15 2}
{cmd:ewgroup} {it:beta_vars} {ifin},
{cmd:se(}{it:sevar}{cmd:)} | {cmd:variance(}{it:varvar}{cmd:)} |
{cmd:sigma(}{it:covariance_vars}{cmd:)}
[{cmd:sigma2(}{it:#}{cmd:)} {cmd:gamma(}{it:#}{cmd:)}
 {cmd:factor(}{it:#}{cmd:)} {cmd:generate(}{it:newvar}{cmd:)}
 {cmd:prefix(}{it:stub}{cmd:)} {cmd:tildeprefix(}{it:stub}{cmd:)}
 {cmd:weightsprefix(}{it:stub}{cmd:)} {cmd:returnweights}
 {cmd:nomatrices} {cmd:replace}]

{pstd}Or, for raw observations:{p_end}

{p 8 15 2}
{cmd:ewgroup} {it:depvar} [{it:xvars}] {ifin},
{cmd:group(}{it:varname}{cmd:)}
[{cmd:noconstant}
 {cmd:sigma2(}{it:#}{cmd:)}
 {cmd:gamma(}{it:#}{cmd:)}
 {cmd:factor(}{it:#}{cmd:)}
 {cmd:generate(}{it:newvar}{cmd:)}
 {cmd:prefix(}{it:stub}{cmd:)}
 {cmd:tildeprefix(}{it:stub}{cmd:)}
 {cmd:nomatrices} {cmd:replace}]

{title:Description}

{pstd}
For precomputed estimates, supply one row per cell and exactly one uncertainty
option: {cmd:se()} or {cmd:variance()} for scalar estimates, or {cmd:sigma()}
with {cmd:sigma2()} for scalar or vector estimates. Do not also supply
{cmd:group()}. This is the usual estimate-and-standard-error workflow familiar
from empirical Bayes commands. The legacy {help ewgroup_core} accepts the same
summary data and continues to return {cmd:r()} results.

{phang2}{cmd:. ewgroup beta_hat, se(se_hat) generate(theta)}{p_end}

{pstd}
{cmd:ewgroup} estimates cell-specific least-squares coefficients by values of
{cmd:group()}, then combines information across cells whose preliminary
estimates are close.  The command first computes the usual cell-by-cell
estimate, uses exponential weights to borrow strength across similar cells,
applies the debiasing correction from the paper, and finally chooses how much
to use the weighted estimate by a feasible SURE rule.

{pstd}
If no covariates are supplied, {cmd:ewgroup} estimates cell means.  If
covariates are supplied, a constant is included by default, following Stata's
regression convention.  Specify {cmd:noconstant} to estimate the model without
an intercept.  Every cell needs more usable observations than coefficients, and its design
matrix must have full column rank to estimate residual uncertainty.

{pstd}
The command is useful when there are many values of a discrete covariate and
the researcher expects some cells to have similar coefficients, but does not
want to impose a known grouping of cells.

{title:Empirical Workflow}

{pstd}
The raw-data command does the following:

{phang2}
1. For each value of {cmd:group()}, run a separate least-squares regression.

{phang2}
2. Estimate the scaled covariance matrix of each cell-specific coefficient
using the HC0 residual formula.

{phang2}
3. Compare the preliminary coefficient estimates across cells, with the
comparison adjusted for their covariance estimates.

{phang2}
4. Construct the exponentially weighted estimate and combine it with the
cell-by-cell estimate using SURE.

{pstd}
The final cell-level estimates are stored in {cmd:e(theta)}.  If {cmd:generate()}
or {cmd:prefix()} is specified, the estimates are also written back to the data,
repeated for all observations in the same cell.

{title:Typical Use}

{pstd}
For cell means, where the only variable that changes by cell is the mean of
{it:y}, run:

{phang2}
{cmd:. ewgroup y, group(w) generate(theta)}

{pstd}
For cell-specific regressions, include the regressors after the dependent
variable:

{phang2}
{cmd:. ewgroup y x1 x2, group(w) prefix(theta_)}

{pstd}
In the first command, {cmd:theta} is the adjusted cell mean.  In the second
command, {cmd:theta_1}, {cmd:theta_2}, and so on are the adjusted regression
coefficients.  These generated variables repeat the same cell-level estimate
for every observation in the same cell.

{title:Options}

{phang}
{cmd:group(}{it:varname}{cmd:)} specifies the discrete cell variable.  Numeric
and string variables are allowed.  String groups are internally encoded; the
stored {cmd:e(cells)} matrix then contains encoded cell ids.

{phang}
{cmd:noconstant} omits the intercept when covariates are supplied.

{phang}
{cmd:sigma2(}{it:#}{cmd:)} sets the scale parameter.  If omitted, the default
with {cmd:group()} is {it:J}/{it:N}, where {it:J} is the number of observed cells and {it:N} is the
estimation sample size. With {cmd:se()} or {cmd:variance()}, the default is 1;
with {cmd:sigma()}, this option is required. Covariances are interpreted as
{cmd:sigma2()} times the supplied scaled covariance matrices.

{phang}
{cmd:gamma(}{it:#}{cmd:)} sets the exponential-weight tuning parameter.  If
omitted, the default is {cmd:factor()}/({it:d}*{cmd:max_lambda}), where {it:d}
is the number of coefficients and {cmd:max_lambda} is the largest eigenvalue
among the scaled covariance estimates.  Thus the default is proportional to
1/{it:d} when the covariance eigenvalues remain bounded.  Larger values make
the weights place more emphasis on cells with similar preliminary estimates.
The command checks that the value is small enough for the required matrix
inverses to exist.

{phang}
{cmd:factor(}{it:#}{cmd:)} sets the default-gamma multiplier.  The default is
{cmd:0.2}.

{phang}
{cmd:generate(}{it:newvar}{cmd:)} stores the final estimate in {it:newvar}.
This option is allowed only in the scalar case.

{phang}
{cmd:prefix(}{it:stub}{cmd:)} stores final estimates in variables
{it:stub}{cmd:1}, ..., {it:stub}{it:d}.  Values are repeated for observations
belonging to the same cell.

{pmore}
The order follows the regression coefficient order.  With one regressor and
the default intercept, {it:stub}{cmd:1} is the slope and {it:stub}{cmd:2} is the
intercept.

{phang}
{cmd:tildeprefix(}{it:stub}{cmd:)} stores the exponentially weighted estimates
before SURE recombination in variables {it:stub}{cmd:1}, ...,
{it:stub}{it:d}.

{phang}
{cmd:se(}{it:varname}{cmd:)} supplies standard errors for one scalar estimate
per cell. {cmd:variance()} instead supplies variances. Both must be nonnegative.

{phang}
{cmd:sigma(}{it:varlist}{cmd:)} supplies scaled covariances: one variable for
scalar estimates, {it:d} variables for diagonal vector covariances, or {it:d*d}
variables in row order for full covariances (for example, s11 s12 s21 s22).
Matrices must be symmetric and positive semidefinite; only relative roundoff
errors are corrected. See {help ewgroup_core} for detailed covariance examples.

{phang}
{cmd:weightsprefix(}{it:stub}{cmd:)} generates one weight variable per cell.
{cmd:returnweights} returns {cmd:e(weights)}. These options are available only
for summary inputs and explicitly request quadratic storage.

{phang}
{cmd:nomatrices} omits returned matrices while retaining scalar results and
{cmd:e(sample)}. It requires {cmd:generate()} or {cmd:prefix()} and cannot be
combined with {cmd:returnweights}. Use it for large jobs. Without it, result
matrices must fit within the current Stata edition's matrix-operation limits.

{phang}
{cmd:replace} permits existing numeric outputs and promotes them to double.
Values outside the estimation sample remain unchanged; new outputs are missing
there. Names must be unique and cannot overlap any input. Errors and Break
leave the data and previous estimation results unchanged.

{pstd}
The SURE mixing weight is invariant to a common change of units. When the
weighted correction is exactly zero, {cmd:alpha} is set to zero. If every
summary covariance is zero, specify {cmd:gamma()} explicitly. Custom gamma
must satisfy gamma times the largest scaled covariance eigenvalue less than
one; the default uses the draft's 0.2/(dimension times maximum eigenvalue).

{title:Stored results}

{pstd}
{cmd:ewgroup} stores the following in {cmd:e()}. Matrices are omitted with
{cmd:nomatrices}; {cmd:e(sample)} marks the observations used in either mode:

{synoptset 24 tabbed}{...}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(J)}}number of cells{p_end}
{synopt:{cmd:e(d)}}coefficient dimension{p_end}
{synopt:{cmd:e(alpha)}}SURE mixing weight{p_end}
{synopt:{cmd:e(alpha_unconstrained)}}unconstrained SURE mixing weight{p_end}
{synopt:{cmd:e(mode)}}summary or raw{p_end}
{synopt:{cmd:e(factor)}}default-gamma multiplier{p_end}
{synopt:{cmd:e(weights)}}weights when returnweights is requested{p_end}
{synopt:{cmd:e(gamma)}}tuning parameter used{p_end}
{synopt:{cmd:e(sigma2)}}scale parameter used{p_end}
{synopt:{cmd:e(sure_A)}}quadratic SURE coefficient{p_end}
{synopt:{cmd:e(sure_D)}}linear SURE coefficient divided by two{p_end}
{synopt:{cmd:e(max_lambda)}}largest covariance eigenvalue{p_end}
{synopt:{cmd:e(theta)}}final cell-level estimates{p_end}
{synopt:{cmd:e(beta_hat)}}naive cell-level estimates{p_end}
{synopt:{cmd:e(tilde)}}exponentially weighted estimates{p_end}
{synopt:{cmd:e(cells)}}raw group ids, or original observation numbers for summary data{p_end}
{synopt:{cmd:e(ncell)}}raw cell sample sizes, or ones for summary data{p_end}
{synopt:{cmd:e(derivative_trace)}}cell-level divergence terms used by SURE{p_end}
{synopt:{cmd:e(cmd)}}{cmd:ewgroup}{p_end}

{title:Examples}

{pstd}Cell means:{p_end}
{phang2}{cmd:. ewgroup y, group(w) generate(theta) replace}{p_end}

{pstd}Cell-specific slopes and intercepts:{p_end}
{phang2}{cmd:. ewgroup y x1 x2, group(w) prefix(theta_) tildeprefix(tilde_) replace}{p_end}

{pstd}Without intercept:{p_end}
{phang2}{cmd:. ewgroup y x1 x2, group(w) noconstant prefix(theta_) replace}{p_end}

{pstd}Inspecting the cell-level estimates:{p_end}
{phang2}{cmd:. matrix list e(beta_hat)}{p_end}
{phang2}{cmd:. matrix list e(theta)}{p_end}

{pstd}Running the bundled example file from the package directory:{p_end}
{phang2}{cmd:. do examples/basic_usage.do}{p_end}

{title:Low-level interface}

{pstd}
Use {help ewgroup_core} when existing code expects {cmd:r()} results for
preliminary estimates and their uncertainties. It computes the same estimator:

{phang2}{cmd:. ewgroup_core beta_hat, se(se_hat) generate(theta)}{p_end}

{title:Remarks}

{pstd}
The exact calculation compares every pair of cells.  This is usually fine for
large jobs with {cmd:nomatrices}, but runtime remains quadratic in cell count.
Working memory is linear at fixed dimension unless full weights are requested.

{pstd}
Rows of {cmd:e(theta)}, {cmd:e(beta_hat)}, and {cmd:e(tilde)} follow original
observation order for summary data and sorted cell ids for raw data. For string
groups, Stata internally assigns numeric ids in sorted string order.

{pstd}
Extremely small {cmd:sure_A} or {cmd:sure_D} diagnostics can underflow to zero;
the mixing ratio is calculated with scaled arithmetic. An unconstrained alpha
outside Stata's numeric range is returned as missing, while the constrained
{cmd:alpha} remains zero or one.

{title:Author}

{pstd}
Manu Navjeevan, Denis Chetverikov, and Andrei Voronin.
