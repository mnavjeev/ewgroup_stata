{smcl}
{* *! version 0.1.0 14may2026}{...}
{vieweralsosee "ewgroup_core" "help ewgroup_core"}{...}
{title:Title}

{p2colset 5 18 20 2}{...}
{p2col:{hi:ewgroup} {hline 2}}Exponentially weighted grouped-heterogeneity estimator from raw data{p_end}
{p2colreset}{...}

{title:Syntax}

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
 {cmd:replace}]

{title:Description}

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
an intercept.  The command requires the within-cell design matrix to have full
column rank in every cell.

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
is {it:J}/{it:N}, where {it:J} is the number of observed cells and {it:N} is the
estimation sample size.  This is the normalization used in the paper.

{phang}
{cmd:gamma(}{it:#}{cmd:)} sets the exponential-weight tuning parameter.  If
omitted, the default is {cmd:factor()}/{cmd:max_lambda}, where
{cmd:max_lambda} is the largest eigenvalue among the scaled covariance
estimates.  Larger values make the weights place more emphasis on cells with
similar preliminary estimates.  The command checks that the value is small
enough for the required matrix inverses to exist.

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
{cmd:replace} allows generated output variables to already exist.

{title:Stored results}

{pstd}
{cmd:ewgroup} stores the following in {cmd:e()}:

{synoptset 24 tabbed}{...}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(J)}}number of cells{p_end}
{synopt:{cmd:e(d)}}coefficient dimension{p_end}
{synopt:{cmd:e(alpha)}}SURE mixing weight{p_end}
{synopt:{cmd:e(alpha_unconstrained)}}unconstrained SURE mixing weight{p_end}
{synopt:{cmd:e(gamma)}}tuning parameter used{p_end}
{synopt:{cmd:e(sigma2)}}scale parameter used{p_end}
{synopt:{cmd:e(sure_A)}}quadratic SURE coefficient{p_end}
{synopt:{cmd:e(sure_D)}}linear SURE coefficient divided by two{p_end}
{synopt:{cmd:e(max_lambda)}}largest covariance eigenvalue{p_end}
{synopt:{cmd:e(theta)}}final cell-level estimates{p_end}
{synopt:{cmd:e(beta_hat)}}naive cell-level estimates{p_end}
{synopt:{cmd:e(tilde)}}exponentially weighted estimates{p_end}
{synopt:{cmd:e(cells)}}cell ids corresponding to rows of the matrices{p_end}
{synopt:{cmd:e(ncell)}}cell sample sizes{p_end}
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
If you already have one row per cell with preliminary estimates and scaled
covariance estimates, use {help ewgroup_core}.  For scalar estimates,
{cmd:ewgroup_core} can also take the familiar empirical-Bayes-style input of a
point estimate and a standard error:

{phang2}{cmd:. ewgroup_core beta_hat, se(se_hat) generate(theta)}{p_end}

{title:Remarks}

{pstd}
The exact calculation compares every pair of cells.  This is usually fine for
hundreds or a few thousand cells.  Very large numbers of cells can require more
time and memory because the command must compute many pairwise comparisons.

{pstd}
Rows of {cmd:e(theta)}, {cmd:e(beta_hat)}, and {cmd:e(tilde)} are ordered by the
numeric cell id in {cmd:e(cells)}.  For string groups, Stata internally assigns
numeric ids before estimation.

{title:Author}

{pstd}
Manu Navjeevan, Denis Chetverikov, and Andrei Voronin.
