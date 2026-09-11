{smcl}
{* *! version 0.2.0 11sep2026}{...}
{vieweralsosee "ewgroup" "help ewgroup"}{...}
{title:Title}

{p2colset 5 22 24 2}{...}
{p2col:{hi:ewgroup_core} {hline 2}}Exponentially weighted estimator from cell-level summaries{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 15 2}
{cmd:ewgroup_core} {it:beta_vars} {ifin},
[{cmd:sigma(}{it:sigma_vars}{cmd:)}
 {cmd:se(}{it:varname}{cmd:)}
 {cmd:variance(}{it:varname}{cmd:)}
 {cmd:sigma2(}{it:#}{cmd:)}
 {cmd:gamma(}{it:#}{cmd:)}
 {cmd:factor(}{it:#}{cmd:)}
 {cmd:generate(}{it:newvar}{cmd:)}
 {cmd:prefix(}{it:stub}{cmd:)}
 {cmd:tildeprefix(}{it:stub}{cmd:)}
 {cmd:weightsprefix(}{it:stub}{cmd:)}
 {cmd:returnweights}
 {cmd:nomatrices} {cmd:replace}]

{title:Description}

{pstd}
This is the compatible {cmd:rclass} summary interface. New workflows may use
{cmd:ewgroup beta_hat, se(se_hat) generate(theta)} for {cmd:e()} results and
{cmd:e(sample)}. Both commands use the same numerical estimator.

{pstd}
{cmd:ewgroup_core} is the cell-summary interface for the exponentially weighted
grouped-heterogeneity estimator.  Use it when the first-stage cell estimates
have already been computed outside this command.  The data must contain one
observation per cell.  The variables in {it:beta_vars} contain the preliminary
cell estimates.

{pstd}
For scalar estimates, the most familiar workflow is to pass a point estimate
and a standard error, as in many empirical Bayes commands:

{phang2}
{cmd:. ewgroup_core beta_hat, se(se_hat) generate(theta)}

{pstd}
If you have variances rather than standard errors, use {cmd:variance()}.
For vector-valued estimates, use {cmd:sigma()} to supply diagonal or full
covariance estimates.

{title:Typical Use}

{pstd}
The common scalar workflow is:

{phang2}
1. Create a dataset with one row per unit.

{phang2}
2. Put the original estimate in one variable, for example {cmd:beta_hat}.

{phang2}
3. Put its standard error in another variable, for example {cmd:se_hat}.

{phang2}
4. Run:

{phang2}
{cmd:. ewgroup_core beta_hat, se(se_hat) generate(theta)}

{pstd}
The new variable {cmd:theta} is the adjusted estimate.  The original estimate
is not changed.

{pstd}
When {cmd:sigma()} is used, the covariance of the preliminary estimate in cell
{it:j} is interpreted as {cmd:sigma2()} times {it:Sigma_hat_j}.  This matches
the notation used in the paper and in the companion R package.

{title:Input Formats}

{pstd}
The command accepts three covariance formats:

{phang2}
1. Scalar standard error: one coefficient variable and {cmd:se(}{it:sevar}{cmd:)}.

{phang2}
2. Scalar variance: one coefficient variable and
{cmd:variance(}{it:varvar}{cmd:)}.

{phang2}
3. Scaled covariance: {cmd:sigma()} plus {cmd:sigma2()}.  For scalar estimates,
{cmd:sigma()} contains one variable.

{phang2}
4. Diagonal vector case: {it:d} coefficient variables and {it:d} covariance
variables.  The covariance variables are the diagonal entries.

{phang2}
5. Full vector case: {it:d} coefficient variables and {it:d*d} covariance
variables.  The covariance variables are read in row-major order.

{pstd}
For example, if {it:d}=2 and the covariance matrix is

{pmore}
({it:S11}, {it:S12}; {it:S21}, {it:S22}),

{pstd}
then specify {cmd:sigma(S11 S12 S21 S22)}.

{title:Options}

{phang}
{cmd:se(}{it:varname}{cmd:)} specifies standard errors for scalar preliminary
estimates.  If {cmd:sigma2()} is omitted, the command uses {cmd:sigma2(1)}.

{phang}
{cmd:variance(}{it:varname}{cmd:)} specifies variances for scalar preliminary
estimates.  If {cmd:sigma2()} is omitted, the command uses {cmd:sigma2(1)}.

{phang}
{cmd:sigma(}{it:sigma_vars}{cmd:)} specifies scaled covariance estimates.
When {cmd:sigma()} is used, {cmd:sigma2()} must also be specified.

{phang}
{cmd:sigma2(}{it:#}{cmd:)} specifies the positive scale parameter.  With
{cmd:se()} or {cmd:variance()}, this option is optional.  With {cmd:sigma()},
it is required.

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
{it:stub}{cmd:1}, ..., {it:stub}{it:d}.

{pmore}
The variables follow the same order as {it:beta_vars}.

{phang}
{cmd:tildeprefix(}{it:stub}{cmd:)} stores the exponentially weighted estimates
before SURE recombination in variables {it:stub}{cmd:1}, ...,
{it:stub}{it:d}.

{phang}
{cmd:weightsprefix(}{it:stub}{cmd:)} stores the {it:J} by {it:J} matrix of
exponential weights as variables {it:stub}{cmd:1}, ..., {it:stub}{it:J}.  Row
{it:j}, column {it:k} is the weight placed by cell {it:j} on cell {it:k}.

{phang}
{cmd:returnweights} stores the weight matrix in {cmd:r(weights)}.

{phang}
{cmd:nomatrices} omits all returned matrices while keeping scalar diagnostics.
It requires {cmd:generate()} or {cmd:prefix()} and cannot be combined with
{cmd:returnweights}. Use this option for large jobs. Without it, matrices must
fit within the current Stata edition's matrix-operation limits.

{phang}
{cmd:replace} permits existing numeric outputs and promotes them to double.
Values outside {cmd:if/in} and rows excluded for missing inputs are preserved.
New outputs are missing outside that sample. Output names must be distinct and
cannot overlap estimate or uncertainty inputs. Errors and Break preserve data.

{pstd}
Tuning options must be finite and positive. Covariances must be symmetric and
positive semidefinite; only relative roundoff errors are corrected. If all
covariances are zero, choose {cmd:gamma()} explicitly.

{title:Stored results}

{pstd}
{cmd:ewgroup_core} stores the following in {cmd:r()}; matrices are omitted with
{cmd:nomatrices}:

{synoptset 28 tabbed}{...}
{synopt:{cmd:r(alpha)}}SURE mixing weight{p_end}
{synopt:{cmd:r(alpha_unconstrained)}}unconstrained SURE mixing weight{p_end}
{synopt:{cmd:r(gamma)}}tuning parameter used{p_end}
{synopt:{cmd:r(sigma2)}}scale parameter used{p_end}
{synopt:{cmd:r(sure_A)}}quadratic SURE coefficient{p_end}
{synopt:{cmd:r(sure_D)}}linear SURE coefficient divided by two{p_end}
{synopt:{cmd:r(max_lambda)}}largest covariance eigenvalue{p_end}
{synopt:{cmd:r(J)}}number of cells{p_end}
{synopt:{cmd:r(d)}}coefficient dimension{p_end}
{synopt:{cmd:r(theta)}}final estimates{p_end}
{synopt:{cmd:r(tilde)}}exponentially weighted estimates{p_end}
{synopt:{cmd:r(derivative_trace)}}cell-level divergence terms used by SURE{p_end}
{synopt:{cmd:r(weights)}}exponential weights, if {cmd:returnweights} is specified{p_end}

{title:Examples}

{pstd}Scalar estimates with standard errors:{p_end}
{phang2}{cmd:. ewgroup_core bhat, se(sehat) generate(theta)}{p_end}

{pstd}Scalar estimates with variances:{p_end}
{phang2}{cmd:. ewgroup_core bhat, variance(vhat) generate(theta)}{p_end}

{pstd}Scalar estimates with scaled covariance notation:{p_end}
{phang2}{cmd:. ewgroup_core bhat, sigma(Sigmahat) sigma2(.05) gamma(.2) generate(theta)}{p_end}

{pstd}Vector estimates with diagonal covariance variables:{p_end}
{phang2}{cmd:. ewgroup_core b1 b2, sigma(s11 s22) sigma2(.1) gamma(.05) prefix(theta_)}{p_end}

{pstd}Vector estimates with full covariance variables in row-major order:{p_end}
{phang2}{cmd:. ewgroup_core b1 b2, sigma(S11 S12 S21 S22) sigma2(.1) gamma(.05) prefix(theta_)}{p_end}

{pstd}Return the full matrix of exponential weights:{p_end}
{phang2}{cmd:. ewgroup_core bhat, sigma(vhat) sigma2(.05) returnweights}{p_end}
{phang2}{cmd:. matrix list r(weights)}{p_end}

{title:Remarks}

{pstd}
Rows are treated as cells.  If your data contain observation-level records, use
{help ewgroup} instead; it will compute the first-stage cell estimates before
calling the same Mata estimator.

{pstd}
The exact calculation compares every pair of cells.  This is usually fine for
tens of thousands of cells with {cmd:nomatrices}, but runtime is quadratic.
Working memory is linear in cell count at fixed dimension unless full weights
are requested. Full weights require quadratic memory.

{pstd}
Extremely small SURE diagnostics can underflow to zero without changing the
scaled mixing ratio. An unconstrained alpha beyond Stata's numeric range is
returned as missing; the constrained alpha remains zero or one.

{title:Author}

{pstd}
Manu Navjeevan, Denis Chetverikov, and Andrei Voronin.
