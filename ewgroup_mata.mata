version 16.0

// Drop old copies of these Mata functions if this file is run more than once in
// the same Stata session. `capture` prevents an error when a function was not
// already loaded.
capture mata: mata drop _ewgroup_version()
capture mata: mata drop _ewgroup_core_compute()
capture mata: mata drop _ewgroup_core_stata()
capture mata: mata drop _ewgroup_raw_stata()
capture mata: mata drop _ewgroup_unflatten()
capture mata: mata drop _ewgroup_flatten()
capture mata: mata drop _ewgroup_trace()
capture mata: mata drop _ewgroup_sym()
capture mata: mata drop _ewgroup_check_inverse()
capture mata: mata drop _ewgroup_max_lambda()
capture mata: mata drop _ewgroup_make_sigma_flat()
capture mata: mata drop _ewgroup_store_scalar()
capture mata: mata drop ewgroup_fit()

mata:

// Versioned loader sentinel prevents an installed update from reusing older
// Mata adapters that are still resident in the current Stata session.
real scalar _ewgroup_version()
{
    return(200)
}

// A struct is a named bundle of results. This one keeps all values produced by
// the estimator together so helper functions can pass one object around instead
// of many separate matrices and scalars.
struct ewgroup_fit {
    real matrix theta
    real matrix tilde
    real matrix weights
    real colvector derivative_trace
    real colvector trace_sigma
    real scalar alpha
    real scalar alpha_unconstrained
    real scalar gamma
    real scalar sigma2
    real scalar sure_A
    real scalar sure_D
    real scalar max_lambda
}

// Return a symmetric version of A by averaging it with its transpose. This
// removes tiny roundoff differences between mirrored entries.
real matrix _ewgroup_sym(real matrix A)
{
    return(A / 2 + A' / 2)
}

// Add up the diagonal entries of a matrix. Mata has matrix commands, but this
// small loop makes the intended calculation explicit.
real scalar _ewgroup_trace(real matrix A)
{
    real scalar i, out

    out = 0
    for (i = 1; i <= min((rows(A), cols(A))); i++) {
        out = out + A[i, i]
    }
    return(out)
}

// Turn one row of d*d numbers into a d by d matrix. The Stata wrappers store
// covariance matrices in rows before sending them to the core routine.
real matrix _ewgroup_unflatten(real rowvector x, real scalar d)
{
    real matrix A
    real scalar r, c, p

    A = J(d, d, 0)
    p = 1
    for (r = 1; r <= d; r++) {
        for (c = 1; c <= d; c++) {
            A[r, c] = x[p]
            p = p + 1
        }
    }
    return(A)
}

// Turn a d by d matrix into one row of d*d numbers. This is the reverse of
// _ewgroup_unflatten().
real rowvector _ewgroup_flatten(real matrix A)
{
    real rowvector x
    real scalar r, c, p, d

    d = rows(A)
    x = J(1, d * d, 0)
    p = 1
    for (r = 1; r <= d; r++) {
        for (c = 1; c <= d; c++) {
            x[p] = A[r, c]
            p = p + 1
        }
    }
    return(x)
}

// Check that a matrix can be inverted, then return its inverse. Inverting means
// finding a matrix that undoes the original matrix. The correction matrices
// must be positive definite; cholinv() reports failure instead of silently
// supplying a generalized inverse for a singular matrix.
real matrix _ewgroup_check_inverse(real matrix A, string scalar message)
{
    real matrix inverse

    inverse = cholinv(_ewgroup_sym(A))
    if (hasmissing(inverse)) {
        errprintf("%s\n", message)
        exit(506)
    }
    return(inverse)
}

// Compute the largest eigenvalue of a covariance matrix. An eigenvalue is a
// number that describes the size of the matrix in some direction; the largest
// one is used to keep gamma in a safe range.
real scalar _ewgroup_max_lambda(real matrix S)
{
    real colvector evals

    evals = symeigenvalues(_ewgroup_sym(S))
    return(max(evals))
}

// Store a Mata scalar in a Stata scalar, but only when a name was provided. This
// lets callers skip outputs they do not need.
void _ewgroup_store_scalar(string scalar name, real scalar value)
{
    if (strlen(name) > 0) {
        st_numscalar(name, value)
    }
}

// Convert the user's covariance input into one common row format. Each output
// row contains the d by d covariance matrix for one cell, flattened into d*d
// columns.
real matrix _ewgroup_make_sigma_flat(real matrix Sin, real scalar d,
                                     real scalar covtype)
{
    real scalar r
    real matrix Sflat

    if (d < 1 | d != floor(d) | !anyof((1, 2, 3), covtype)) {
        errprintf("invalid covariance dimensions or covariance type.\n")
        exit(503)
    }
    if ((covtype == 1 & (d != 1 | cols(Sin) != 1)) |
        (covtype == 2 & cols(Sin) != d) |
        (covtype == 3 & cols(Sin) != d*d)) {
        errprintf("covariance dimensions do not match beta estimates.\n")
        exit(503)
    }
    // Preserve full entries: validation must see any asymmetry before
    // roundoff-level symmetrization is allowed.
    if (covtype != 2) return(Sin)
    Sflat = J(rows(Sin), d*d, 0)
    for (r = 1; r <= d; r++) {
        Sflat[., (r-1)*d+r] = Sin[., r]
    }
    return(Sflat)
}

// Core estimator used by both Stata commands. B is the table of preliminary
// estimates, one row per cell. Sflat is the matching covariance information in
// flattened form. The function returns theta, tilde, alpha, weights, and
// diagnostic numbers in an ewgroup_fit struct.
struct ewgroup_fit scalar _ewgroup_core_compute(real matrix B,
                                                real matrix Sflat,
                                                real scalar covtype,
                                                real scalar sigma2,
                                                real scalar gamma,
                                                real scalar factor,
                                                real scalar want_weights)
{
    struct ewgroup_fit scalar fit
    real scalar Jn, d, j, k, r, c, is_diagonal, scale_S, max_lambda, tol
    real scalar rowsum_w, trace_update, adjustment_scale, norm_scaled
    real scalar trace_scale, trace_sum_scaled, log_A, log_abs_D, log_ratio
    real matrix Id, S, S_work, eigenvectors, Omega, H, Sigma_H, Q, contraction
    real matrix Cflat, Omega_flat, Sdiag, Cdiag, Odiag
    real matrix diff, scores, corrected_diff, centered_scores, adjustment
    real rowvector evals, N, s, score_bar, h, sigma_h
    real colvector logw, w, trace_terms, trace_sigma

    Jn = rows(B)
    d = cols(B)
    if (Jn < 1 | d < 1) {
        errprintf("beta estimates must have at least one row and one column.\n")
        exit(198)
    }
    if (rows(Sflat) != Jn | cols(Sflat) != d*d |
        !anyof((1, 2, 3), covtype) | (covtype == 1 & d != 1)) {
        errprintf("covariance dimensions do not match beta estimates.\n")
        exit(503)
    }
    if (hasmissing(B) | hasmissing(Sflat)) {
        errprintf("beta estimates and covariances must be finite.\n")
        exit(198)
    }
    if (sigma2 <= 0 | missing(sigma2) | factor <= 0 | missing(factor)) {
        errprintf("sigma2 and factor must be positive finite scalars.\n")
        exit(198)
    }
    // Check byte-count arithmetic before an optional quadratic allocation.
    // Allocation failures below this limit are reported by Mata itself.
    if (want_weights & Jn*Jn > (2^53-1)/8) {
        errprintf("requested weights matrix is too large; omit weight outputs.\n")
        exit(3900)
    }

    Id = I(d)
    S_work = Sflat
    Sdiag = J(Jn, d, 0)
    trace_sigma = J(Jn, 1, 0)
    max_lambda = 0
    is_diagonal = 1
    tol = sqrt(epsilon(1))
    // Validate covariance matrices in their own units, so a change of units
    // cannot hide asymmetry or negative eigenvalues behind an absolute cutoff.
    for (k = 1; k <= Jn; k++) {
        S = _ewgroup_unflatten(S_work[k, .], d)
        scale_S = max(abs(S))
        if (scale_S > 0) {
            if (max(abs(S/scale_S - S'/scale_S)) > tol) {
                errprintf("Sigma_hat must be symmetric (cell %g).\n", k)
                exit(506)
            }
            S = _ewgroup_sym(S)
        }
        if (sum(abs(S - diag(diagonal(S)))) == 0) {
            if (min(diagonal(S)) < 0) {
                errprintf("Sigma_hat has a negative variance (cell %g).\n", k)
                exit(506)
            }
            max_lambda = max((max_lambda, max(diagonal(S))))
        }
        else {
            is_diagonal = 0
            symeigensystem(S/scale_S, eigenvectors, evals)
            if (hasmissing(evals) | min(evals) < -tol) {
                errprintf("Sigma_hat must be positive semidefinite (cell %g).\n", k)
                exit(506)
            }
            // Remove only negative eigenvalues within roundoff tolerance.
            if (min(evals) < 0) {
                evals = evals :* (evals :> 0)
                S = _ewgroup_sym(eigenvectors * diag(evals) * eigenvectors') * scale_S
            }
            max_lambda = max((max_lambda, max(evals)*scale_S))
        }
        S_work[k, .] = _ewgroup_flatten(S)
        Sdiag[k, .] = diagonal(S)'
        trace_sigma[k] = sum(Sdiag[k, .])
    }
    if (hasmissing(trace_sigma) | missing(max_lambda)) {
        errprintf("covariance magnitudes exceed finite numerical range.\n")
        exit(430)
    }
    if (missing(gamma)) {
        if (max_lambda <= 0) {
            errprintf("Sigma_hat must have a positive eigenvalue to compute gamma.\n")
            exit(198)
        }
        gamma = (factor/d)/max_lambda
    }
    if (gamma <= 0 | missing(gamma)) {
        errprintf("gamma must be a positive finite scalar.\n")
        exit(198)
    }
    if (gamma*max_lambda >= 1) {
        errprintf("gamma must satisfy gamma * max_j lambda_max(Sigma_hat_j) < 1.\n")
        exit(198)
    }

    // Cache all candidate covariance transforms once. All paths process one
    // target row at a time; full weights are the only quadratic allocation.
    if (is_diagonal) {
        Cdiag = gamma * Sdiag
        Odiag = gamma :/ (1 :- Cdiag)
        if (hasmissing(Odiag)) {
            errprintf("covariance transform exceeds finite numerical range.\n")
            exit(430)
        }
    }
    else {
        Cflat = gamma * S_work
        Omega_flat = J(Jn, d*d, 0)
        for (k = 1; k <= Jn; k++) {
            S = _ewgroup_unflatten(Cflat[k, .], d)
            Omega = gamma * _ewgroup_check_inverse(Id - S,
                "encountered non-invertible distance matrix.")
            if (hasmissing(Omega)) {
                errprintf("covariance transform exceeds finite numerical range.\n")
                exit(430)
            }
            Omega_flat[k, .] = _ewgroup_flatten(Omega)
        }
    }
    adjustment = J(Jn, d, 0)
    trace_terms = J(Jn, 1, 0)
    if (want_weights) fit.weights = J(Jn, Jn, .)
    else fit.weights = J(0, 0, .)

    // Mata's normal break-key processing remains active in these loops.
    for (j = 1; j <= Jn; j++) {
        diff = B :- B[j, .]
        if (is_diagonal) scores = (diff :* Odiag) / sigma2
        else {
            scores = J(Jn, d, 0)
            for (r = 1; r <= d; r++) {
                for (c = 1; c <= d; c++) {
                    scores[., r] = scores[., r] +
                        Omega_flat[., (r-1)*d+c] :* diff[., c]
                }
            }
            scores = scores / sigma2
        }
        contraction = diff :* scores
        // Mata reductions omit missing entries; reject overflow before reducing.
        if (hasmissing(contraction) | hasmissing(scores)) {
            errprintf("distances exceed finite numerical range; rescale inputs.\n")
            exit(430)
        }
        logw = -0.5 * rowsum(contraction)
        if (hasmissing(logw)) {
            errprintf("distances exceed finite numerical range; rescale inputs.\n")
            exit(430)
        }
        w = exp(logw :- max(logw))
        rowsum_w = sum(w)
        if (rowsum_w <= 0 | missing(rowsum_w)) {
            errprintf("encountered invalid exponential weights.\n")
            exit(430)
        }
        w = w / rowsum_w
        if (want_weights) fit.weights[j, .] = w'

        // Accumulate differences, avoiding cancellation between two large means.
        N = colsum(w :* diff)
        score_bar = colsum(w :* scores)
        if (is_diagonal) {
            h = 1 :/ (1 :- colsum(w :* Cdiag))
            s = N :* h
            sigma_h = Sdiag[j, .] :* h
            corrected_diff = diff + Cdiag :* s
            if (hasmissing(corrected_diff) | hasmissing(sigma_h)) {
                errprintf("smoothing exceeds finite numerical range; rescale inputs.\n")
                exit(430)
            }
            centered_scores = scores :- score_bar
            contraction = (w :* centered_scores) :* corrected_diff
            if (hasmissing(centered_scores) | hasmissing(contraction)) {
                errprintf("SURE contraction exceeds finite numerical range; rescale inputs.\n")
                exit(430)
            }
            // Contract a displacement with a score before multiplying by a
            // covariance. This intermediate is dimensionless even in extreme units.
            contraction = sigma_h :* (colsum(contraction) :- (1-w[j]))
            if (hasmissing(contraction)) {
                errprintf("SURE exceeds finite numerical range; rescale inputs.\n")
                exit(430)
            }
            trace_update = sum(contraction)
        }
        else {
            H = _ewgroup_check_inverse(Id - _ewgroup_unflatten(colsum(w :* Cflat), d),
                "encountered non-invertible debiasing matrix.")
            s = (H * N')'
            corrected_diff = diff
            for (r = 1; r <= d; r++) {
                for (c = 1; c <= d; c++) {
                    corrected_diff[., r] = corrected_diff[., r] +
                        Cflat[., (r-1)*d+c] * s[c]
                }
            }
            Sigma_H = _ewgroup_unflatten(S_work[j, .], d) * H
            centered_scores = scores :- score_bar
            if (hasmissing(corrected_diff) | hasmissing(Sigma_H) | hasmissing(centered_scores)) {
                errprintf("smoothing exceeds finite numerical range; rescale inputs.\n")
                exit(430)
            }
            Q = corrected_diff' * (w :* centered_scores) - (1-w[j])*Id
            contraction = Sigma_H * Q
            if (hasmissing(Q) | hasmissing(contraction)) {
                errprintf("SURE contraction exceeds finite numerical range; rescale inputs.\n")
                exit(430)
            }
            trace_update = _ewgroup_trace(contraction)
        }
        if (hasmissing(s) | missing(trace_update)) {
            errprintf("smoothing or SURE exceeds finite numerical range; rescale inputs.\n")
            exit(430)
        }
        adjustment[j, .] = s
        trace_terms[j] = trace_update
    }

    // Evaluate the SURE ratio on a logarithmic scale. Diagnostics may underflow
    // to zero (or overflow to Stata missing) without changing the finite blend.
    adjustment_scale = max(abs(adjustment))
    trace_scale = max(abs(trace_terms))
    fit.sure_A = 0
    fit.sure_D = 0
    fit.alpha_unconstrained = 0
    fit.alpha = 0
    if (trace_scale > 0) {
        trace_sum_scaled = sum(trace_terms / trace_scale)
        if (trace_sum_scaled != 0) {
            log_abs_D = log(sigma2) + log(trace_scale) + log(abs(trace_sum_scaled))
            fit.sure_D = sign(trace_sum_scaled) * exp(log_abs_D)
        }
    }
    else trace_sum_scaled = 0
    if (adjustment_scale > 0) {
        norm_scaled = sum((adjustment / adjustment_scale):^2)
        log_A = 2*log(adjustment_scale) + log(norm_scaled)
        fit.sure_A = exp(log_A)
        if (trace_sum_scaled != 0) {
            log_ratio = log_abs_D - log_A
            fit.alpha_unconstrained = -sign(trace_sum_scaled)*exp(log_ratio)
            if (trace_sum_scaled < 0) {
                fit.alpha = (log_ratio >= 0 ? 1 : exp(log_ratio))
            }
        }
    }
    fit.theta = B + fit.alpha * adjustment
    fit.tilde = B + adjustment
    fit.derivative_trace = trace_sigma + trace_terms
    if (hasmissing(fit.theta) | hasmissing(fit.tilde) | hasmissing(fit.derivative_trace)) {
        errprintf("estimates or derivatives exceed finite numerical range; rescale inputs.\n")
        exit(430)
    }
    fit.trace_sigma = trace_sigma
    fit.gamma = gamma
    fit.sigma2 = sigma2
    fit.max_lambda = max_lambda
    return(fit)
}

// Adapter used by ewgroup_core.ado. It reads cell-level estimates from the
// Stata data set, calls the core estimator, writes optional output variables,
// and stores requested matrices/scalars back in Stata.
void _ewgroup_core_stata(string scalar bvars,
                         string scalar sigvars,
                         string scalar tousevar,
                         real scalar sigma2,
                         real scalar gamma,
                         real scalar factor,
                         real scalar covtype,
                         string scalar theta_vars,
                         string scalar tilde_vars,
                         string scalar weight_vars,
                         string scalar theta_mat,
                         string scalar tilde_mat,
                         string scalar weights_mat,
                         string scalar deriv_mat,
                         string scalar alpha_s,
                         string scalar alpha_u_s,
                         string scalar gamma_s,
                         string scalar sure_A_s,
                         string scalar sure_D_s,
                         string scalar max_lambda_s)
{
    struct ewgroup_fit scalar fit
    real colvector idx
    real matrix B, Sin, Sflat

    // `st_data()` reads variables from the Stata data set into Mata. `idx`
    // selects only observations marked for use by the ado wrapper.
    idx = selectindex(st_data(., tousevar) :!= 0)
    B = st_data(idx, tokens(bvars))
    Sin = st_data(idx, tokens(sigvars))

    // Convert the covariance variables to the common flattened matrix format.
    Sflat = _ewgroup_make_sigma_flat(Sin, cols(B), covtype)

    // Run the estimator. We request weights if either output variables or a
    // returned matrix need them.
    fit = _ewgroup_core_compute(
        B,
        Sflat,
        covtype,
        sigma2,
        gamma,
        factor,
        strlen(weight_vars) > 0 | strlen(weights_mat) > 0
    )

    // Write optional generated variables back into the Stata data set.
    if (strlen(theta_vars) > 0) {
        st_store(idx, tokens(theta_vars), fit.theta)
    }
    if (strlen(tilde_vars) > 0) {
        st_store(idx, tokens(tilde_vars), fit.tilde)
    }
    if (strlen(weight_vars) > 0) {
        st_store(idx, tokens(weight_vars), fit.weights)
    }

    // Store optional matrices for return results.
    if (strlen(theta_mat) > 0) {
        st_matrix(theta_mat, fit.theta)
    }
    if (strlen(tilde_mat) > 0) {
        st_matrix(tilde_mat, fit.tilde)
    }
    if (strlen(weights_mat) > 0) {
        st_matrix(weights_mat, fit.weights)
    }
    if (strlen(deriv_mat) > 0) {
        st_matrix(deriv_mat, fit.derivative_trace)
    }

    // Store scalar return values. Empty scalar names are ignored by
    // _ewgroup_store_scalar().
    _ewgroup_store_scalar(alpha_s, fit.alpha)
    _ewgroup_store_scalar(alpha_u_s, fit.alpha_unconstrained)
    _ewgroup_store_scalar(gamma_s, fit.gamma)
    _ewgroup_store_scalar(sure_A_s, fit.sure_A)
    _ewgroup_store_scalar(sure_D_s, fit.sure_D)
    _ewgroup_store_scalar(max_lambda_s, fit.max_lambda)
}

// Adapter used by ewgroup.ado. It starts from raw observations, estimates one
// regression inside each group, builds the covariance inputs, then calls the
// same core estimator used by ewgroup_core.ado.
void _ewgroup_raw_stata(string scalar depvar,
                        string scalar xvars,
                        string scalar groupvar,
                        string scalar tousevar,
                        real scalar noconstant,
                        real scalar sigma2_arg,
                        real scalar gamma,
                        real scalar factor,
                        string scalar theta_vars,
                        string scalar tilde_vars,
                        string scalar theta_mat,
                        string scalar beta_mat,
                        string scalar tilde_mat,
                        string scalar cells_mat,
                        string scalar ncell_mat,
                        string scalar deriv_mat,
                        string scalar alpha_s,
                        string scalar alpha_u_s,
                        string scalar gamma_s,
                        string scalar sigma2_s,
                        string scalar sure_A_s,
                        string scalar sure_D_s,
                        string scalar max_lambda_s,
                        string scalar J_s)
{
    struct ewgroup_fit scalar fit
    real colvector idx, y, g, cells, ord, e, e2, cell_n, ycell
    real matrix X, X0, Xcell, panels, B, Sflat, G, Ginv, meat, Sigma, output
    real rowvector xscale, bscaled
    real scalar N, Jn, d, j, sigma2, covtype
    string scalar group_label_name, cell_label

    // Read the selected observations from Stata into Mata.
    idx = selectindex(st_data(., tousevar) :!= 0)
    y = st_data(idx, depvar)
    g = st_data(idx, groupvar)
    N = rows(y)

    if (N < 1) {
        errprintf("no observations.\n")
        exit(2000)
    }

    // Build the regression design matrix X. With no covariates, X is just a
    // column of ones so each cell estimate is a mean. With covariates, append a
    // column of ones unless noconstant was requested.
    if (strlen(xvars) == 0) {
        if (noconstant) {
            errprintf("noconstant is not allowed when no covariates are specified.\n")
            exit(198)
        }
        X = J(N, 1, 1)
    }
    else {
        X0 = st_data(idx, tokens(xvars))
        if (noconstant) {
            X = X0
        }
        else {
            X = X0, J(N, 1, 1)
        }
    }

    // Identify the unique cells and choose sigma2. If the user did not provide
    // sigma2, use J/N, matching the package's default for the raw-data command.
    d = cols(X)
    // Sort only Mata copies. idx retains each observation's original Stata
    // row, and panel boundaries avoid rescanning the dataset for every cell.
    ord = order(g, 1)
    idx = idx[ord]
    y = y[ord]
    g = g[ord]
    X = X[ord, .]
    panels = panelsetup(g, 1)
    cells = g[panels[., 1]]
    Jn = rows(cells)
    if (sigma2_arg > 0 & !missing(sigma2_arg)) {
        sigma2 = sigma2_arg
    }
    else {
        sigma2 = Jn / N
    }

    B = J(Jn, d, .)
    Sflat = J(Jn, d * d, .)
    cell_n = J(Jn, 1, 0)
    group_label_name = st_varvaluelabel(groupvar)

    // HC0 covariance is (X'X)^-1 X'diag(e^2)X (X'X)^-1.
    // Scale design columns for the solve, then undo this scale in both outputs.
    for (j = 1; j <= Jn; j++) {
        cell_label = ""
        if (strlen(group_label_name) > 0) cell_label = st_vlmap(group_label_name, cells[j])
        if (strlen(cell_label) == 0) cell_label = strtrim(strofreal(cells[j], "%18.0g"))
        cell_n[j] = panels[j, 2] - panels[j, 1] + 1
        if (cell_n[j] <= d) {
            errprintf("cell %g (group %s) has %g observations for %g coefficients; more observations are required.\n",
                j, cell_label, cell_n[j], d)
            exit(2001)
        }
        Xcell = panelsubmatrix(X, j, panels)
        ycell = panelsubmatrix(y, j, panels)
        xscale = colmax(abs(Xcell))
        if (min(xscale) == 0) {
            errprintf("within-cell design matrix is singular in cell %g (group %s).\n", j, cell_label)
            exit(506)
        }
        Xcell = Xcell :/ xscale
        G = quadcross(Xcell, Xcell)
        Ginv = cholinv(_ewgroup_sym(G))
        if (hasmissing(Ginv)) {
            errprintf("within-cell design matrix is singular in cell %g (group %s).\n", j, cell_label)
            exit(506)
        }
        bscaled = (Ginv * quadcross(Xcell, ycell))'
        B[j, .] = bscaled :/ xscale
        e = ycell - Xcell * bscaled'
        e2 = e:^2
        if (hasmissing(e) | hasmissing(e2)) {
            errprintf("regression residuals exceed finite numerical range in cell %g (group %s).\n", j, cell_label)
            exit(430)
        }
        meat = quadcross(Xcell, e2, Xcell)
        Sigma = _ewgroup_sym(Ginv * meat * Ginv)
        Sigma = (Sigma :/ xscale :/ xscale') / sigma2
        if (hasmissing(B[j, .]) | hasmissing(Sigma)) {
            errprintf("regression exceeds finite numerical range in cell %g (group %s).\n", j, cell_label)
            exit(430)
        }
        Sflat[j, .] = _ewgroup_flatten(_ewgroup_sym(Sigma))
    }

    // Raw-data estimates always produce scalar covariance input when d is one
    // and full covariance input when there is more than one coefficient.
    covtype = (d == 1 ? 1 : 3)
    fit = _ewgroup_core_compute(B, Sflat, covtype, sigma2, gamma, factor, 0)

    // Expand output using cached panel ranges and write to original row indices.
    // Allocate just one observation-level output at a time.
    if (strlen(theta_vars) > 0) {
        output = J(N, d, .)
        for (j = 1; j <= Jn; j++) {
            output[|panels[j, 1], 1 \ panels[j, 2], d|] =
                J(cell_n[j], 1, 1) * fit.theta[j, .]
        }
        st_store(idx, tokens(theta_vars), output)
    }
    if (strlen(tilde_vars) > 0) {
        output = J(N, d, .)
        for (j = 1; j <= Jn; j++) {
            output[|panels[j, 1], 1 \ panels[j, 2], d|] =
                J(cell_n[j], 1, 1) * fit.tilde[j, .]
        }
        st_store(idx, tokens(tilde_vars), output)
    }

    // Store result matrices requested by the ado wrapper.
    if (strlen(theta_mat) > 0) {
        st_matrix(theta_mat, fit.theta)
    }
    if (strlen(beta_mat) > 0) {
        st_matrix(beta_mat, B)
    }
    if (strlen(tilde_mat) > 0) {
        st_matrix(tilde_mat, fit.tilde)
    }
    if (strlen(cells_mat) > 0) {
        st_matrix(cells_mat, cells)
    }
    if (strlen(ncell_mat) > 0) {
        st_matrix(ncell_mat, cell_n)
    }
    if (strlen(deriv_mat) > 0) {
        st_matrix(deriv_mat, fit.derivative_trace)
    }

    // Store scalar results requested by the ado wrapper.
    _ewgroup_store_scalar(alpha_s, fit.alpha)
    _ewgroup_store_scalar(alpha_u_s, fit.alpha_unconstrained)
    _ewgroup_store_scalar(gamma_s, fit.gamma)
    _ewgroup_store_scalar(sigma2_s, sigma2)
    _ewgroup_store_scalar(sure_A_s, fit.sure_A)
    _ewgroup_store_scalar(sure_D_s, fit.sure_D)
    _ewgroup_store_scalar(max_lambda_s, fit.max_lambda)
    _ewgroup_store_scalar(J_s, Jn)
}

end
