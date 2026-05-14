version 16.0

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

mata:

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

real matrix _ewgroup_sym(real matrix A)
{
    return((A + A') / 2)
}

real scalar _ewgroup_trace(real matrix A)
{
    real scalar i, out

    out = 0
    for (i = 1; i <= min((rows(A), cols(A))); i++) {
        out = out + A[i, i]
    }
    return(out)
}

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

real matrix _ewgroup_check_inverse(real matrix A, string scalar message)
{
    real matrix S

    S = _ewgroup_sym(A)
    if (rank(S) < rows(S)) {
        errprintf("%s\n", message)
        exit(506)
    }
    return(invsym(S))
}

real scalar _ewgroup_max_lambda(real matrix S)
{
    real colvector evals

    evals = symeigenvalues(_ewgroup_sym(S))
    return(max(evals))
}

void _ewgroup_store_scalar(string scalar name, real scalar value)
{
    if (strlen(name) > 0) {
        st_numscalar(name, value)
    }
}

real matrix _ewgroup_make_sigma_flat(real matrix Sin, real scalar d,
                                     real scalar covtype)
{
    real scalar Jn, j, r, c, p
    real matrix Sflat, S

    Jn = rows(Sin)
    Sflat = J(Jn, d * d, 0)

    if (covtype == 1) {
        for (j = 1; j <= Jn; j++) {
            Sflat[j, 1] = Sin[j, 1]
        }
        return(Sflat)
    }

    if (covtype == 2) {
        for (j = 1; j <= Jn; j++) {
            for (r = 1; r <= d; r++) {
                Sflat[j, (r - 1) * d + r] = Sin[j, r]
            }
        }
        return(Sflat)
    }

    for (j = 1; j <= Jn; j++) {
        S = J(d, d, 0)
        p = 1
        for (r = 1; r <= d; r++) {
            for (c = 1; c <= d; c++) {
                S[r, c] = Sin[j, p]
                p = p + 1
            }
        }
        Sflat[j, .] = _ewgroup_flatten(_ewgroup_sym(S))
    }

    return(Sflat)
}

struct ewgroup_fit scalar _ewgroup_core_compute(real matrix B,
                                                real matrix Sflat,
                                                real scalar covtype,
                                                real scalar sigma2,
                                                real scalar gamma,
                                                real scalar factor,
                                                real scalar want_weights)
{
    struct ewgroup_fit scalar fit
    real scalar Jn, d, j, k, r, rowmax, rowsum, qform
    real scalar trace_update, sure_A, sure_D, max_lambda
    real matrix Id, S, Omega, H, Gamma, Sigma_H, weights_out
    real matrix Omega_flat, theta, tilde
    real rowvector bj, diff, logw, w, weighted_mean, N, s
    real rowvector score_bar, q, transformed_row
    real matrix scores, Sigma_s_all
    real colvector diff_col, transformed, evals
    real colvector trace_terms, trace_sigma

    Jn = rows(B)
    d = cols(B)

    if (Jn < 1 | d < 1) {
        errprintf("beta estimates must have at least one row and one column.\n")
        exit(198)
    }
    if (sigma2 <= 0 | missing(sigma2)) {
        errprintf("sigma2 must be a positive scalar.\n")
        exit(198)
    }
    if (factor <= 0 | missing(factor)) {
        errprintf("factor must be a positive scalar.\n")
        exit(198)
    }

    Id = I(d)
    Omega_flat = J(Jn, d * d, 0)
    trace_sigma = J(Jn, 1, 0)
    max_lambda = .

    for (k = 1; k <= Jn; k++) {
        S = _ewgroup_unflatten(Sflat[k, .], d)
        S = _ewgroup_sym(S)
        evals = symeigenvalues(S)
        if (min(evals) < -sqrt(epsilon(1)) * max((1, max(abs(S))))) {
            errprintf("Sigma_hat must be positive semidefinite.\n")
            exit(506)
        }
        trace_sigma[k] = _ewgroup_trace(S)
        if (k == 1) {
            max_lambda = max(evals)
        }
        else {
            max_lambda = max((max_lambda, max(evals)))
        }
    }

    if (gamma <= 0 | missing(gamma)) {
        if (max_lambda <= 0 | missing(max_lambda)) {
            errprintf("Sigma_hat must have a positive eigenvalue to compute gamma.\n")
            exit(198)
        }
        gamma = factor / max_lambda
    }
    if (gamma <= 0 | missing(gamma)) {
        errprintf("gamma must be a positive scalar.\n")
        exit(198)
    }
    if (gamma * max_lambda >= 1) {
        errprintf("gamma must satisfy gamma * max_j lambda_max(Sigma_hat_j) < 1.\n")
        exit(198)
    }

    for (k = 1; k <= Jn; k++) {
        S = _ewgroup_unflatten(Sflat[k, .], d)
        Omega = gamma * _ewgroup_check_inverse(
            Id - gamma * S,
            "gamma times the largest covariance eigenvalue must be below one."
        )
        Omega_flat[k, .] = _ewgroup_flatten(_ewgroup_sym(Omega))
    }

    theta = J(Jn, d, .)
    tilde = J(Jn, d, .)
    trace_terms = J(Jn, 1, 0)
    if (want_weights) {
        weights_out = J(Jn, Jn, .)
    }
    else {
        weights_out = J(0, 0, .)
    }

    for (j = 1; j <= Jn; j++) {
        bj = B[j, .]
        logw = J(1, Jn, 0)
        scores = J(Jn, d, 0)

        for (k = 1; k <= Jn; k++) {
            diff = B[k, .] - bj
            diff_col = diff'
            Omega = _ewgroup_unflatten(Omega_flat[k, .], d)
            transformed = Omega * diff_col
            qform = diff * transformed
            logw[k] = -0.5 * qform / sigma2
            scores[k, .] = (transformed / sigma2)'
        }

        rowmax = max(logw)
        w = exp(logw :- rowmax)
        rowsum = sum(w)
        if (rowsum <= 0 | missing(rowsum)) {
            errprintf("encountered invalid exponential weights.\n")
            exit(430)
        }
        w = w / rowsum
        if (want_weights) {
            weights_out[j, .] = w
        }

        weighted_mean = J(1, d, 0)
        Gamma = J(d, d, 0)
        for (k = 1; k <= Jn; k++) {
            S = _ewgroup_unflatten(Sflat[k, .], d)
            weighted_mean = weighted_mean + w[k] * B[k, .]
            Gamma = Gamma + gamma * w[k] * S
        }

        H = _ewgroup_check_inverse(Id - Gamma, "Encountered non-invertible debiasing matrix.")
        N = weighted_mean - bj
        s = (H * N')'
        tilde[j, .] = bj + s

        score_bar = w * scores
        Sigma_s_all = J(Jn, d, 0)
        for (k = 1; k <= Jn; k++) {
            S = _ewgroup_unflatten(Sflat[k, .], d)
            Sigma_s_all[k, .] = (S * s')'
        }

        S = _ewgroup_unflatten(Sflat[j, .], d)
        Sigma_H = S * H
        trace_update = -(1 - w[j]) * _ewgroup_trace(Sigma_H)
        for (k = 1; k <= Jn; k++) {
            diff = B[k, .] - bj + gamma * Sigma_s_all[k, .]
            transformed = Sigma_H * diff'
            transformed_row = transformed'
            q = w[k] * (scores[k, .] - score_bar)
            for (r = 1; r <= d; r++) {
                trace_update = trace_update + q[r] * transformed_row[r]
            }
        }
        trace_terms[j] = trace_update
    }

    sure_A = 0
    for (j = 1; j <= Jn; j++) {
        diff = tilde[j, .] - B[j, .]
        sure_A = sure_A + diff * diff'
    }
    sure_D = sigma2 * sum(trace_terms)

    fit.alpha_unconstrained = (sure_A <= epsilon(1) ? 0 : -sure_D / sure_A)
    fit.alpha = min((1, max((0, fit.alpha_unconstrained))))
    fit.theta = fit.alpha * tilde + (1 - fit.alpha) * B
    fit.tilde = tilde
    fit.weights = weights_out
    fit.derivative_trace = trace_sigma + trace_terms
    fit.trace_sigma = trace_sigma
    fit.gamma = gamma
    fit.sigma2 = sigma2
    fit.sure_A = sure_A
    fit.sure_D = sure_D
    fit.max_lambda = max_lambda

    return(fit)
}

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

    idx = selectindex(st_data(., tousevar) :!= 0)
    B = st_data(idx, tokens(bvars))
    Sin = st_data(idx, tokens(sigvars))
    Sflat = _ewgroup_make_sigma_flat(Sin, cols(B), covtype)

    fit = _ewgroup_core_compute(
        B,
        Sflat,
        covtype,
        sigma2,
        gamma,
        factor,
        strlen(weight_vars) > 0 | strlen(weights_mat) > 0
    )

    if (strlen(theta_vars) > 0) {
        st_store(idx, tokens(theta_vars), fit.theta)
    }
    if (strlen(tilde_vars) > 0) {
        st_store(idx, tokens(tilde_vars), fit.tilde)
    }
    if (strlen(weight_vars) > 0) {
        st_store(idx, tokens(weight_vars), fit.weights)
    }

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

    _ewgroup_store_scalar(alpha_s, fit.alpha)
    _ewgroup_store_scalar(alpha_u_s, fit.alpha_unconstrained)
    _ewgroup_store_scalar(gamma_s, fit.gamma)
    _ewgroup_store_scalar(sure_A_s, fit.sure_A)
    _ewgroup_store_scalar(sure_D_s, fit.sure_D)
    _ewgroup_store_scalar(max_lambda_s, fit.max_lambda)
}

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
                        string scalar max_lambda_s)
{
    struct ewgroup_fit scalar fit
    real colvector idx, y, g, cells, sel, e, cell_n
    real matrix X, X0, B, Sflat, G, Ginv, meat, Sigma, theta_obs, tilde_obs
    real scalar N, Jn, d, j, i, sigma2, covtype

    idx = selectindex(st_data(., tousevar) :!= 0)
    y = st_data(idx, depvar)
    g = st_data(idx, groupvar)
    N = rows(y)

    if (N < 1) {
        errprintf("no observations.\n")
        exit(2000)
    }

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

    d = cols(X)
    cells = uniqrows(sort(g, 1))
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

    for (j = 1; j <= Jn; j++) {
        sel = selectindex(g :== cells[j])
        cell_n[j] = rows(sel)
        if (cell_n[j] < d) {
            errprintf("cell has fewer observations than coefficients.\n")
            exit(2001)
        }

        G = quadcross(X[sel, .], X[sel, .])
        if (rank(G) < d) {
            errprintf("within-cell design matrix is singular.\n")
            exit(506)
        }
        Ginv = invsym(_ewgroup_sym(G))
        B[j, .] = (Ginv * quadcross(X[sel, .], y[sel]))'
        e = y[sel] - X[sel, .] * B[j, .]'

        meat = J(d, d, 0)
        for (i = 1; i <= rows(sel); i++) {
            meat = meat + e[i]^2 * X[sel[i], .]' * X[sel[i], .]
        }

        Sigma = (Ginv * meat * Ginv) / sigma2
        Sflat[j, .] = _ewgroup_flatten(_ewgroup_sym(Sigma))
    }

    covtype = (d == 1 ? 1 : 3)
    fit = _ewgroup_core_compute(B, Sflat, covtype, sigma2, gamma, factor, 0)

    if (strlen(theta_vars) > 0 | strlen(tilde_vars) > 0) {
        theta_obs = J(N, d, .)
        tilde_obs = J(N, d, .)
        for (j = 1; j <= Jn; j++) {
            sel = selectindex(g :== cells[j])
            theta_obs[sel, .] = J(rows(sel), 1, 1) * fit.theta[j, .]
            tilde_obs[sel, .] = J(rows(sel), 1, 1) * fit.tilde[j, .]
        }
        if (strlen(theta_vars) > 0) {
            st_store(idx, tokens(theta_vars), theta_obs)
        }
        if (strlen(tilde_vars) > 0) {
            st_store(idx, tokens(tilde_vars), tilde_obs)
        }
    }

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

    _ewgroup_store_scalar(alpha_s, fit.alpha)
    _ewgroup_store_scalar(alpha_u_s, fit.alpha_unconstrained)
    _ewgroup_store_scalar(gamma_s, fit.gamma)
    _ewgroup_store_scalar(sigma2_s, sigma2)
    _ewgroup_store_scalar(sure_A_s, fit.sure_A)
    _ewgroup_store_scalar(sure_D_s, fit.sure_D)
    _ewgroup_store_scalar(max_lambda_s, fit.max_lambda)
}

end
