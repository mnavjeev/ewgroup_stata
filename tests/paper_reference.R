# Deliberately small reference implementation of the draft's equations.
# It uses full matrices and solves, and is only used for small test fixtures.
paper_reference <- function(beta_hat, Sigma_hat, sigma2, gamma) {
  B <- as.matrix(beta_hat)
  J <- nrow(B)
  d <- ncol(B)
  S <- if (is.list(Sigma_hat)) Sigma_hat else if (length(dim(Sigma_hat)) == 3L) {
    lapply(seq_len(J), function(j) matrix(Sigma_hat[, , j], d, d))
  } else lapply(Sigma_hat, function(s) matrix(s, 1L, 1L))
  Id <- diag(d)
  Omega <- lapply(S, function(s) gamma * solve(Id - gamma * s))
  W <- matrix(0, J, J)
  displacement <- matrix(0, J, d)
  P <- vector("list", J)
  trace <- numeric(J)
  for (j in seq_len(J)) {
    delta <- sweep(B, 2L, B[j, ], "-")
    scores <- t(matrix(vapply(seq_len(J), function(k) {
      as.numeric(Omega[[k]] %*% delta[k, ]) / sigma2
    }, numeric(d)), nrow = d, ncol = J))
    logits <- vapply(seq_len(J), function(k) {
      -as.numeric(crossprod(delta[k, ], Omega[[k]] %*% delta[k, ])) / (2 * sigma2)
    }, numeric(1L))
    w <- exp(logits - max(logits)); w <- w / sum(w)
    W[j, ] <- w
    Gamma <- Reduce(`+`, Map(function(wk, sk) gamma * wk * sk, w, S))
    H <- solve(Id - Gamma)
    move <- as.numeric(H %*% as.numeric(crossprod(w, delta)))
    displacement[j, ] <- move
    centered_score <- sweep(scores, 2L, as.numeric(crossprod(w, scores)), "-")
    update <- matrix(0, d, d)
    for (k in seq_len(J)) {
      update <- update + w[k] * tcrossprod(
        delta[k, ] + as.numeric(gamma * S[[k]] %*% move), centered_score[k, ])
    }
    P[[j]] <- Id - (1 - w[j]) * H + H %*% update
    trace[j] <- sum(diag(S[[j]] %*% P[[j]]))
  }
  A <- sum(displacement^2)
  D <- sigma2 * sum(vapply(seq_len(J), function(j) {
    sum(diag(S[[j]] %*% (P[[j]] - Id)))
  }, numeric(1L)))
  au <- if (A == 0) 0 else -D / A
  alpha <- max(0, min(1, au))
  list(theta = B + alpha * displacement, tilde = B + displacement,
       alpha = alpha, alpha_unconstrained = au, sure_A = A, sure_D = D,
       weights = W, derivative = P, derivative_trace = trace)
}
