get_script_path <- function() {
  # `commandArgs(FALSE)` returns the command-line arguments used to start R. When
  # this script is run with Rscript, one argument usually starts with "--file=".
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) == 0L) {
    # Fallback for interactive runs from the tests directory.
    return(normalizePath("make_reference.R"))
  }

  # Remove the "--file=" prefix and turn the result into a full path.
  normalizePath(sub("^--file=", "", file_arg[1L]))
}

# This file is an R script, so comments use # below.
script_path <- get_script_path()
test_dir <- dirname(script_path)
repo_root <- normalizePath(file.path(test_dir, "..", ".."))

# Run against the installed compiled package and the pure-R implementation.
# --library= selects an isolated installation; --output= selects a fixture directory.
args <- commandArgs(trailingOnly = TRUE)
arg <- function(name, default) {
  hit <- args[startsWith(args, paste0("--", name, "="))]
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[[1L]]) else default
}
lib <- arg("library", "")
if (nzchar(lib)) .libPaths(c(normalizePath(lib), .libPaths()))
if (!requireNamespace("ewgroup", quietly = TRUE)) {
  stop("Install the companion R package first; optionally pass --library=/path/to/library.")
}
out_dir <- arg("output", test_dir)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fallback <- new.env(parent = baseenv())
sys.source(file.path(repo_root, "ewgroup", "R", "ewgroup.R"), envir = fallback)
source(file.path(test_dir, "paper_reference.R"))

# Every saved fixture must agree across installed C++, pure R, and a small
# implementation written directly from the draft's equations.
ewgroup <- function(beta_hat, Sigma_hat, sigma2, gamma = NULL,
                    return_weights = FALSE, return_derivative = FALSE) {
  compiled <- ewgroup::ewgroup(beta_hat, Sigma_hat, sigma2, gamma,
                              return_weights, return_derivative)
  pure_r <- fallback$ewgroup(beta_hat, Sigma_hat, sigma2, gamma,
                            return_weights, return_derivative)
  oracle <- paper_reference(beta_hat, Sigma_hat, sigma2, compiled$gamma)
  for (key in c("theta", "tilde", "alpha", "alpha_unconstrained", "sure_A", "sure_D")) {
    actual <- as.numeric(compiled[[key]])
    for (expected in list(pure_r[[key]], oracle[[key]])) {
      if (!isTRUE(all.equal(actual, as.numeric(expected), tolerance = 1e-9))) {
        stop("Fixture disagreement for ", key, ": installed C++, R, or paper reference.")
      }
    }
  }
  if (return_weights) stopifnot(isTRUE(all.equal(unname(compiled$weights), oracle$weights,
                                               tolerance = 1e-10)))
  if (return_derivative) {
    got <- if (is.list(compiled$derivative)) unlist(compiled$derivative) else compiled$derivative
    stopifnot(isTRUE(all.equal(as.numeric(got), unlist(oracle$derivative), tolerance = 1e-9)))
  }
  compiled
}

write_csv <- function(x, file) {
  # `utils::write.csv()` is the standard R function for writing a data frame to a
  # CSV file. The `utils::` prefix makes clear that the function comes from R's
  # built-in utils package.
  utils::write.csv(x, file = file.path(out_dir, file), row.names = FALSE)
}

write_scalars <- function(x, file) {
  # Stata reads tabular data most easily, so turn a named numeric vector into two
  # columns: the scalar name and its value.
  out <- data.frame(
    name = names(x),
    value = as.numeric(x),
    stringsAsFactors = FALSE
  )
  write_csv(out, file)
}

# Scalar reference fixture: five cell estimates and one scaled covariance number
# per cell.
beta_hat <- c(-1.1, -0.95, 0.2, 0.27, 1.3)
Sigma_hat <- c(0.8, 1.1, 0.9, 1.2, 1.0)
fit_scalar <- ewgroup(
  beta_hat,
  Sigma_hat,
  sigma2 = 0.05,
  gamma = 0.2,
  return_weights = TRUE,
  return_derivative = TRUE
)

# Store the values that Stata should reproduce for the scalar fixture.
scalar_data <- data.frame(
  beta_hat = beta_hat,
  sigma_hat = Sigma_hat,
  theta = as.numeric(fit_scalar$theta),
  tilde = as.numeric(fit_scalar$tilde),
  derivative = as.numeric(fit_scalar$derivative),
  fit_scalar$weights,
  check.names = FALSE
)
names(scalar_data)[6:ncol(scalar_data)] <- paste0("w", seq_along(beta_hat))
write_csv(scalar_data, "reference_scalar.csv")

# Store scalar diagnostics separately because the certification script compares
# them through returned r() scalars.
write_scalars(
  c(
    alpha = fit_scalar$alpha,
    alpha_unconstrained = fit_scalar$alpha_unconstrained,
    gamma = fit_scalar$gamma,
    sure_A = fit_scalar$sure_A,
    sure_D = fit_scalar$sure_D
  ),
  "reference_scalar_scalars.csv"
)

# Vector-valued reference fixture: four cells and two coefficients per cell.
B <- rbind(
  c(-1.0, 0.2),
  c(-0.8, 0.1),
  c(0.6, -0.3),
  c(0.7, -0.4)
)
Sigma <- list(
  matrix(c(1.0, 0.1, 0.1, 1.3), 2),
  matrix(c(1.1, 0.0, 0.0, 0.9), 2),
  matrix(c(0.8, 0.2, 0.2, 1.4), 2),
  matrix(c(1.2, 0.1, 0.1, 1.0), 2)
)
fit_vector <- ewgroup(
  B,
  Sigma,
  sigma2 = 0.1,
  gamma = 0.05,
  return_weights = FALSE,
  return_derivative = FALSE
)

# Pull covariance entries into named columns because Stata's sigma() option
# receives variables, not R list objects.
vector_data <- data.frame(
  b1 = B[, 1],
  b2 = B[, 2],
  s11 = vapply(Sigma, function(x) x[1, 1], numeric(1)),
  s12 = vapply(Sigma, function(x) x[1, 2], numeric(1)),
  s21 = vapply(Sigma, function(x) x[2, 1], numeric(1)),
  s22 = vapply(Sigma, function(x) x[2, 2], numeric(1)),
  theta1 = fit_vector$theta[, 1],
  theta2 = fit_vector$theta[, 2],
  tilde1 = fit_vector$tilde[, 1],
  tilde2 = fit_vector$tilde[, 2]
)
write_csv(vector_data, "reference_vector.csv")

# Store diagnostics for the vector fixture.
write_scalars(
  c(
    alpha = fit_vector$alpha,
    alpha_unconstrained = fit_vector$alpha_unconstrained,
    gamma = fit_vector$gamma,
    sure_A = fit_vector$sure_A,
    sure_D = fit_vector$sure_D
  ),
  "reference_vector_scalars.csv"
)

# Raw-data reference fixture. Stata's ewgroup command starts from observations,
# so this fixture includes the original group and outcome columns.
raw <- data.frame(
  group = rep(1:4, each = 5),
  y = c(
    -1.25, -1.10, -1.05, -0.95, -0.90,
    -0.85, -0.78, -0.72, -0.70, -0.65,
    0.35, 0.42, 0.55, 0.60, 0.63,
    1.05, 1.10, 1.22, 1.28, 1.35
  )
)
J <- length(unique(raw$group))
N <- nrow(raw)
sigma2_raw <- J / N

# Recreate the cell-level mean and scaled covariance calculations used by the
# Stata raw-data command, then run the R estimator on those cell-level values.
beta_raw <- as.numeric(tapply(raw$y, raw$group, mean))
Sigma_raw <- as.numeric(tapply(raw$y, raw$group, function(z) {
  sum((z - mean(z))^2) / length(z)^2 / sigma2_raw
}))
fit_raw <- ewgroup(beta_raw, Sigma_raw, sigma2 = sigma2_raw)

# Expand the cell-level estimates back to the observation rows so Stata can
# compare generated variables directly.
raw$theta <- fit_raw$theta[raw$group]
raw$tilde <- fit_raw$tilde[raw$group]
write_csv(raw, "reference_raw.csv")

# Store diagnostics for the raw-data fixture.
write_scalars(
  c(
    alpha = fit_raw$alpha,
    alpha_unconstrained = fit_raw$alpha_unconstrained,
    gamma = fit_raw$gamma,
    sigma2 = sigma2_raw,
    sure_A = fit_raw$sure_A,
    sure_D = fit_raw$sure_D
  ),
  "reference_raw_scalars.csv"
)

# General fixtures include all covariance paths, SURE boundaries, and changes
# of measurement units. Explicitly store all covariances in row-major order.
set.seed(904)
cases <- list()
add_case <- function(name, B, S, sigma2 = .15, gamma = NULL) {
  B <- as.matrix(B)
  d <- ncol(B)
  fit <- ewgroup(B, S, sigma2, gamma, return_derivative = TRUE)
  oracle <- paper_reference(B, S, sigma2, fit$gamma)
  out <- data.frame(case = name, cell = seq_len(nrow(B)), d = d,
                    sigma2 = sigma2, gamma = fit$gamma, alpha = fit$alpha,
                    sure_A = fit$sure_A, sure_D = fit$sure_D,
                    derivative_trace = oracle$derivative_trace)
  for (k in 1:3) {
    out[[paste0("b", k)]] <- if (k <= d) B[, k] else 0
    out[[paste0("theta", k)]] <- if (k <= d) fit$theta[, k] else 0
    out[[paste0("tilde", k)]] <- if (k <= d) fit$tilde[, k] else 0
    for (ell in 1:3) out[[paste0("s", k, ell)]] <- if (k <= d && ell <= d) {
      vapply(S, function(s) s[k, ell], numeric(1))
    } else 0
  }
  cases[[length(cases) + 1L]] <<- out
}
y <- c(-1.1, -.95, .2, .27, 1.3)
s <- c(.8, 1.1, .9, 1.2, 1)
add_case("scalar_interior", y, lapply(s, matrix), .05)
add_case("scalar_alpha_one", y, lapply(s, matrix), .05, .2)
add_case("scalar_alpha_zero", c(-3, 0, 3), rep(list(matrix(1)), 3), 1, .2)
add_case("tiny_units", y * 1e-8, lapply(s * 1e-16, matrix), .05)
add_case("large_units", y * 1e8, lapply(s * 1e16, matrix), .05)
add_case("single", matrix(c(.3, -.7), 1), list(diag(c(.5, 1))))
add_case("identical", matrix(rep(c(.25, -.5), each = 5), 5), rep(list(diag(2)), 5))
B3 <- matrix(rnorm(24), 8, 3)
S3 <- lapply(1:8, function(j) { A <- matrix(rnorm(9), 3); crossprod(A) / 3 + diag(.2, 3) })
add_case("full_three", B3, S3)
add_case("diagonal_three", B3, lapply(S3, function(s) diag(diag(s))))
add_case("zero_covariance", y, rep(list(matrix(0)), 5), 1, .2)
write_csv(do.call(rbind, cases), "reference_cases.csv")

# Independently estimate raw-data OLS and HC0 with QR least squares. Include
# unsorted, nonconsecutive groups and both intercept conventions.
raw_reg <- data.frame(group = rep(c(17, 3, 40, 9), each = 9),
                      x = rep(seq(-1, 1, length.out = 9), 4))
raw_reg$y <- .2 * raw_reg$group + (.4 + raw_reg$group / 20) * raw_reg$x +
  .3 * sin(seq_len(nrow(raw_reg)))
raw_reg <- raw_reg[sample(seq_len(nrow(raw_reg))), ]
raw_reg$row_id <- seq_len(nrow(raw_reg))
for (intercept in c(TRUE, FALSE)) {
  ids <- sort(unique(raw_reg$group)); d <- if (intercept) 2L else 1L
  sigma2 <- length(ids) / nrow(raw_reg)
  X <- if (intercept) cbind(raw_reg$x, 1) else matrix(raw_reg$x, ncol = 1)
  B <- matrix(0, length(ids), d); S <- vector("list", length(ids))
  for (j in seq_along(ids)) {
    sel <- raw_reg$group == ids[j]
    Xj <- X[sel, , drop = FALSE]
    fit <- lm.fit(Xj, raw_reg$y[sel])
    B[j, ] <- fit$coefficients
    invG <- solve(crossprod(Xj))
    S[[j]] <- invG %*% crossprod(Xj * fit$residuals) %*% invG / sigma2
  }
  fit <- ewgroup(B, S, sigma2)
  out <- raw_reg
  pos <- match(out$group, ids)
  for (k in seq_len(d)) {
    out[[paste0("theta", k)]] <- fit$theta[pos, k]
    out[[paste0("tilde", k)]] <- fit$tilde[pos, k]
    out[[paste0("beta", k)]] <- B[pos, k]
  }
  out$alpha <- fit$alpha
  out$sigma2 <- sigma2
  write_csv(out, if (intercept) "reference_raw_regression.csv" else "reference_raw_noconstant.csv")
}
cat("Fixtures checked against installed C++, pure R, and paper equations; written to ",
    normalizePath(out_dir), "\n", sep = "")
