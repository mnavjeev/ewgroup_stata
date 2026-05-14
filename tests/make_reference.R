get_script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) == 0L) {
    return(normalizePath("make_reference.R"))
  }
  normalizePath(sub("^--file=", "", file_arg[1L]))
}

script_path <- get_script_path()
test_dir <- dirname(script_path)
repo_root <- normalizePath(file.path(test_dir, "..", ".."))

source(file.path(repo_root, "ewgroup", "R", "ewgroup.R"))

write_csv <- function(x, file) {
  utils::write.csv(x, file = file.path(test_dir, file), row.names = FALSE)
}

write_scalars <- function(x, file) {
  out <- data.frame(
    name = names(x),
    value = as.numeric(x),
    stringsAsFactors = FALSE
  )
  write_csv(out, file)
}

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
beta_raw <- as.numeric(tapply(raw$y, raw$group, mean))
Sigma_raw <- as.numeric(tapply(raw$y, raw$group, function(z) {
  sum((z - mean(z))^2) / length(z)^2 / sigma2_raw
}))
fit_raw <- ewgroup(beta_raw, Sigma_raw, sigma2 = sigma2_raw)
raw$theta <- fit_raw$theta[raw$group]
raw$tilde <- fit_raw$tilde[raw$group]
write_csv(raw, "reference_raw.csv")
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
