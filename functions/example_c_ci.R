# Monte Carlo confidence intervals for Example C

simulate_example_c_data <- function(n) {
  B <- matrix(
    c(
      0, 0, 0,
      0, 0, -1,
      1, 0, 0
    ),
    nrow = 3, byrow = TRUE
  )
  M <- solve(diag(3) - B)

  epsilon_e1 <- matrix(rnorm(n * 3), ncol = 3)
  A_e1 <- cbind(
    rnorm(n, 0, 1),
    rnorm(n, 0, sqrt(3)),
    rep(0, n)
  )
  Z_e1 <- t(M %*% t(epsilon_e1 + A_e1))

  epsilon_e2 <- matrix(rnorm(n * 3), ncol = 3)
  A_e2 <- cbind(
    rep(0, n),
    rnorm(n, 0, 1),
    rep(0, n)
  )
  Z_e2 <- t(M %*% t(epsilon_e2 + A_e2))

  list(
    Xe = Z_e1[, 1:2, drop = FALSE],
    ye = Z_e1[, 3],
    Xo = Z_e2[, 1:2, drop = FALSE],
    yo = Z_e2[, 3]
  )
}

estimate_gamma_path <- function(data, gamma_grid) {
  m <- build_matrices(data)
  beta_cd <- drop(compute_cd(m))
  path <- vapply(gamma_grid, function(gamma) {
    drop(solve(
      m$Gplus + gamma * m$Gdelta_plus,
      m$Zplus + gamma * m$Gdelta_plus %*% beta_cd
    ))
  }, numeric(m$p))
  list(path = path, moments = m, beta_cd = beta_cd)
}

cv_gamma_path <- function(data, gamma_grid, folds = 5) {
  sizes <- validate_cv_data(data, folds)
  fold_e <- split_idx(folds, sizes[["environment_e"]])
  fold_o <- split_idx(folds, sizes[["environment_o"]])
  fold_risk <- matrix(NA_real_, nrow = folds, ncol = length(gamma_grid))

  for (i in seq_len(folds)) {
    training <- subset_data(data, fold_e != i, fold_o != i)
    validation <- subset_data(data, fold_e == i, fold_o == i)
    training_path <- estimate_gamma_path(training, gamma_grid)$path
    validation_m <- build_matrices(validation)
    validation_cd <- drop(compute_cd(validation_m))
    beta_difference <- sweep(training_path, 1, validation_cd, "-")
    fold_risk[i, ] <- colSums(
      beta_difference * (validation_m$Gdelta_plus %*% beta_difference)
    )
  }
  colMeans(fold_risk)
}

full_sample_rdelta <- function(path_fit, data) {
  m <- path_fit$moments
  response_difference <- mean(data$ye^2) - mean(data$yo^2)
  quadratic <- colSums(path_fit$path * (m$Gdelta %*% path_fit$path))
  linear <- 2 * drop(crossprod(m$Zdelta, path_fit$path))
  quadratic - linear + response_difference
}

example_c_confidence_intervals <- function(
    n, gamma_grid, repetitions = 500, folds = 5,
    level = 0.95, seed = 2026) {
  set.seed(seed)
  p <- 2L
  n_gamma <- length(gamma_grid)
  beta_draws <- array(NA_real_, c(repetitions, p, n_gamma))
  cv_draws <- matrix(NA_real_, repetitions, n_gamma)
  rdelta_draws <- matrix(NA_real_, repetitions, n_gamma)

  for (b in seq_len(repetitions)) {
    simulated <- simulate_example_c_data(n)
    path_fit <- estimate_gamma_path(simulated, gamma_grid)
    beta_draws[b, , ] <- path_fit$path
    cv_draws[b, ] <- cv_gamma_path(simulated, gamma_grid, folds)
    rdelta_draws[b, ] <- full_sample_rdelta(path_fit, simulated)
  }

  alpha <- (1 - level) / 2
  interval <- function(x) {
    apply(x, 2, quantile, probs = c(alpha, 1 - alpha), na.rm = TRUE)
  }
  beta_lower <- beta_upper <- matrix(NA_real_, p, n_gamma)
  for (j in seq_len(p)) {
    bounds <- interval(beta_draws[, j, ])
    beta_lower[j, ] <- bounds[1, ]
    beta_upper[j, ] <- bounds[2, ]
  }
  cv_bounds <- interval(cv_draws)
  rdelta_bounds <- interval(rdelta_draws)

  list(
    repetitions = repetitions,
    level = level,
    beta_lower = beta_lower,
    beta_upper = beta_upper,
    cv_lower = cv_bounds[1, ],
    cv_upper = cv_bounds[2, ],
    rdelta_lower = rdelta_bounds[1, ],
    rdelta_upper = rdelta_bounds[2, ]
  )
}
