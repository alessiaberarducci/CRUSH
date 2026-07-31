# Functions for the robust RPE analysis

# Retain environments with enough cells for a stable comparison.
filter_environments <- function(data, environment_column, minimum_cells = 100) {
  environment <- as.character(data[[environment_column]])
  counts <- table(environment)
  keep <- names(counts)[counts >= minimum_cells]
  data[environment %in% setdiff(keep, "excluded"), , drop = FALSE]
}

make_two_environment_data <- function(
    data, environment_column, first_environments, second_environments,
    predictor_genes, target_gene) {
  # Pool several distinct interventions into each training environment.
  overlap <- intersect(first_environments, second_environments)
  if (length(overlap) > 0L) {
    stop("The two training environments contain overlapping interventions.")
  }

  available <- unique(as.character(data[[environment_column]]))
  missing <- setdiff(c(first_environments, second_environments), available)
  if (length(missing) > 0L) {
    stop("Training interventions not found: ", paste(missing, collapse = ", "))
  }

  first <- data[
    data[[environment_column]] %in% first_environments, , drop = FALSE
  ]
  second <- data[
    data[[environment_column]] %in% second_environments, , drop = FALSE
  ]
  list(
    data = list(
      Xe = as.matrix(first[, predictor_genes, drop = FALSE]),
      ye = as.numeric(first[[target_gene]]),
      Xo = as.matrix(second[, predictor_genes, drop = FALSE]),
      yo = as.numeric(second[[target_gene]])
    ),
    first = first,
    second = second,
    first_environments = first_environments,
    second_environments = second_environments
  )
}

compute_common_centering <- function(data) {
  # Equal environment weights match the definition of Rplus.
  x_center <- 0.5 * (
    colMeans(data$Xe) + colMeans(data$Xo)
  )
  y_center <- 0.5 * (mean(data$ye) + mean(data$yo))
  list(
    data = list(
      Xe = sweep(data$Xe, 2, x_center, "-"),
      ye = data$ye - y_center,
      Xo = sweep(data$Xo, 2, x_center, "-"),
      yo = data$yo - y_center
    ),
    x_center = x_center,
    y_center = y_center
  )
}

apply_common_centering <- function(data, centering) {
  list(
    Xe = sweep(data$Xe, 2, centering$x_center, "-"),
    ye = data$ye - centering$y_center,
    Xo = sweep(data$Xo, 2, centering$x_center, "-"),
    yo = data$yo - centering$y_center
  )
}

centered_cross_validation <- function(
    data, estimators, fold_assignment,
    validation_loss = c("ours", "kania")) {
  # Estimate centering moments inside each fold to avoid validation leakage.
  validation_loss <- match.arg(validation_loss)
  fold_e <- fold_assignment$environment_e
  fold_o <- fold_assignment$environment_o
  folds <- length(unique(fold_e))
  gamma_values <- as.numeric(names(estimators))
  foldwise <- matrix(
    NA_real_, nrow = folds, ncol = length(gamma_values),
    dimnames = list(NULL, names(estimators))
  )

  for (fold in seq_len(folds)) {
    raw_training <- subset_data(data, fold_e != fold, fold_o != fold)
    raw_validation <- subset_data(data, fold_e == fold, fold_o == fold)
    centering <- compute_common_centering(raw_training)
    training <- centering$data
    validation <- apply_common_centering(raw_validation, centering)

    # Compute fold moments once and reuse them along the complete gamma path.
    training_moments <- build_matrices(training)
    beta_cp <- drop(compute_cd(training_moments))
    beta_path <- vapply(gamma_values, function(gamma) {
      drop(solve(
        training_moments$Gplus +
          gamma * training_moments$Gdelta_plus,
        training_moments$Zplus +
          gamma * training_moments$Gdelta_plus %*% beta_cp
      ))
    }, numeric(training_moments$p))

    if (validation_loss == "ours") {
      validation_moments <- build_matrices(validation)
      validation_cp <- drop(compute_cd(validation_moments))
      beta_difference <- sweep(beta_path, 1, validation_cp, "-")
      foldwise[fold, ] <- colSums(
        beta_difference *
          (validation_moments$Gdelta_plus %*% beta_difference)
      )
    } else {
      residual_e <- sweep(
        validation$Xe %*% beta_path, 1, validation$ye, "-"
      )
      residual_o <- sweep(
        validation$Xo %*% beta_path, 1, validation$yo, "-"
      )
      foldwise[fold, ] <- abs(
        colMeans(residual_e^2) - colMeans(residual_o^2)
      )
    }
  }

  mean_loss <- colMeans(foldwise)
  best <- which.min(mean_loss)
  list(
    best_gamma = names(estimators)[best],
    mean_loss = mean_loss,
    foldwise_loss = foldwise,
    fold_assignment = fold_assignment
  )
}

recover_intercept <- function(beta, centering) {
  drop(centering$y_center - crossprod(centering$x_center, beta))
}

compute_environment_risk <- function(
    data, beta, predictor_genes, target_gene, intercept = 0) {
  # Missing observations are removed within each test environment.
  X <- as.matrix(data[, predictor_genes, drop = FALSE])
  y <- as.numeric(data[[target_gene]])
  complete <- complete.cases(X, y)
  if (!any(complete)) return(NA_real_)
  prediction <- intercept +
    X[complete, , drop = FALSE] %*% beta[predictor_genes]
  mean((y[complete] - prediction)^2)
}

compute_test_risks <- function(
    data, environments, environment_column, beta_ols, beta_cr,
    predictor_genes, target_gene, intercept_ols = 0, intercept_cr = 0) {
  # Evaluate OLS and CR on the same unseen environments.
  result <- lapply(environments, function(environment) {
    current <- data[data[[environment_column]] == environment, , drop = FALSE]
    data.frame(
      environment = environment,
      n_cells = nrow(current),
      OLS = compute_environment_risk(
        current, beta_ols, predictor_genes, target_gene, intercept_ols
      ),
      CR = compute_environment_risk(
        current, beta_cr, predictor_genes, target_gene, intercept_cr
      )
    )
  })
  result <- do.call(rbind, result)
  result[is.finite(result$OLS) & is.finite(result$CR), , drop = FALSE]
}

summarise_test_risks <- function(risks) {
  # Include the worst observed environment together with standard summaries.
  data.frame(
    method = c("OLS", "CR"),
    mean = c(mean(risks$OLS), mean(risks$CR)),
    median = c(median(risks$OLS), median(risks$CR)),
    sd = c(sd(risks$OLS), sd(risks$CR)),
    worst = c(max(risks$OLS), max(risks$CR))
  )
}

compare_training_risks <- function(data, estimates, gamma) {
  # Compare Rdelta+, Rplus and their weighted objective Rgamma.
  rows <- lapply(names(estimates), function(method) {
    beta <- estimates[[method]]
    r_delta <- calc_Rdelta(data, beta)
    r_plus <- in_sample_risk(data, beta)
    data.frame(
      method = method,
      Rdelta_plus = r_delta,
      Rplus = r_plus,
      Rgamma = 0.5 * r_plus + 0.5 * gamma * r_delta
    )
  })
  do.call(rbind, rows)
}

nested_cross_validation <- function(
    data, estimators, outer_folds = 25, inner_folds = 5) {
  # Outer folds assess the complete procedure; inner folds select gamma.
  sample_sizes <- validate_cv_data(data, outer_folds)
  outer_e <- split_idx(outer_folds, sample_sizes[["environment_e"]])
  outer_o <- split_idx(outer_folds, sample_sizes[["environment_o"]])

  fold_results <- lapply(seq_len(outer_folds), function(fold) {
    outer_training <- subset_data(data, outer_e != fold, outer_o != fold)
    outer_validation <- subset_data(data, outer_e == fold, outer_o == fold)

    inner_cv <- cross_validation(inner_folds, outer_training, estimators)
    selected_gamma <- as.numeric(inner_cv$best_gamma)
    training_moments <- moments(outer_training)
    estimates <- list(
      OLS = drop(compute_ols(training_moments)),
      CR = drop(make_estimator(selected_gamma)(outer_training)),
      CD = drop(compute_cd(training_moments))
    )

    do.call(rbind, lapply(names(estimates), function(method) {
      beta <- estimates[[method]]
      r_delta <- calc_Rdelta(outer_validation, beta)
      r_plus <- in_sample_risk(outer_validation, beta)
      data.frame(
        outer_fold = fold,
        method = method,
        selected_gamma = selected_gamma,
        Rdelta_plus = r_delta,
        Rplus = r_plus,
        Rgamma = 0.5 * r_plus + 0.5 * selected_gamma * r_delta
      )
    }))
  })

  list(
    fold_results = do.call(rbind, fold_results),
    fold_assignment = list(environment_e = outer_e, environment_o = outer_o)
  )
}

summarise_nested_cv <- function(results) {
  # Average each outer validation criterion across the 25 held-out folds.
  methods <- unique(results$method)
  do.call(rbind, lapply(methods, function(method) {
    current <- results[results$method == method, , drop = FALSE]
    data.frame(
      method = method,
      mean_Rdelta_plus = mean(current$Rdelta_plus),
      mean_Rplus = mean(current$Rplus),
      mean_Rgamma = mean(current$Rgamma),
      sd_Rgamma = sd(current$Rgamma)
    )
  }))
}

plot_test_risks <- function(
    risks, gamma, file = NULL, second_label = "CR", title = NULL) {
  # Save the distribution of unseen-environment risks.
  all_risks <- c(risks$OLS, risks$CR)
  risk_range <- range(all_risks, finite = TRUE)
  y_padding <- max(0.02, 0.05 * diff(risk_range))
  y_limits <- c(max(0, risk_range[1] - y_padding), risk_range[2] + y_padding)

  if (!is.null(file)) {
    grDevices::pdf(file, width = 6, height = 4.5, useDingbats = FALSE)
    on.exit(grDevices::dev.off())
  } else {
    if (grDevices::dev.cur() == 1L) {
      grDevices::dev.new(width = 6, height = 4.5)
    }
    layout(1)
    par(mfrow = c(1, 1))
  }
  method_risks <- list(risks$OLS, risks$CR)
  names(method_risks) <- c("OLS", second_label)
  boxplot(
    method_risks,
    ylab = expression(R[e](beta)),
    ylim = y_limits,
    main = title,
    outline = FALSE
  )
  stripchart(
    method_risks, vertical = TRUE, method = "jitter",
    add = TRUE, pch = 16, cex = 0.6, col = grDevices::adjustcolor("black", 0.5)
  )
  invisible(file)
}

plot_cv_method_comparison <- function(risks, file = NULL) {
  # Compare test risks after selecting gamma with the two CV criteria.
  method_risks <- list(
    OLS = risks$OLS,
    CR = risks$CR,
    CR_Kania = risks$CR_Kania
  )
  all_risks <- unlist(method_risks, use.names = FALSE)
  risk_range <- range(all_risks, finite = TRUE)
  y_padding <- max(0.02, 0.05 * diff(risk_range))
  y_limits <- c(max(0, risk_range[1] - y_padding), risk_range[2] + y_padding)

  if (!is.null(file)) {
    grDevices::pdf(file, width = 7, height = 4.5, useDingbats = FALSE)
    on.exit(grDevices::dev.off())
  } else {
    if (grDevices::dev.cur() == 1L) {
      grDevices::dev.new(width = 7, height = 4.5)
    }
    layout(1)
    par(mfrow = c(1, 1))
  }
  boxplot(
    method_risks,
    names = c("OLS", "CR", "CR (Kania CV)"),
    ylab = expression(R[e](beta)),
    ylim = y_limits, outline = FALSE
  )
  stripchart(
    method_risks, vertical = TRUE, method = "jitter",
    add = TRUE, pch = 16, cex = 0.5,
    col = grDevices::adjustcolor("black", 0.4)
  )
  invisible(file)
}

compute_training_risk_path <- function(training_data, gamma_values, gamma_star) {
  # Follow beta_gamma from OLS towards the CD limit.
  beta_path <- vapply(
    gamma_values,
    function(gamma) drop(make_estimator(gamma)(training_data)),
    numeric(ncol(training_data$Xe))
  )
  m <- build_matrices(training_data)
  beta_cd <- drop(compute_cd(m))
  beta_path <- cbind(beta_path, CD = beta_cd)

  r_delta <- apply(beta_path, 2L, function(beta) {
    difference <- beta - beta_cd
    drop(t(difference) %*% m$Gdelta_plus %*% difference)
  })
  r_plus <- apply(beta_path, 2L, function(beta) {
    in_sample_risk(training_data, beta)
  })
  r_gamma <- 0.5 * r_plus + 0.5 * gamma_star * r_delta

  max_gamma <- max(gamma_values)
  data.frame(
    gamma = c(gamma_values, Inf),
    path_position = c(log1p(gamma_values) / log1p(max_gamma), 1),
    Rdelta_plus = r_delta,
    Rplus = r_plus,
    Rgamma = r_gamma
  )
}

plot_training_risk_path <- function(path, gamma_star, file) {
  # The transformed x-axis includes both gamma = 0 and the CD limit.
  selected <- which.min(abs(path$gamma - gamma_star))
  gamma_position <- path$path_position[selected]
  y_range <- range(path$Rgamma, finite = TRUE)
  padding <- max(0.02, 0.08 * diff(y_range))

  grDevices::pdf(file, width = 7, height = 5, family = "serif", useDingbats = FALSE)
  on.exit(grDevices::dev.off())
  old_par <- graphics::par(
    family = "serif", mar = c(5.2, 5.2, 1.2, 1.2),
    mgp = c(3.1, 0.8, 0), lend = "round", ljoin = "round"
  )
  on.exit(graphics::par(old_par), add = TRUE)

  plot(
    path$path_position, path$Rgamma,
    type = "l", lwd = 1.5, axes = FALSE,
    xlab = expression(hat(beta)[gamma]),
    ylab = expression(hat(R)[gamma[CV]](hat(beta)[gamma])),
    ylim = y_range + c(-padding, padding), xaxs = "i"
  )
  axis(
    1, at = c(0, gamma_position, 1),
    labels = expression(hat(beta)[OLS], hat(beta)[CR], hat(beta)[RL])
  )
  axis(2, las = 1)
  box(bty = "l")
  abline(v = gamma_position, lty = "dotdash", col = "grey35")
  points(
    path$path_position[selected], path$Rgamma[selected],
    pch = 16, cex = 0.9
  )
  invisible(file)
}
