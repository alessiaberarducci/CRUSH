# Causal Chamber application

# Load estimation, cross-validation and plotting functions.
source("functions/datautil.R")
source("functions/cd.R")
source("functions/ols.R")
source("functions/measurements.R")
source("functions/cross_valid_fun_n.R")
source("functions/visualizations.R")

# Read interventions on two distinct direct causes of the target sensor.
data_folder <- "lt_interventions_standard_v1"
environment_file_1 <- Sys.getenv(
  "CHAMBER_ENVIRONMENT_1", "uniform_l_21_mid.csv"
)
environment_file_2 <- Sys.getenv(
  "CHAMBER_ENVIRONMENT_2", "uniform_l_22_mid.csv"
)
chamber_output_suffix <- Sys.getenv("CHAMBER_OUTPUT_SUFFIX", "_1")
direct_path_y_limits <- c(
  as.numeric(Sys.getenv("CHAMBER_DIRECT_Y_MIN", "-600")),
  as.numeric(Sys.getenv("CHAMBER_DIRECT_Y_MAX", "5000"))
)
chamber_output_file <- function(directory, stem, extension) {
  file.path(
    directory,
    paste0(stem, chamber_output_suffix, ".", extension)
  )
}
environment_l21 <- read.csv(
  file.path(data_folder, environment_file_1)
)
environment_l22 <- read.csv(
  file.path(data_folder, environment_file_2)
)

# Retain every varying numerical covariate and remove record metadata.
metadata <- c("timestamp", "config", "counter", "flag", "intervention")
manually_excluded_predictors <- trimws(strsplit(
  Sys.getenv("CHAMBER_EXCLUDE_PREDICTORS", ""), ",", fixed = TRUE
)[[1]])
manually_excluded_predictors <- manually_excluded_predictors[
  nzchar(manually_excluded_predictors)
]
common_columns <- intersect(names(environment_l21), names(environment_l22))
candidate_predictors <- setdiff(
  common_columns,
  c(metadata, "ir_2", manually_excluded_predictors)
)
is_numeric <- vapply(
  candidate_predictors,
  function(variable) {
    is.numeric(environment_l21[[variable]]) &&
      is.numeric(environment_l22[[variable]])
  },
  logical(1)
)
candidate_predictors <- candidate_predictors[is_numeric]

pooled_covariates <- rbind(
  environment_l21[, candidate_predictors, drop = FALSE],
  environment_l22[, candidate_predictors, drop = FALSE]
)
is_varying <- vapply(
  pooled_covariates,
  function(values) length(unique(values[is.finite(values)])) > 1,
  logical(1)
)
predictors <- candidate_predictors[is_varying]
constant_predictors <- candidate_predictors[!is_varying]
required_columns <- c(predictors, "ir_2")

# Standardization is needed because the raw moment matrix is numerically singular.
predictor_center <- colMeans(pooled_covariates[, predictors, drop = FALSE])
predictor_scale <- vapply(
  pooled_covariates[, predictors, drop = FALSE],
  sd, numeric(1)
)
response_center <- mean(c(environment_l21$ir_2, environment_l22$ir_2))

# Use the continuous infrared intensity as the response.
make_environment <- function(chamber_data) {
  complete <- complete.cases(chamber_data[, required_columns])
  X <- scale(
    as.matrix(chamber_data[complete, predictors, drop = FALSE]),
    center = predictor_center,
    scale = predictor_scale
  )
  list(
    X = X,
    y = chamber_data$ir_2[complete] - response_center
  )
}

env_l21 <- make_environment(environment_l21)
env_l22 <- make_environment(environment_l22)
data <- list(
  Xe = env_l21$X,
  ye = env_l21$y,
  Xo = env_l22$X,
  yo = env_l22$y
)

# Select gamma using five-fold cross-validation.
set.seed(123)
gamma_max <- 10
gamma_grid <- seq(0, gamma_max, by = 0.1)
beta_path_gamma_max <- 5
estimators <- setNames(lapply(gamma_grid, make_estimator), gamma_grid)
cv_result <- cross_validation(5, data, estimators)
best_gamma <- as.numeric(cv_result$best_gamma)

# Refit OLS, CR and CD using all observations from both environments.
m <- build_matrices(data)
beta_ols <- drop(compute_ols(m))
beta_cd <- drop(compute_cd(m))
beta_gamma <- vapply(
  gamma_grid,
  function(gamma) drop(make_estimator(gamma)(data)),
  numeric(m$p)
)
best_gamma_index <- which.min(abs(gamma_grid - best_gamma))
beta_cr <- beta_gamma[, best_gamma_index]
direct_causes <- c("red", "green", "blue", "l_21", "l_22")

# Estimate the conditional HC0 covariance matrix of beta_OLS.
# The two environments receive equal weight in Gplus and Zplus.
G_plus_inverse <- solve(m$Gplus)
residual_e <- drop(data$ye - data$Xe %*% beta_ols)
residual_o <- drop(data$yo - data$Xo %*% beta_ols)
score_e <- data$Xe * residual_e
score_o <- data$Xo * residual_o
variance_Ze <- crossprod(score_e) / nrow(data$Xe)^2
variance_Zo <- crossprod(score_o) / nrow(data$Xo)^2
variance_Zplus <- variance_Ze + variance_Zo
variance_beta_ols <- G_plus_inverse %*% variance_Zplus %*%
  t(G_plus_inverse)
variance_beta_ols <- 0.5 * (variance_beta_ols + t(variance_beta_ols))
dimnames(variance_beta_ols) <- list(
  paste0("beta_", predictors), paste0("beta_", predictors)
)

# Restore the conditional covariance matrix of beta_CV at gamma_CV.
identity_p <- diag(m$p)
G_delta_inverse <- qr.solve(m$Gdelta, identity_p)
A_gamma <- best_gamma * m$Gdelta_plus %*% G_delta_inverse
G_gamma <- m$Gplus + best_gamma * m$Gdelta_plus
G_gamma_inverse <- solve(G_gamma)
residual_e_cv <- drop(data$ye - data$Xe %*% beta_cr)
residual_o_cv <- drop(data$yo - data$Xo %*% beta_cr)
score_e_cv <- data$Xe * residual_e_cv
score_o_cv <- data$Xo * residual_o_cv
variance_Ze_cv <- crossprod(score_e_cv) / nrow(data$Xe)^2
variance_Zo_cv <- crossprod(score_o_cv) / nrow(data$Xo)^2
variance_Zgamma <-
  (identity_p + A_gamma) %*% variance_Ze_cv %*% t(identity_p + A_gamma) +
  (identity_p - A_gamma) %*% variance_Zo_cv %*% t(identity_p - A_gamma)
variance_beta_cv <- G_gamma_inverse %*% variance_Zgamma %*%
  t(G_gamma_inverse)
variance_beta_cv <- 0.5 * (variance_beta_cv + t(variance_beta_cv))
dimnames(variance_beta_cv) <- list(
  paste0("beta_", predictors), paste0("beta_", predictors)
)

build_conditional_variance_table <- function(beta_hat, variance_matrix) {
  coefficient_variance <- pmax(diag(variance_matrix), 0)
  coefficient_se <- sqrt(coefficient_variance)
  result <- data.frame(
    parameter = paste0("beta_", predictors),
    variable = predictors,
    beta_hat = beta_hat,
    variance = coefficient_variance,
    standard_error = coefficient_se,
    beta_over_standard_error = beta_hat / coefficient_se,
    row.names = NULL
  )
  result$is_direct_cause <- result$variable %in% direct_causes
  result$absolute_beta_over_standard_error <- abs(
    result$beta_over_standard_error
  )
  result[
    order(-result$absolute_beta_over_standard_error), , drop = FALSE
  ]
}

conditional_variance_table_cv <- build_conditional_variance_table(
  beta_cr, variance_beta_cv
)
conditional_variance_table_ols <- build_conditional_variance_table(
  beta_ols, variance_beta_ols
)

# Save the cross-validation and coefficient-path visualizations.
chamber_output_directory <- "visualization/chamber_output"
dir.create(chamber_output_directory, recursive = TRUE, showWarnings = FALSE)
show_plots <- interactive()
write.csv(
  conditional_variance_table_cv,
  chamber_output_file(
    chamber_output_directory, "chamber_conditional_variance", "csv"
  ),
  row.names = FALSE
)
write.csv(
  conditional_variance_table_ols,
  chamber_output_file(
    chamber_output_directory, "chamber_OLS_conditional_variance", "csv"
  ),
  row.names = FALSE
)
write.csv(
  data.frame(
    parameter = rownames(variance_beta_ols),
    variance_beta_ols,
    row.names = NULL,
    check.names = FALSE
  ),
  chamber_output_file(
    chamber_output_directory, "chamber_OLS_variance_covariance", "csv"
  ),
  row.names = FALSE
)
write.csv(
  data.frame(
    parameter = rownames(variance_beta_cv),
    variance_beta_cv,
    row.names = NULL,
    check.names = FALSE
  ),
  chamber_output_file(
    chamber_output_directory, "chamber_CV_variance_covariance", "csv"
  ),
  row.names = FALSE
)

# Correlations of the response and covariates pooled across the two training
# environments used in the current Chamber analysis.
pooled_model_variables <- rbind(
  cbind(ir_2 = data$ye, data$Xe),
  cbind(ir_2 = data$yo, data$Xo)
)
colnames(pooled_model_variables) <- c("ir_2", predictors)
pooled_correlation <- stats::cor(
  pooled_model_variables, use = "pairwise.complete.obs"
)
plot_correlation_matrix(
  pooled_correlation,
  file = chamber_output_file(
    chamber_output_directory, "chamber_correlation", "pdf"
  ),
  variable_names = colnames(pooled_correlation)
)

# Create readable OLS and CV tables and highlight direct causes in red.
write_conditional_variance_pdf <- function(
    variance_table, file, estimator_label) {
  grDevices::pdf(
    file, width = 11, height = 8.5,
    family = "serif", useDingbats = FALSE
  )
  on.exit(grDevices::dev.off())
  graphics::par(mar = c(0.5, 0.5, 1.8, 0.5), family = "serif")
  graphics::plot.new()
  n_table_rows <- nrow(variance_table)
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, n_table_rows + 2.2))
  column_x <- c(0.02, 0.34, 0.50, 0.64, 0.80, 0.96)
  column_alignment <- c(0, 1, 1, 1, 1, 1)
  column_headers <- c(
    "Parameter", "Beta hat", "Variance", "Std. error",
    "Beta / std. error", "Absolute value"
  )
  header_y <- n_table_rows + 1.25
  for (j in seq_along(column_x)) {
    graphics::text(
      column_x[j], header_y, column_headers[j],
      adj = c(column_alignment[j], 0.5), font = 2, cex = 0.82
    )
  }
  graphics::segments(0.02, header_y - 0.48, 0.98, header_y - 0.48)
  for (i in seq_len(n_table_rows)) {
    current <- variance_table[i, ]
    row_y <- n_table_rows - i + 0.55
    row_colour <- if (current$is_direct_cause) "red3" else "black"
    row_values <- c(
      current$parameter,
      formatC(current$beta_hat, format = "f", digits = 3),
      formatC(current$variance, format = "f", digits = 3),
      formatC(current$standard_error, format = "f", digits = 3),
      formatC(current$beta_over_standard_error, format = "f", digits = 3),
      formatC(
        current$absolute_beta_over_standard_error,
        format = "f", digits = 3
      )
    )
    for (j in seq_along(column_x)) {
      graphics::text(
        column_x[j], row_y, row_values[j],
        adj = c(column_alignment[j], 0.5),
        col = row_colour, cex = 0.75
      )
    }
  }
  graphics::mtext(
    paste("Conditional variance of the Chamber", estimator_label, "estimator"),
    side = 3, line = 0.45, font = 2, cex = 1.05
  )
  graphics::mtext(
    paste(
      "Rows ordered by absolute", estimator_label,
      "beta hat / standard error; direct causes in red"
    ),
    side = 3, line = -0.65, cex = 0.75
  )
}

write_conditional_variance_pdf(
  conditional_variance_table_cv,
  chamber_output_file(
    chamber_output_directory,
    "chamber_conditional_variance_table_absolute_standardized", "pdf"
  ),
  "CV"
)
write_conditional_variance_pdf(
  conditional_variance_table_ols,
  chamber_output_file(
    chamber_output_directory,
    "chamber_OLS_conditional_variance_table_absolute_standardized", "pdf"
  ),
  "OLS"
)

plot_cv_risk(
  cv_result, gamma_grid,
  file = chamber_output_file(chamber_output_directory, "chamber_risk", "pdf"),
  x_max = gamma_max
)
if (show_plots) {
  plot_cv_risk(cv_result, gamma_grid, x_max = gamma_max)
}

# Display coefficient paths only up to gamma = 5.
gamma_grid_wide <- seq(0, beta_path_gamma_max, length.out = 201)
beta_gamma_wide <- vapply(
  gamma_grid_wide,
  function(gamma) drop(make_estimator(gamma)(data)),
  numeric(m$p)
)
full_path_range <- range(beta_gamma_wide, finite = TRUE)
full_path_padding <- 0.08 * diff(full_path_range)
if (!is.finite(full_path_padding) || full_path_padding == 0) {
  full_path_padding <- max(1, 0.08 * max(abs(full_path_range)))
}
full_path_y_limits <- full_path_range + c(-1, 1) * full_path_padding

plot_beta_path(
  beta_gamma_wide, gamma_grid_wide, beta_cd, best_gamma,
  file = chamber_output_file(chamber_output_directory, "chamber", "pdf"),
  x_max = beta_path_gamma_max,
  legend_ncol = 7, legend_cex = 0.5,
  pdf_width = 11, pdf_height = 7, bottom_margin = 12,
  y_limits = full_path_y_limits, coefficient_names = predictors,
  log_x = FALSE
)
if (show_plots) {
  plot_beta_path(
    beta_gamma_wide, gamma_grid_wide, beta_cd, best_gamma,
    x_max = beta_path_gamma_max, legend_ncol = 7, legend_cex = 0.5,
    bottom_margin = 12, y_limits = full_path_y_limits,
    coefficient_names = predictors, log_x = FALSE
  )
}

# Display the paths of the five ground-truth causal predictors separately.
direct_cause_index <- match(direct_causes, predictors)
plot_beta_path(
  beta_gamma_wide[direct_cause_index, , drop = FALSE],
  gamma_grid_wide, beta_cd[direct_cause_index], best_gamma,
  file = chamber_output_file(
    chamber_output_directory, "chamber_direct_causes", "pdf"
  ),
  x_max = beta_path_gamma_max,
  legend_ncol = 5, legend_cex = 0.7,
  y_limits = direct_path_y_limits, coefficient_names = direct_causes,
  log_x = FALSE
)
if (show_plots) {
  plot_beta_path(
    beta_gamma_wide[direct_cause_index, , drop = FALSE],
    gamma_grid_wide, beta_cd[direct_cause_index], best_gamma,
    x_max = beta_path_gamma_max, legend_ncol = 5, legend_cex = 0.7,
    y_limits = direct_path_y_limits, coefficient_names = direct_causes,
    log_x = FALSE
  )
}

# Report the standardized CR effects used in the TikZ causal-effect graph.
response_sd <- sd(c(environment_l21$ir_2, environment_l22$ir_2))
direct_cause_graph <- data.frame(
  cause = direct_causes,
  beta_hat = beta_cr[direct_cause_index],
  standardized_effect = beta_cr[direct_cause_index] / response_sd,
  row.names = NULL
)

# Report the selected value and the three coefficient estimates.
coefficient_table <- data.frame(
  variable = predictors,
  OLS = beta_ols,
  CR = beta_cr,
  CD = beta_cd,
  row.names = NULL
)
cat("Environments:", environment_file_1, "and", environment_file_2, "\n")
cat("Sample sizes:", nrow(data$Xe), "and", nrow(data$Xo), "\n")
cat("Pooled response center:", response_center, "\n")
cat("Number of varying covariates:", length(predictors), "\n")
cat(
  "Constant covariates excluded:",
  paste(constant_predictors, collapse = ", "), "\n"
)
cat(
  "Manually excluded covariates:",
  paste(manually_excluded_predictors, collapse = ", "), "\n"
)
cat("Best gamma:", best_gamma, "\n")
print(coefficient_table, row.names = FALSE)
cat("Conditional variance of beta_CV at gamma_CV:\n")
print(conditional_variance_table_cv, row.names = FALSE)
cat("Conditional variance of beta_OLS:\n")
print(conditional_variance_table_ols, row.names = FALSE)
cat("Direct-cause standardized CR effects (beta_CR / sd(Y)):\n")
print(
  direct_cause_graph[, c("cause", "beta_hat", "standardized_effect")],
  row.names = FALSE
)
