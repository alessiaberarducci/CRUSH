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
environment_l21 <- read.csv(
  file.path(data_folder, "uniform_l_21_mid.csv")
)
environment_l22 <- read.csv(
  file.path(data_folder, "uniform_l_22_mid.csv")
)

# Retain every varying numerical covariate and remove record metadata.
metadata <- c("timestamp", "config", "counter", "flag", "intervention")
common_columns <- intersect(names(environment_l21), names(environment_l22))
candidate_predictors <- setdiff(common_columns, c(metadata, "ir_2"))
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
gamma_grid <- seq(0, 30, by = 0.1)
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

# Save the cross-validation and coefficient-path visualizations.
dir.create("visualization", showWarnings = FALSE)
show_plots <- interactive()

plot_cv_risk(
  cv_result, gamma_grid,
  file = "visualization/chamber_risk.pdf", x_max = 30
)
if (show_plots) {
  plot_cv_risk(cv_result, gamma_grid, x_max = 30)
}

plot_beta_path(
  beta_gamma, gamma_grid, beta_cd, best_gamma,
  file = "visualization/chamber.pdf", x_max = 30,
  legend_ncol = 7, legend_cex = 0.5,
  pdf_width = 11, pdf_height = 7, bottom_margin = 12,
  y_limits = c(-700, 700), coefficient_names = predictors
)
if (show_plots) {
  plot_beta_path(
    beta_gamma, gamma_grid, beta_cd, best_gamma,
    x_max = 30, legend_ncol = 7, legend_cex = 0.5,
    bottom_margin = 12, y_limits = c(-700, 700),
    coefficient_names = predictors
  )
}

# Display the paths of the five ground-truth causal predictors separately.
direct_causes <- c("red", "green", "blue", "l_21", "l_22")
direct_cause_index <- match(direct_causes, predictors)
plot_beta_path(
  beta_gamma[direct_cause_index, , drop = FALSE],
  gamma_grid, beta_cd[direct_cause_index], best_gamma,
  file = "visualization/chamber_direct_causes.pdf", x_max = 30,
  legend_ncol = 5, legend_cex = 0.7,
  y_limits = c(-700, 700), coefficient_names = direct_causes
)
if (show_plots) {
  plot_beta_path(
    beta_gamma[direct_cause_index, , drop = FALSE],
    gamma_grid, beta_cd[direct_cause_index], best_gamma,
    x_max = 30, legend_ncol = 5, legend_cex = 0.7,
    y_limits = c(-700, 700), coefficient_names = direct_causes
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
cat("Environments: l_21_mid and l_22_mid\n")
cat("Sample sizes:", nrow(data$Xe), "and", nrow(data$Xo), "\n")
cat("Pooled response center:", response_center, "\n")
cat("Number of varying covariates:", length(predictors), "\n")
cat(
  "Constant covariates excluded:",
  paste(constant_predictors, collapse = ", "), "\n"
)
cat("Best gamma:", best_gamma, "\n")
print(coefficient_table, row.names = FALSE)
cat("Direct-cause standardized CR effects (beta_CR / sd(Y)):\n")
print(
  direct_cause_graph[, c("cause", "beta_hat", "standardized_effect")],
  row.names = FALSE
)
