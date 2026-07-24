# RPE comparison with and without a non-penalized intercept

source("functions/datautil.R")
source("functions/cd.R")
source("functions/ols.R")
source("functions/measurements.R")
source("functions/cross_valid_fun_n.R")
source("functions/kania_function.R")
source("functions/rpe_robust_functions.R")

environment_column <- "interventions"
observed_genes <- c(
  "ENSG00000187514", "ENSG00000075624", "ENSG00000147604",
  "ENSG00000110700", "ENSG00000172757", "ENSG00000133112",
  "ENSG00000067225", "ENSG00000108518", "ENSG00000125691",
  "ENSG00000173812"
)
target_gene <- "ENSG00000173812"
predictor_genes <- setdiff(observed_genes, target_gene)
first_training <- c(
  observed_genes[1:5], "ENSG00000125484", "ENSG00000101057"
)
second_training <- c(
  observed_genes[6:10], "ENSG00000037637", "ENSG00000215712"
)
training_interventions <- c(first_training, second_training)

dat <- read.csv("dataset_rpe.csv", check.names = FALSE)
dat <- filter_environments(dat, environment_column, minimum_cells = 100)
training <- make_two_environment_data(
  dat, environment_column, first_training, second_training,
  predictor_genes, target_gene
)
raw_data <- training$data
centering <- compute_common_centering(raw_data)
centered_data <- centering$data

gamma_grid <- seq(0, 40, by = 0.01)
estimators <- setNames(lapply(gamma_grid, make_estimator), gamma_grid)

fit_raw_specification <- function(data, seed) {
  set.seed(seed)
  our_cv <- cross_validation(5, data, estimators)
  kania_cv <- kania_cross_validation(
    data, estimators, our_cv$fold_assignment
  )
  gamma_ours <- as.numeric(our_cv$best_gamma)
  gamma_kania <- as.numeric(kania_cv$best_gamma)
  m <- moments(data)
  list(
    gamma_ours = gamma_ours,
    gamma_kania = gamma_kania,
    fold_assignment = our_cv$fold_assignment,
    beta = list(
      OLS = drop(compute_ols(m)),
      CR = drop(make_estimator(gamma_ours)(data)),
      CR_Kania = drop(make_estimator(gamma_kania)(data))
    )
  )
}

# Resetting the seed gives both specifications identical row assignments.
raw_fit <- fit_raw_specification(raw_data, seed = 2029)
fold_assignment <- raw_fit$fold_assignment
centered_our_cv <- centered_cross_validation(
  raw_data, estimators, fold_assignment, validation_loss = "ours"
)
centered_kania_cv <- centered_cross_validation(
  raw_data, estimators, fold_assignment, validation_loss = "kania"
)
centered_gamma_ours <- as.numeric(centered_our_cv$best_gamma)
centered_gamma_kania <- as.numeric(centered_kania_cv$best_gamma)
centered_moments <- moments(centered_data)
centered_fit <- list(
  gamma_ours = centered_gamma_ours,
  gamma_kania = centered_gamma_kania,
  beta = list(
    OLS = drop(compute_ols(centered_moments)),
    CR = drop(make_estimator(centered_gamma_ours)(centered_data)),
    CR_Kania = drop(make_estimator(centered_gamma_kania)(centered_data))
  )
)
raw_fit$beta <- lapply(raw_fit$beta, setNames, predictor_genes)
centered_fit$beta <- lapply(centered_fit$beta, setNames, predictor_genes)

raw_intercepts <- c(OLS = 0, CR = 0, CR_Kania = 0)
centered_intercepts <- vapply(
  centered_fit$beta,
  recover_intercept,
  numeric(1),
  centering = centering
)

test_environments <- setdiff(
  unique(as.character(dat[[environment_column]])),
  c(training_interventions, "excluded")
)

evaluate_specification <- function(fit, intercepts, specification) {
  ours <- compute_test_risks(
    dat, test_environments, environment_column,
    fit$beta$OLS, fit$beta$CR, predictor_genes, target_gene,
    intercept_ols = intercepts[["OLS"]],
    intercept_cr = intercepts[["CR"]]
  )
  kania <- compute_test_risks(
    dat, test_environments, environment_column,
    fit$beta$OLS, fit$beta$CR_Kania, predictor_genes, target_gene,
    intercept_ols = intercepts[["OLS"]],
    intercept_cr = intercepts[["CR_Kania"]]
  )
  comparison <- merge(
    ours,
    setNames(kania[, c("environment", "CR")], c("environment", "CR_Kania")),
    by = "environment", all = FALSE, sort = FALSE
  )
  rbind(
    data.frame(specification, environment = comparison$environment,
               method = "OLS", risk = comparison$OLS),
    data.frame(specification, environment = comparison$environment,
               method = "CR", risk = comparison$CR),
    data.frame(specification, environment = comparison$environment,
               method = "CR_Kania", risk = comparison$CR_Kania)
  )
}

comparison_risks <- rbind(
  evaluate_specification(raw_fit, raw_intercepts, "No intercept"),
  evaluate_specification(
    centered_fit, centered_intercepts, "Common centering"
  )
)

gamma_lookup <- c(
  "No intercept.OLS" = 0,
  "No intercept.CR" = raw_fit$gamma_ours,
  "No intercept.CR_Kania" = raw_fit$gamma_kania,
  "Common centering.OLS" = 0,
  "Common centering.CR" = centered_fit$gamma_ours,
  "Common centering.CR_Kania" = centered_fit$gamma_kania
)
groups <- split(
  comparison_risks,
  interaction(
    comparison_risks$specification, comparison_risks$method,
    drop = TRUE, sep = "."
  )
)
comparison_summary <- do.call(rbind, lapply(names(groups), function(group) {
  current <- groups[[group]]
  data.frame(
    specification = current$specification[1],
    method = current$method[1],
    selected_gamma = gamma_lookup[[group]],
    mean_risk = mean(current$risk),
    median_risk = median(current$risk),
    sd_risk = sd(current$risk),
    worst_risk = max(current$risk)
  )
}))

dir.create("visualization", showWarnings = FALSE)
write.csv(
  comparison_summary,
  "visualization/rpe_centering_comparison_summary.csv",
  row.names = FALSE
)

plot_centering_comparison <- function(file = NULL) {
  if (!is.null(file)) {
    grDevices::pdf(file, width = 9, height = 4.8, useDingbats = FALSE)
    on.exit(grDevices::dev.off())
  } else if (grDevices::dev.cur() == 1L) {
    grDevices::dev.new(width = 9, height = 4.8)
  }
  old_par <- par(mfrow = c(1, 2), mar = c(4.4, 4.4, 2.2, 0.8))
  on.exit(par(old_par), add = TRUE)
  specifications <- c("No intercept", "Common centering")
  limits <- range(unlist(lapply(
    split(
      comparison_risks$risk,
      interaction(comparison_risks$specification, comparison_risks$method)
    ),
    function(values) boxplot.stats(values)$stats
  )))
  padding <- 0.04 * diff(limits)
  for (specification in specifications) {
    current <- comparison_risks[
      comparison_risks$specification == specification, ,
      drop = FALSE
    ]
    values <- split(
      current$risk,
      factor(current$method, levels = c("OLS", "CR", "CR_Kania"))
    )
    boxplot(
      values, names = c("OLS", "CR", "CR\n(Kania CV)"),
      outline = FALSE, ylim = limits + c(-padding, padding),
      ylab = expression(R[e](beta)), main = specification
    )
  }
}

plot_centering_comparison(
  "visualization/rpe_centering_comparison_boxplots.pdf"
)
if (interactive()) {
  plot_centering_comparison()
}

print(comparison_summary, row.names = FALSE)
print(centered_intercepts)
