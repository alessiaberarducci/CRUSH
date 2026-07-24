# Paper Application 1. Robust analysis of the RPE dataset

# Load the estimators, risk measures and RPE helper functions.
source("functions/datautil.R")
source("functions/cd.R")
source("functions/ols.R")
source("functions/measurements.R")
source("functions/cross_valid_fun_n.R")
source("functions/kania_function.R")
source("functions/rpe_robust_functions.R")

set.seed(123)

# Choose the target and predictors.
environment_column <- "interventions"
observed_genes <- c(
  "ENSG00000187514", "ENSG00000075624", "ENSG00000147604",
  "ENSG00000110700", "ENSG00000172757", "ENSG00000133112",
  "ENSG00000067225", "ENSG00000108518", "ENSG00000125691",
  "ENSG00000173812"
)
target_gene <- "ENSG00000173812"
predictor_genes <- setdiff(observed_genes, target_gene)

# Pool seven distinct interventions into each training environment.
# The additional groups have similar sizes and keep the totals balanced.
first_training_interventions <- c(
  observed_genes[1:5],
  "ENSG00000125484", "ENSG00000101057"
)
second_training_interventions <- c(
  observed_genes[6:10],
  "ENSG00000037637", "ENSG00000215712"
)
training_interventions <- c(
  first_training_interventions, second_training_interventions
)

# Load the data and discard small intervention groups.
dat <- read.csv("dataset_rpe.csv", check.names = FALSE)
dat <- filter_environments(dat, environment_column, minimum_cells = 100)
environment_counts <- sort(table(dat[[environment_column]]), decreasing = TRUE)

# Build the two-environment object used for training.
training <- make_two_environment_data(
  dat, environment_column,
  first_training_interventions, second_training_interventions,
  predictor_genes, target_gene
)
data_train <- training$data

# Select gamma with five-fold cross-validation.
gamma_grid <- seq(0, 15, by = 0.01)
estimators <- setNames(lapply(gamma_grid, make_estimator), gamma_grid)
cv_results <- cross_validation(5, data_train, estimators)
gamma_star <- as.numeric(cv_results$best_gamma)
kania_cv_results <- kania_cross_validation(
  data_train, estimators, cv_results$fold_assignment
)
gamma_kania <- as.numeric(kania_cv_results$best_gamma)

# Assess the full selection procedure with 25 outer validation folds.
nested_cv <- nested_cross_validation(
  data_train, estimators, outer_folds = 25, inner_folds = 5
)
nested_cv_summary <- summarise_nested_cv(nested_cv$fold_results)
nested_cr <- nested_cv$fold_results[
  nested_cv$fold_results$method == "CR", , drop = FALSE
]
print(nested_cv_summary, row.names = FALSE)
cat(
  "Nested-CV gamma: median", median(nested_cr$selected_gamma),
  "| range", range(nested_cr$selected_gamma), "\n"
)

# Refit OLS, CR and CD using all training cells.
m_train <- moments(data_train)
beta_ols <- drop(compute_ols(m_train))
beta_cr <- drop(make_estimator(gamma_star)(data_train))
beta_cr_kania <- drop(make_estimator(gamma_kania)(data_train))
beta_cd <- drop(compute_cd(m_train))
names(beta_ols) <- names(beta_cr) <- names(beta_cr_kania) <-
  names(beta_cd) <- predictor_genes

coefficient_table <- data.frame(
  gene = predictor_genes, OLS = beta_ols, CR = beta_cr,
  CR_Kania = beta_cr_kania
)
print(coefficient_table, row.names = FALSE)

# Compare the three training objectives at the selected gamma.
training_risk_table <- compare_training_risks(
  data_train,
  list(
    OLS = beta_ols, CR = beta_cr, CR_Kania = beta_cr_kania, CD = beta_cd
  ),
  gamma_star
)
print(training_risk_table, row.names = FALSE)

# Test only on intervention datasets not pooled into the training environments.
all_environments <- unique(as.character(dat[[environment_column]]))
test_environments <- setdiff(
  all_environments,
  c(training_interventions, "excluded")
)
risk_by_environment <- compute_test_risks(
  dat, test_environments, environment_column, beta_ols, beta_cr,
  predictor_genes, target_gene
)
risk_by_environment_kania <- compute_test_risks(
  dat, test_environments, environment_column, beta_ols, beta_cr_kania,
  predictor_genes, target_gene
)
kania_risk <- setNames(
  risk_by_environment_kania[, c("environment", "CR")],
  c("environment", "CR_Kania")
)
risk_comparison <- merge(
  risk_by_environment, kania_risk,
  by = "environment", all = FALSE, sort = FALSE
)
# Summarise marginal and paired differences between OLS and CR.
risk_summary <- summarise_test_risks(risk_by_environment)
paired_comparison <- data.frame(
  mean_CR_minus_OLS = mean(risk_by_environment$CR - risk_by_environment$OLS),
  median_CR_minus_OLS = median(risk_by_environment$CR - risk_by_environment$OLS),
  environments_where_CR_is_better = sum(risk_by_environment$CR < risk_by_environment$OLS),
  total_environments = nrow(risk_by_environment),
  proportion_where_CR_is_better = mean(risk_by_environment$CR < risk_by_environment$OLS)
)

print(risk_summary)
print(paired_comparison)
# Save the risk distribution and print the main sample information.
dir.create("visualization", showWarnings = FALSE)
plot_test_risks(
  risk_by_environment, gamma_star,
  "visualization/rpe1_unseen_environment_risk_boxplot.pdf"
)
plot_cv_method_comparison(
  risk_comparison,
  "visualization/rpe1_kania_cv_risk_boxplot.pdf"
)
if (interactive()) {
  plot_cv_method_comparison(risk_comparison)
}

# Evaluate Rgamma_star for every beta_gamma from OLS to CD.
path_gamma <- unique(c(
  gamma_grid,
  exp(seq(log(max(gamma_grid) + 0.01), log(1e10), length.out = 400))
))
training_risk_path <- compute_training_risk_path(
  data_train, path_gamma, gamma_star
)
plot_training_risk_path(
  training_risk_path, gamma_star,
  "visualization/rpe1_training_risk_path.pdf"
)

path_summary <- data.frame(
  estimator = c("OLS", "CR", "CD"),
  gamma = c(0, gamma_star, Inf),
  training_risk_path[c(
    1,
    which.min(abs(training_risk_path$gamma - gamma_star)),
    nrow(training_risk_path)
  ), c("Rdelta_plus", "Rplus", "Rgamma")]
)
print(path_summary, row.names = FALSE)

cat("Selected gamma:", gamma_star, "\n")
cat("Selected gamma with Kania CV:", gamma_kania, "\n")
cat("Training cells:", nrow(data_train$Xe), "and", nrow(data_train$Xo), "\n")
cat(
  "Training environment 1 interventions:",
  paste(first_training_interventions, collapse = ", "), "\n"
)
cat(
  "Training environment 2 interventions:",
  paste(second_training_interventions, collapse = ", "), "\n"
)
cat("Unseen environments:", nrow(risk_by_environment), "\n")
