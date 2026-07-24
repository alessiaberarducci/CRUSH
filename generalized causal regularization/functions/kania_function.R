# Cross-validation criterion from Kania and Wit

kania_validation_loss <- function(data, beta_hat) {
  # Compare the two validation MSEs directly and take their absolute difference.
  residual_e <- data$ye - data$Xe %*% beta_hat
  residual_o <- data$yo - data$Xo %*% beta_hat
  abs(mean(residual_e^2) - mean(residual_o^2))
}

kania_cross_validation <- function(data, estimators, fold_assignment) {
  # Reuse the same folds as the main CV so only the validation loss changes.
  fold_e <- fold_assignment$environment_e
  fold_o <- fold_assignment$environment_o
  folds <- length(unique(fold_e))

  if (length(unique(fold_o)) != folds) {
    stop("The two environments must use the same number of folds.")
  }
  validate_cv_data(data, folds)

  foldwise <- vapply(estimators, function(estimator) {
    vapply(seq_len(folds), function(fold) {
      training <- subset_data(data, fold_e != fold, fold_o != fold)
      validation <- subset_data(data, fold_e == fold, fold_o == fold)
      kania_validation_loss(validation, estimator(training))
    }, numeric(1))
  }, numeric(folds))

  mean_loss <- colMeans(foldwise)
  best <- which.min(mean_loss)
  list(
    best_gamma = names(estimators)[best],
    mean_abs_rdelta = mean_loss,
    foldwise_rdelta = foldwise,
    fold_assignment = fold_assignment
  )
}

