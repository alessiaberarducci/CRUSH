#######################################
###### Cross Validation Function  #####
#######################################

source("./datautil.R")
source("./cd.R")
source("./ols.R")
source("./measurements.R")
source("./libraries.R")



set.seed(123)

split_idx <- function(k, n)
  sample(rep(seq_len(k), length.out = n))

subset_data <- function(dat, sele, selo) {
  list(
    Xe = dat$Xe[sele, , drop = FALSE],
    ye = dat$ye[sele],
    Xo = dat$Xo[selo, , drop = FALSE],
    yo = dat$yo[selo]
  )
}


build_matrices <- function(split) {
  m <- moments(split)
  eg <- eigen(m$Gdelta, symmetric = TRUE)
  G_delta_plus <- eg$vectors %*%
    diag(abs(eg$values)) %*%
    t(eg$vectors)
  list(
    Gplus        = m$Gplus,
    Gdelta_plus  = G_delta_plus,
    Zplus        = m$Zplus
  )
}

estimate_beta_cp <- function(split) {
  m <- moments(split)
  compute_cd(m)
}



calc_Rdelta <- function(split_val, beta_hat) {
  m <- build_matrices(split_val)
  t(beta_hat - estimate_beta_cp(split_val)) %*% m$Gdelta_plus %*% (beta_hat -
                                                                     estimate_beta_cp(split_val))
}


make_estimator <- function(gam) {
  function(split) {
    m       <- build_matrices(split)
    beta_cp <- estimate_beta_cp(split)
    
    Gtemp   <- m$Gplus + gam * m$Gdelta_plus
    Ztemp   <- m$Zplus + gam * m$Gdelta_plus %*% beta_cp
    
    solve(Gtemp, Ztemp)
  }
}


fit_cross_validation <- function(folds, f, data_train, data_val) {
  vapply(seq_len(folds), function(i) {
    beta_hat <- f(data_train[[i]])
    calc_Rdelta(data_val[[i]], beta_hat)
  }, numeric(1))
}


cross_validation <- function(folds, data, estimators) {
  idx <- split_idx(folds, nrow(data$Xe))
  ge <- go <- idx
  data_train <- vector("list", folds)
  data_val   <- vector("list", folds)
  
  for (i in seq_len(folds)) {
    data_train[[i]] <- subset_data(data, sele = (ge != i), selo = (go != i))
    data_val[[i]]   <- subset_data(data, sele = (ge == i), selo = (go == i))
  }
  cv_mat <- vapply(
    estimators,
    fit_cross_validation,
    folds      = folds,
    data_train = data_train,
    data_val   = data_val,
    FUN.VALUE  = numeric(folds)
  )
  
  mean_abs_rdelta <- colMeans((cv_mat))
  best_idx        <- which.min(mean_abs_rdelta)
  
  list(
    best_gamma      = names(estimators)[best_idx],
    mean_abs_rdelta = mean_abs_rdelta,
    foldwise_rdelta = cv_mat
  )
}
