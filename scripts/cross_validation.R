#################################
#### Cross Validation ###########
#################################

source("./datautil.R")
source("./cd.R")
source("./ols.R")
source("./measurements.R")
source("./libraries.R")



set.seed(123)                                 

split_idx <- function(k, n) sample(rep(seq_len(k), length.out = n))

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



calc_Rdelta<-function(split_val, beta_hat) {
  m<-build_matrices(split_val)
  t(beta_hat-estimate_beta_cp(split_val)) %*% m$Gdelta_plus %*% (beta_hat-estimate_beta_cp(split_val))
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
  idx<-split_idx(folds,nrow(data$Xe))
  ge<-go<-idx
  data_train <- vector("list", folds)
  data_val   <- vector("list", folds)
  
  for (i in seq_len(folds)) {
    data_train[[i]] <- subset_data(data, sele = (ge != i), selo = (go != i))
    data_val[[i]]   <- subset_data(data, sele = (ge == i),  selo = (go == i))
  }
  cv_mat <- vapply(estimators, fit_cross_validation,
                   folds      = folds,
                   data_train = data_train,
                   data_val   = data_val,
                   FUN.VALUE  = numeric(folds))
  
  mean_abs_rdelta <- colMeans((cv_mat))
  best_idx        <- which.min(mean_abs_rdelta)
  
  list(
    best_gamma      = names(estimators)[best_idx],
    mean_abs_rdelta = mean_abs_rdelta,
    foldwise_rdelta = cv_mat
  )
}

###################################
###### Example usage  #############                        
###################################


n <- 700

a11 <- rnorm(n)            
x1e <- rnorm(n) + a11      
ye  <- x1e + rnorm(n)      
x2e <- -ye + rnorm(n)      
Xe  <- cbind(x1e, x2e)


x1o <- rnorm(n)  
yo  <- x1o + rnorm(n)  
a22 <- rnorm(n)           
x2o <- -yo + rnorm(n) + a22  
Xo  <- cbind(x1o, x2o)
data <- list(Xe = Xe, ye = ye, Xo = Xo, yo = yo)

gamma_grid <- seq(0,100, by = 0.5)
estimators <- setNames(lapply(gamma_grid, make_estimator),
        gamma_grid)

res <- cross_validation(folds = 5, data = data, estimators = estimators)

cat("Best γ:", res$best_gamma, "\n")
best_gamma<-res$best_gamma
cat("Mean |Rδ|:", res$mean_abs_rdelta, "\n")





mean_abs_rdelta     <- res$mean_abs_rdelta
foldwise_rdelta_mat <- res$foldwise_rdelta            # rows = γ, cols = fold


library(reshape2)
library(ggplot2)

## “Average” curve
df_avg <- data.frame(
  gamma    = gamma_grid,
  abs_risk = mean_abs_rdelta,
  what     = "Average"
)

m         <- res$foldwise_rdelta   

df_fold   <- melt(m, varnames = c("fold", "gamma"), value.name = "abs_risk")

library(ggplot2)

ggplot() +
  geom_line(data = df_fold,
            aes(x = gamma, y = abs_risk, group = fold, linetype = "Fold"),
            size = 0.8) +
  
  geom_line(data = df_avg,
            aes(x = gamma, y = abs_risk, linetype = "Average"),
            size = 0.8) +
  geom_point(data = df_avg,
             aes(x = gamma, y = abs_risk)) +
  
  geom_vline(xintercept = as.numeric(best_gamma),
             linetype   = "dashed",
             colour     = "red") +
  scale_linetype_manual(name = NULL,           
                        values = c(Average = "solid", Fold = "dotted")) +
  labs(x = "γ",
       y = "RΔG") +
  theme_minimal(base_size = 12)+
  ylim(0.,0.3)





