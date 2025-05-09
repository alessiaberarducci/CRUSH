########################################################
# Measurements #########################################
########################################################

set.seed(123)
 
risk_vector <- function(X, y, beta) {
  return((y - X %*% beta)^2)
}

risk <- function(X, y, beta) {
  n <- dim(X)[1]
  return(sum((y - X %*% beta)^2) / n)
}

# Rplus
in_sample_risk <- function(data, beta) {
  return(risk(X = data$Xe, y = data$ye, beta = beta) + risk(X = data$Xo, y = data$yo, beta = beta))
}


# Rdelta
difference <- function(data, beta) {
  return(risk(data$Xe, data$ye, beta) - risk(data$Xo, data$yo, beta))
}


abs_difference <- function(data, beta) {
  return(abs(difference(data, beta)))
}

measure <- function(data, beta, metric) {
  
  n <- dim(beta)[1]
  measurements <- numeric(n)
  
  for (i in seq(1, n, 1)) {
    
    b <- t(t(beta[i,]))
    
    measurements[[i]] <- metric(data, b)
  }
  
  return(measurements)
  
}

