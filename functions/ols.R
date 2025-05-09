##########################################
# OLS ####################################
##########################################

library(CVXR)

compute_ols <- function(m) {
  
  beta <- Variable(rows = m$p)
  objective <- Minimize(sum_squares(m$Gplus %*% beta - m$Zplus))
  problem <- Problem(objective)
  sol <- solve(problem)
  
  return(sol$getValue(beta))
}
