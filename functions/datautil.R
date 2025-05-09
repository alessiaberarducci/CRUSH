##################################################
# Compute moments ################################
##################################################

moments <- function(data) {
  
  ne <- dim(data$Xe)[1]
  XYe <- t(data$Xe) %*% data$ye / ne
  XXe <- t(data$Xe) %*% data$Xe / ne
  
  no <- dim(data$Xo)[1]
  XYo <- t(data$Xo) %*% data$yo / no
  XXo <- t(data$Xo) %*% data$Xo / no
  
  Zplus <- XYe + XYo
  Gplus <- XXe + XXo
  
  Zdelta <- XYe - XYo
  Gdelta <- XXe - XXo
  
  return(list(XYe = XYe,
              XXe = XXe,
              ne = ne,
              XYo = XYo,
              XXo = XXo,
              no = no,
              Zdelta = Zdelta,
              Gdelta= Gdelta,
              Gplus = Gplus,
              Zplus = Zplus,
              p = dim(data$Xe)[2]))
}
