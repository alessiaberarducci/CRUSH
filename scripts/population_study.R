############################################
###### POPULATION STUDY ####################
##### NON ALIGNED SHIFTS ###################
############################################

Zdelta<- matrix(c(1, -1 ), nrow =2, byrow = TRUE)
Gdelta<-matrix(c(1,-1,-1,0),nrow=2,byrow=TRUE)
Zplus<-matrix(c(3,-5),ncol=1)
Gplus<-matrix(c(3,-3,-3,8),nrow=2,byrow=TRUE)



beta1 <- seq(0, 2, length.out = 100)
beta2 <- seq(-1, 1, length.out = 100)
beta <- expand.grid(beta1 = beta1, beta2 = beta2)

Rdelta <- matrix(0, nrow = length(beta1), ncol = length(beta2))
F_Delta_<- matrix(c(1,0,0,0,-1,0,0,0,0),byrow=TRUE,nrow=3)
F_Delta_plus<-matrix(c(1,0,0,0,1,0,0,0,0),byrow=TRUE,nrow=3)
R_delta_plus<- matrix(0,nrow=length(beta1),ncol=length(beta2))

p<-2
B <- matrix(c(0, 0, 0,
              0, 0, -1,
              1, 0, 0),
            nrow = p+1, ncol = p+1, byrow = TRUE)
I <- diag(p+1)
Px <- matrix(c(1, 0, 0,
               0, 1, 0),  
             nrow = p, ncol = p+1, byrow = TRUE)
Py <- matrix(c(0, 0, 1), nrow = 1, ncol = p+1)

for (i in 1:length(beta1)) {
  for (j in 1:length(beta2)) {
    beta_vec <- c(beta1[i], beta2[j])
    term1 <- 2 * sum(Zdelta * beta_vec) 
    term2 <- t(beta_vec) %*% Gdelta %*% beta_vec  
    Rdelta[i, j] <- 1 - term1 + term2
    ub<-(Py-beta_vec%*%Px)%*%solve(I-B)
    R_delta_plus[i,j]<-ub%*%F_Delta_plus%*%t(ub) 
  }
}

# Abs(Rdelta) & Rdelta+

library(plotly)
p <- plot_ly()
p <- p %>% add_surface(   x = ~beta1,   y = ~beta2,   z = ~abs(Rdelta),  colorscale = list(c(0, "orange"), c(1, "red")),   showscale = TRUE,   colorbar = list(title = "|RΔ| ", tickfont = list(color =  "red"),                   titlefont = list(color =  'red')),   contours = list(     x = list(show = TRUE, color = "black"),     y = list(show = TRUE, color = "black"),     z = list(show = TRUE, color = "black")   ) ) 
p
p<-p %>% add_surface(   x = ~beta1,   y = ~beta2,   z = ~(R_delta_plus),   colorscale = list(c(0, "lightblue"), c(1,'blue' )),   showscale = TRUE,   colorbar = list(title = "RΔ+ ", tickfont = list(color ='blue'  ),    titlefont = list(color =  'blue')),   contours = list(     x = list(show = TRUE, color = "black"),     y = list(show = TRUE, color = "black"),     z = list(show = TRUE, color = "black" )   ) ) 
p <- p %>% add_markers(
  x = 1,
  y = 0,
  z = 0,
  marker = list(size = 8, color = "blue"),
  name = paste("βCP ")
)
p <- p %>% layout(
  title = list(
    text = "Surface Plot of |RΔ| and RΔ+",
    x = 0.5,                  
    y = 0.92,                 
    xanchor = "center",       
    yanchor = "top",          
    font = list(size = 18)   
  ),
  scene = list(
    xaxis = list(title = "β1"),
    yaxis = list(title = "	β2"),
    zaxis = list(title = "R-values (|RΔ| & RΔ+)")
  ), 
  legend = list(
    orientation = "h",   
    x = 0.3,             
    y = -0.1             
  ))
p



# Theoretical Contour 

contour(beta1, beta2, Rdelta, level=0,xlab="β1",ylab="β2", main="Contour Plot of Population RΔ")


# Curve
gamma<-seq(0,500,length=30000)


beta_min<-matrix(NA,nrow=2,ncol=length(gamma))
beta_plu<-matrix(NA,nrow=2,ncol=length(gamma))



R_gamma_plu1 <- rep(NA_real_, length(gamma))
R_gamma_min1 <- rep(NA_real_, length(gamma))
R_gamma_plu2 <- rep(NA_real_, length(gamma))
R_gamma_min2 <- rep(NA_real_, length(gamma))
R_delta_min <- rep(NA_real_, length(gamma))

for(i in 1:length(gamma)) {
  beta_min[,i] <- solve(Gplus - gamma[i] * Gdelta, Zplus - gamma[i] * Zdelta)
  beta_plu[,i] <- solve(Gplus + gamma[i] * Gdelta, Zplus + gamma[i] * Zdelta)
  
  R_plus <- 5 - 2 * t(Zplus) %*% beta_plu[,i] + t(beta_plu[,i]) %*% Gplus %*% beta_plu[,i]
  R_plus_min <- 5 - 2 * t(Zplus) %*% beta_min[,i] + t(beta_min[,i]) %*% Gplus %*% beta_min[,i]
  R_delta <- 1 - 2 * t(Zdelta) %*% beta_plu[,i] + t(beta_plu[,i]) %*% Gdelta %*% beta_plu[,i]
  R_delta_min[i] <- as.numeric(1 - 2 * t(Zdelta) %*% beta_min[,i] + t(beta_min[,i]) %*% Gdelta %*% beta_min[,i])
  
  R_gamma_min1[i] <- 0.5 * as.numeric(R_plus_min) - gamma[i] * R_delta_min[i]
  R_gamma_plu1[i] <- 0.5 * as.numeric(R_plus) + gamma[i] * as.numeric(R_delta)
  R_gamma_min2[i] <- 0.5 * as.numeric(R_plus_min) + gamma[i] * R_delta_min[i]
  R_gamma_plu2[i] <- 0.5 * as.numeric(R_plus) - gamma[i] * as.numeric(R_delta)
}

# WORKING
plot(gamma,R_gamma_min1,type='l',ylab='Rγ_',xlab='γ',ylim=c(0.4,1),lwd=3)
abline(h=0.6,col='blue')
abline(h=1,col='red')

plot(gamma,R_delta_min,type='l',ylab='Rdelta_min',xlab='γ',lwd=3)

R_plus_ols <- 5- 2*t(Zplus)%*% t(t(beta_plu[,1]))+ t(beta_plu[,1])%*%Gplus%*%t(t(beta_plu[,1])) 


beta_ols<-solve(Gplus,Zplus)
beta_cp<-solve(Gdelta,Zdelta)


# NOT WORKING
plot(gamma,R_gamma_plu1,type='l',xlab='γ',ylab='Rgamma_plus',ylim=c(0.4,1))
abline(h=0.5,col='blue')
abline(h=1,col='red')

# NOT WORKING
plot(gamma,R_gamma_min2,type='l',ylab='Rgamma_minus',xlab='γ')
abline(h=0.6,col='blue')
abline(h=1,col='red')

# NOT WORKING 
plot(gamma,R_gamma_plu2,type='l',xlab='γ',ylab='Rgamma_plus',ylim=c(0.4,1))
abline(h=0.5,col='blue')
abline(h=1,col='red')



# BETA MINUS WITH THEORETICAL CONTOUR 

plot(beta_min[1,],beta_min[2,],type='p',xlim=c(0.4,1.5),ylim=c(-1,1.2),xlab='βγ1',ylab='βγ2',cex=0.7,main=' βγ-')
contour(beta1, beta2, Rdelta, level=0,xlab="β1",ylab="β2", add=TRUE)



###### SWAP SIGN ####################


beta1 <- seq(-10, 10, length.out = 1000)
beta2 <- seq(-10, 10, length.out = 1000)
beta <- expand.grid(beta1 = beta1, beta2 = beta2)
R_delta <- matrix(0, nrow = length(beta1), ncol = length(beta2))


for (i in 1:length(beta1)) {
  for (j in 1:length(beta2)) {
    beta_vec <- c(beta1[i], beta2[j])
    term1 <- 2 * sum(Zdelta * beta_vec) 
    term2 <- t(beta_vec) %*% Gdelta %*% beta_vec  
    R_delta[i, j] <- -1 - term1 + term2
  }
}


contour(beta1, beta2, R_delta, level=0,xlab="β1",ylab="β2", main="Contour Plot of Population RΔ")


gamma<-seq(0,500,length=30000)
Zdelta<- matrix(c(-1, 1 ), nrow =2, byrow = TRUE)
Gdelta<-matrix(c(-1,1,1,0),nrow=2,byrow=TRUE)
Zplus<-matrix(c(3,-5),ncol=1)
Gplus<-matrix(c(3,-3,-3,8),nrow=2,byrow=TRUE)


beta_min<-matrix(NA,nrow=2,ncol=length(gamma))
beta_plu<-matrix(NA,nrow=2,ncol=length(gamma))
R_gamma_plu1 <- rep(NA_real_, length(gamma))
R_gamma_min1 <- rep(NA_real_, length(gamma))
R_gamma_plu2 <- rep(NA_real_, length(gamma))
R_gamma_min2 <- rep(NA_real_, length(gamma))
R_delta_min <- rep(NA_real_, length(gamma))

for(i in 1:length(gamma)) {
  beta_min[,i] <- solve(Gplus - gamma[i] * Gdelta, Zplus - gamma[i] * Zdelta)
  beta_plu[,i] <- solve(Gplus + gamma[i] * Gdelta, Zplus + gamma[i] * Zdelta)
  
  R_plus <- 5 - 2 * t(Zplus) %*% beta_plu[,i] + t(beta_plu[,i]) %*% Gplus %*% beta_plu[,i]
  R_plus_min <- 5 - 2 * t(Zplus) %*% beta_min[,i] + t(beta_min[,i]) %*% Gplus %*% beta_min[,i]
  R_delta <- -1 - 2 * t(Zdelta) %*% beta_plu[,i] + t(beta_plu[,i]) %*% Gdelta %*% beta_plu[,i]
  R_delta_min[i] <- as.numeric(-1 - 2 * t(Zdelta) %*% beta_min[,i] + t(beta_min[,i]) %*% Gdelta %*% beta_min[,i])
  
  R_gamma_min1[i] <- 0.5 * as.numeric(R_plus_min) - gamma[i] * R_delta_min[i]
  R_gamma_plu1[i] <- 0.5 * as.numeric(R_plus) + gamma[i] * as.numeric(R_delta)
  R_gamma_min2[i] <- 0.5 * as.numeric(R_plus_min) + gamma[i] * R_delta_min[i]
  R_gamma_plu2[i] <- 0.5 * as.numeric(R_plus) - gamma[i] * as.numeric(R_delta)
}
# NOT WORKING
plot(gamma,R_gamma_min1,type='l',ylab='Rγ_',xlab='γ',lwd=3)
abline(h=0.6,col='blue')
abline(h=1,col='red')

# WORKING
plot(gamma,R_gamma_plu1,type='l',xlab='γ',ylab='Rgamma_plus')
abline(h=0.6,col='blue')
abline(h=1,col='red')

# NOT WORKING
plot(gamma,R_gamma_min2,type='l',ylab='Rgamma_minus',xlab='γ')
abline(h=0.6,col='blue')
abline(h=1,col='red')

# NOT WORKING 
plot(gamma,R_gamma_plu2,type='l',xlab='γ',ylab='Rgamma_plus')
abline(h=0.6,col='blue')
abline(h=1,col='red')

# BETA PLUS WITH THEORETICAL CONTOUR

plot(beta_plu[1,],beta_plu[2,],type='p',xlim=c(-2,2),ylim=c(-2,2),xlab='βγ1',ylab='βγ2',cex=0.7,main='Beta_plus and Beta_mins')
contour(beta1, beta2, R_delta, level=0,xlab="β1",ylab="β2", add=TRUE)



