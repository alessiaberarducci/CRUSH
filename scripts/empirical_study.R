###########################################
####### EMPIRICAL STUDY ###################
####### NON ALIGNED SHIFTS ################
###########################################


#################################################
####### beta_gamma related to abs(R_\Delta)######
#################################################

source("./datautil.R")
source("./cd.R")
source("./ols.R")
source("./measurements.R")
source("./libraries.R")

set.seed(123) 
n<-100000
X1<- rnorm(n)
a11<-rnorm(n,0,1)
X1<- X1+a11
ye<- X1+ rnorm(n)
X2<- -ye+rnorm(n)
Xe<-matrix(c(X1,X2),ncol=2)
betaXe<-lm(ye~Xe-1)$coef
betaXe

p <- 2
B <- matrix(c(0, 0, 0,
              0, 0, -1,
              1, 0, 0),
            nrow = p+1, ncol = p+1, byrow = TRUE)


I <- diag(p+1)
Px <- matrix(c(1, 0, 0,
               0, 1, 0),  
             nrow = p, ncol = p+1, byrow = TRUE)
Py <- matrix(c(0, 0, 1), nrow = 1, ncol = p+1)

data_e <- list(Xe = Xe, ye = ye)

X1<- rnorm(n)
yo<- X1+rnorm(n)
a22<-rnorm(n,0,1)
X2<- -yo+rnorm(n)+a22 
Xo<- matrix(c(X1,X2),ncol=2)


betaXo<-lm(yo~Xo-1)$coef
betaXo



data_o <- list(Xo = Xo, yo = yo)

data_combined <- c(data_e, data_o)
m <- moments(data_combined)
beta_pa<-compute_cd(m)
beta_ols<-compute_ols(m)


#gamma<-seq(0,50,length=3000)
gamma<-seq(0,10,length.out=5000)
gamma<-c(gamma,100,1000,10000)
n <- length(gamma)

dims <- dim(m$Gplus)


beta_gamma <- matrix(NA, nrow = dims[1], ncol = n)

if (difference(data_combined, beta_ols) < 0) {
  for (i in seq_along(gamma)) {
    Gtemp <- m$Gplus - (gamma[i] * m$Gdelta)
    Ztemp <- m$Zplus - (gamma[i] * m$Zdelta)
    beta_gamma[, i] <- solve(Gtemp, Ztemp)
  }
} else {
  for (i in seq_along(gamma)) {
    Gtemp <- m$Gplus + (gamma[i] * m$Gdelta)
    Ztemp <- m$Zplus + (gamma[i] * m$Zdelta)
    beta_gamma[, i] <- solve(Gtemp, Ztemp)
  }
}


dims <- dim(m$Gplus)

Rdelta_vec<-rep(0,length(gamma))
Rplus_vec<-rep(0,length(gamma))

m$Gdelta <- rbind(cbind(m$Gdelta, rep(0, m$p)),
                  rep(0, m$p + 1))

m$Gplus <- rbind(cbind(m$Gplus, rep(0, m$p)),
                 rep(0, m$p + 1))


m$Zdelta <- c(m$Zdelta, 0)
m$Zplus  <- c(m$Zplus, 0)

for (i in seq_along(gamma)){
  Rdelta_vec[i]<-difference(data_combined,matrix(beta_gamma[,i],nrow=2,ncol=1))
  Rplus_vec[i]<-in_sample_risk(data_combined,matrix(beta_gamma[,i],nrow=2,ncol=1))
}


######################
####### Results N1 ###
######################

plot(gamma,Rdelta_vec,type='l',xlab='γ',ylab='RΔ',ylim=c(-1,2),xlim=c(-1,6),main='Behaviour of RΔ across γ values')
abline(h=0,col='red')


Rdelta_beta_pa<-difference(data_combined,matrix(beta_pa,nrow=2,ncol=1))
Rplus_beta_pa<-in_sample_risk(data_combined,beta_pa) 




beta1<-seq(-10,10, length.out=100)
beta2<-seq(-10,10,length.out=100)
beta<-expand.grid(beta1,beta2)

R_delta<-matrix(0,nrow=length(beta1),ncol=length(beta2))
for(i in 1:length(beta1)){
  for(j in 1:length(beta2)){
    R_delta[i,j]<-difference(data_combined,matrix(c(beta1[i],beta2[j]),nrow=2,ncol=1))
  }
}

plot(beta_gamma[1,],beta_gamma[2,],type='l',xlab='βγ1',ylab='βγ2',xlim=c(0,1.1),ylim=c(-0.8,1.2),main='βγ-')

contour(beta1,beta2,R_delta,level=0,add=TRUE)





grid_n <- 80
beta1_vals <- seq(-2, 2,length.out = grid_n)
beta2_vals <- seq(-2, 2, length.out = grid_n)

grid <- expand.grid(beta1 = beta1_vals, beta2 = beta2_vals)
n_points <- nrow(grid)

Rplus <- numeric(n_points)
Rdelta <- numeric(n_points)
gamma<-seq(0,50,length=300)

for (i in 1:n_points) {
  current_param <- matrix(c(grid$beta1[i], grid$beta2[i]), nrow = 2, ncol = 1)
  Rplus[i] <- in_sample_risk(data_combined, current_param)
  Rdelta[i] <- difference(data_combined, current_param)
}


palette<-brewer.pal(9,"Dark2")
palette2<-brewer.pal(7,"BrBG")
display.brewer.pal(7,"Accent")


Rdelta_matrix <- matrix(Rdelta, nrow = grid_n, ncol = grid_n, byrow = TRUE)
Rplus_matrix  <- matrix(Rplus, nrow = grid_n, ncol = grid_n, byrow = TRUE)

p <- plot_ly()
p <- p %>% add_surface(   x = ~beta1_vals,   y = ~beta2_vals,   z = ~Rdelta_matrix,   colorscale = list(c(0, "yellow"), c(1, palette[2])),   showscale = TRUE,   colorbar = list(title = "RΔ ", tickfont = list(color =  palette[2]),                   titlefont = list(color =  palette[2])),   contours = list(     x = list(show = TRUE, color = "black"),     y = list(show = TRUE, color = "black"),     z = list(show = TRUE, color = "black")   ) ) 

p <- p %>% add_surface(
  x = beta1_vals,
  y = beta2_vals,
  z = Rplus_matrix,
  colorscale = list(c(0,"cornflowerblue"), c(1, palette2[7])),
  opacity = 0.8,  
  showscale = TRUE,
  colorbar = list(title = "R+", tickfont = list(color = "2543B1"),
                  titlefont = list(color = "2543B1")),
  contours = list(
    x = list(show = TRUE, color = "black"),
    y = list(show = TRUE, color = "black"),
    z = list(show = TRUE, color = "black")
  )
)
p

Rplus_beta_ols<-in_sample_risk(data_combined, beta_ols)
p <- p %>% add_markers(
  x = beta_pa[1],
  y = beta_pa[2],
  z = Rplus_beta_pa,
  marker = list(size = 8, color = "darkblue"),
  name = paste("βCP ")
)


p <- p %>% add_markers(
  x = beta_ols[1],
  y = beta_ols[2],
  z = Rplus_beta_ols,
  marker = list(size = 8, color = "red"),
  name = paste("βOLS")
)

p <- p %>% layout(
  title = list(
    x = 0.5,                 
    y = 0.92,               
    xanchor = "center",       
    yanchor = "top",          
    font = list(size = 18)
  ),
  scene = list(
    xaxis = list(title = "β1"),
    yaxis = list(title = "	β2"),
    zaxis = list(title = "R-values (RΔ & R+)", range = c(-3, 20))
  ), 
  legend = list(
    orientation = "h",   
    x = 0.3,             
    y = -0.1           
  ),
  colorbar = list(
    title = "RΔ",
    tickfont = list(color = "black"),
    titlefont = list(color = "black"),
    x = 0.2,    
    y = 0.2,     
    len = 0.5     
  ),
  colorbar = list(
    title = "RΔ",
    tickfont = list(color = palette[2]),
    titlefont = list(color = palette[2]),
    x = 0.5,   
    y = 0.2,      
    len = 0.5     
  ))
p




################################################
####### beta_\gamma^G related to R_\Delta^G ####
################################################

set.seed(123) 
n<-100000
X1<- rnorm(n)
a11<-rnorm(n,0,1)
X1<- X1+a11
ye<- X1+ rnorm(n)
X2<- -ye+rnorm(n)
Xe<-matrix(c(X1,X2),ncol=2)
betaXe<-lm(ye~Xe-1)$coef
betaXe

m$Fdelta<-matrix(c(var(a11),0,0,0,-var(a22),0,0,0,0),nrow=3,ncol=3)
Fdelta<-m$Fdelta
eig_F <- eigen(Fdelta)
V <- eig_F$vectors
Lambda <- diag(abs(eig_F$values))
F_Delta_plus <- V %*% Lambda %*% t(V)

p <- 2


B <- matrix(c(0, 0, 0,
              0, 0, -1,
              1, 0, 0),
            nrow = p+1, ncol = p+1, byrow = TRUE)


I <- diag(p+1)

Px <- matrix(c(1, 0, 0,
               0, 1, 0),  
             nrow = p, ncol = p+1, byrow = TRUE)


Py <- matrix(c(0, 0, 1), nrow = 1, ncol = p+1)
ub=Py %*% solve(I-B) %*% t(Px)

data_e <- list(Xe = Xe, ye = ye)



X1<- rnorm(n)
yo<- X1+rnorm(n)
a22<-rnorm(n,0,1)
X2<- -yo+rnorm(n)+a22 
Xo<- matrix(c(X1,X2),ncol=2)

betaXo<-lm(yo~Xo-1)$coef
betaXo



data_o <- list(Xo = Xo, yo = yo)

data_combined <- c(data_e, data_o)
m <- moments(data_combined)
beta_cp<-compute_cd(m)

beta_ols<-compute_ols(m)


#gamma<-seq(0,500,length=30000)
gamma<-seq(0,10,length.out=5000)
gamma<-c(gamma,100,1000,10000)

n <- length(gamma)


dims <- dim(m$Gplus)


beta_gammaG <- matrix(NA, nrow = dims[1], ncol = n)
eigen_Gdelta <- eigen(m$Gdelta)

V <- eigen_Gdelta$vectors              
Lambda <- diag(abs(eigen_Gdelta$values)) 

G_delta_plus <- V %*% Lambda %*% t(V)
m$Gdelta_plus <- G_delta_plus

for (i in seq_along(gamma)) {
  Gtemp <- m$Gplus +  gamma[i] * m$Gdelta_plus
  Ztemp <- m$Zplus +  gamma[i] * m$Gdelta_plus %*% beta_cp
  beta_gammaG[, i] <- solve(Gtemp, Ztemp)
}


dims <- dim(m$Gplus)

Rdelta_g<-rep(0,length(gamma))
m$Gdelta <- rbind(cbind(m$Gdelta, rep(0, m$p)),
                  rep(0, m$p + 1))

m$Gplus <- rbind(cbind(m$Gplus, rep(0, m$p)),
                 rep(0, m$p + 1))



m$Zdelta <- c(m$Zdelta, 0)
m$Zplus  <- c(m$Zplus, 0)


for (i in seq_along(gamma)){
  Rdelta_g[i]<-t(beta_gammaG[,i]-beta_cp) %*% m$Gdelta_plus %*% (beta_gammaG[,i]-beta_cp)
}



m$Gdelta_plus<-rbind(cbind(m$Gdelta_plus, rep(0, m$p)),
                     rep(0, m$p + 1))





plot(gamma,Rdelta_vec,type='l',xlab='γ',ylab='RΔ',xlim=c(0,10))
abline(h=0,col='red')

plot(gamma,Rdelta_g,type='l',xlab='γ',ylab='RΔG',ylim=c(-0.5,1.8),xlim=c(0,10),main='Variation of RΔG across γ values' )
abline(h=0,col='red')

plot(beta_gammaG[1,],beta_gammaG[2,],type='l',xlab='βγG1',ylab='βγG2',xlim=c(-0.4,1.5),ylim=c(-1,1),main='βγG')





beta1<-seq(-10,10, length.out=100)
beta2<-seq(-10,10,length.out=100)
beta<-expand.grid(beta1,beta2)
R_delta<-matrix(0,nrow=length(beta1),ncol=length(beta2))
for(i in 1:length(beta1)){
  for(j in 1:length(beta2)){
    R_delta[i,j]<-difference(data_combined,matrix(c(beta1[i],beta2[j]),nrow=2,ncol=1))
  }
}

contour(beta1,beta2,R_delta,level=0,add=TRUE)


##############################################
###### beta_gamma related to R_\Delta^+ ######
##############################################

set.seed(123) 
n<-100000
X1<- rnorm(n)
a11<-rnorm(n,0,1)
X1<- X1+a11
ye<- X1+ rnorm(n)
X2<- -ye+rnorm(n)
Xe<-matrix(c(X1,X2),ncol=2)
betaXe<-lm(ye~Xe-1)$coef
betaXe

m$Fdelta<-matrix(c(var(a11),0,0,0,-var(a22),0,0,0,0),nrow=3,ncol=3)
Fdelta<-m$Fdelta
eig_F <- eigen(Fdelta)
V <- eig_F$vectors
Lambda <- diag(abs(eig_F$values))
F_Delta_plus <- V %*% Lambda %*% t(V)

p <- 2


B <- matrix(c(0, 0, 0,
              0, 0, -1,
              1, 0, 0),
            nrow = p+1, ncol = p+1, byrow = TRUE)


I <- diag(p+1)

Px <- matrix(c(1, 0, 0,
               0, 1, 0),  
             nrow = p, ncol = p+1, byrow = TRUE)


Py <- matrix(c(0, 0, 1), nrow = 1, ncol = p+1)
ub=Py %*% solve(I-B) %*% t(Px)

data_e <- list(Xe = Xe, ye = ye)



X1<- rnorm(n)
yo<- X1+rnorm(n)
a22<-rnorm(n,0,1)
X2<- -yo+rnorm(n)+a22 
Xo<- matrix(c(X1,X2),ncol=2)

betaXo<-lm(yo~Xo-1)$coef
betaXo



data_o <- list(Xo = Xo, yo = yo)

data_combined <- c(data_e, data_o)
m <- moments(data_combined)
beta_cp<-compute_cd(m)

beta_ols<-compute_ols(m)


#gamma<-seq(0,500,length=30000)
gamma<-seq(0,10,length.out=5000)
gamma<-c(gamma,100,1000,10000)

n <- length(gamma)


dims <- dim(m$Gplus)


beta_gamma_all <- matrix(0, nrow = 2, ncol = length(gamma))

for (i in seq_along(gamma)) {
  A <- 2 * gamma[i] * u_x %*% F_Delta_plus %*% t(u_x)
  b <- m$Zplus + 2 * gamma[i]* u_x %*% F_Delta_plus %*% t(u_y)
  beta_gamma_all[, i] <- solve(m$Gplus + A, b)
}


dims <- dim(m$Gplus)


m$Gdelta <- rbind(cbind(m$Gdelta, rep(0, m$p)),
                  rep(0, m$p + 1))

m$Gplus <- rbind(cbind(m$Gplus, rep(0, m$p)),
                 rep(0, m$p + 1))



m$Zdelta <- c(m$Zdelta, 0)
m$Zplus  <- c(m$Zplus, 0)




u_x <- Px %*% solve(I - B)
u_y <- Py %*% solve(I - B)



# Preallocate vector
R_delta_plus_gamma <- numeric(length(gamma))
IB_inv <- solve(I - B)

for (i in seq_along(gamma)) {
  beta_gamma_i <- beta_gamma_all[, i, drop = FALSE]  
  u_beta <- (Py - t(beta_gamma_i) %*% Px) %*% IB_inv
  R_delta_plus_gamma[i] <- u_beta %*% F_Delta_plus %*% t(u_beta)
}



###########################
##### Comparison ##########
###########################


plot(beta_gammaG[1, ], beta_gammaG[2, ],
     type = 'l', xlab = 'βγ1', ylab = 'βγ2',
     xlim = c(-0.4, 1.5), ylim = c(-1, 1),,pch=1,cex=1.2,lwd=3)

contour(beta1, beta2, R_delta, level = 0, add = TRUE)


lines(beta_gamma[1, ], beta_gamma[2, ],
       type = 'l', col = 'darkorange', pch = 3, cex = 0.5,lwd=3)

lines(beta_gamma_all[1, ], beta_gamma_all[2, ],
       type = 'l', col = 'blue', pch = 3, cex = 0.5,lwd=3,lty=22)


legend("topleft",
       legend = c("βγG", "βγ-","βγ"),
       col = c("black", "darkorange","blue"),
       lwd=c(3,3,3),
       lty=c(1,1,22),
       cex = 0.9,
       pt.cex = 0.7,
       x.intersp = 0.3,
       y.intersp = 0.3,
       bty = "o")




plot(gamma, Rdelta_g,
     type = 'l', xlab = 'γ', ylab = 'RΔG, |RΔ| & RΔ+',
     xlim = c(0, 20), ylim = c(-0.01, 0.4),pch=1,cex=1.2,lwd=3)


lines(gamma, abs(Rdelta_vec),
       type = 'l', col = 'darkorange', pch = 3, cex = 0.5,lwd=3)

lines(gamma,R_delta_plus_gamma,
       type = 'l', col = 'blue', pch = 3, cex = 0.5,lwd=3,lty=22)


legend("topright",
       legend = c("RΔG", "|RΔ|","RΔ+"),
       col = c("black", "darkorange","blue"),
       lwd=c(3,3,3),
       lty=c(1,1,22),
       cex = 0.9,
       pt.cex = 0.7,
       x.intersp = 0.3,
       y.intersp = 0.3,
       bty = "o")




