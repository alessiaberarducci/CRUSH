############################################
###### POPULATION STUDY ####################
###### NOT ALIGNED SHIFT ##################
############################################
############################################
###### POPULATION STUDY ####################
############################################

source("./libraries.R")

Zdelta<- matrix(c(1, -1 ), nrow =2, byrow = TRUE)
Gdelta<-matrix(c(1,-1,-1,0),nrow=2,byrow=TRUE)
Zplus<-matrix(c(3,-5),ncol=1)
Gplus<-matrix(c(3,-3,-3,8),nrow=2,byrow=TRUE)


eigen_Gdelta <- eigen(Gdelta)
V <- eigen_Gdelta$vectors              
Lambda <- diag(abs(eigen_Gdelta$values)) 
G_delta_plus <- V %*% Lambda %*% t(V)


beta1 <- seq(0, 2, length.out = 100)
beta2 <- seq(-1, 1, length.out = 100)
beta <- expand.grid(beta1 = beta1, beta2 = beta2)

Rdelta <- matrix(0, nrow = length(beta1), ncol = length(beta2))
F_Delta_<- matrix(c(1,0,0,0,-1,0,0,0,0),byrow=TRUE,nrow=3)
F_Delta_plus<-matrix(c(1,0,0,0,1,0,0,0,0),byrow=TRUE,nrow=3)
R_delta_plus<- matrix(0,nrow=length(beta1),ncol=length(beta2))
R_delta_g<- matrix(0,nrow=length(beta1),ncol=length(beta2))

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
beta_cp<-matrix(c(1,0),nrow=1)
for (i in 1:length(beta1)) {
  for (j in 1:length(beta2)) {
    beta_vec <- c(beta1[i], beta2[j])
    term1 <- 2 * sum(Zdelta * beta_vec) 
    term2 <- t(beta_vec) %*% Gdelta %*% beta_vec  
    Rdelta[i, j] <- 1 - term1 + term2
    ub<-(Py-beta_vec%*%Px)%*%solve(I-B)
    R_delta_plus[i,j]<-ub%*%F_Delta_plus%*%t(ub) 
    R_delta_g[i,j]<-(beta_vec- beta_cp)%*%G_delta_plus%*%t(beta_vec- beta_cp)
  }
}
ub<-(Py-beta_vec%*%Px)%*%solve(I-B)

# Abs(Rdelta) & Rdelta+

p <- plot_ly()
p <- p %>% add_surface(   x = ~beta1,   y = ~beta2,   z = ~abs(Rdelta),  colorscale = list(c(0, "orange"), c(1, "red")),   showscale = TRUE,   colorbar = list(title = "|RΔ| ", tickfont = list(color =  "red"),                   titlefont = list(color =  'red')),   contours = list(     x = list(show = TRUE, color = "black"),     y = list(show = TRUE, color = "black"),     z = list(show = TRUE, color = "black")   ) ) 
p
p<-p %>% add_surface(   x = ~beta1,   y = ~beta2,   z = ~(R_delta_plus),   colorscale = list(c(0, "lightblue"), c(1,'darkblue' )),   showscale = TRUE,   colorbar = list(title = "RΔ+ ", tickfont = list(color ='darkblue'  ),    titlefont = list(color =  'darkblue')),   contours = list(     x = list(show = TRUE, color = "black"),     y = list(show = TRUE, color = "black"),     z = list(show = TRUE, color = "black" )   ) ) 
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




# abs(RDelta) & RDelta

palette<-brewer.pal(9,"Dark2")
palette2<-brewer.pal(7,"BrBG")


p <- plot_ly()
p <- p %>% add_surface(   x = ~beta1,   y = ~beta2,   z = ~Rdelta,   colorscale = list(c(0, "yellow"), c(1, palette[2])),   showscale = TRUE,   colorbar = list(title = "RΔ ", tickfont = list(color =  palette[2]),                   titlefont = list(color =  palette[2])),   contours = list(     x = list(show = TRUE, color = "black"),     y = list(show = TRUE, color = "black"),     z = list(show = TRUE, color = "black")   ) ) 


p
p <- p %>% add_surface(   x = ~beta1,   y = ~beta2,   z = ~abs(Rdelta),  colorscale = list(c(0, "darkorange"), c(1, "red")),   showscale = TRUE,   colorbar = list(title = "|RΔ| ", tickfont = list(color =  "red"),                   titlefont = list(color =  'red')),   contours = list(     x = list(show = TRUE, color = "black"),     y = list(show = TRUE, color = "black"),     z = list(show = TRUE, color = "black")   ) ) 
p

p <- p %>% add_markers(
  x = 1,
  y = 0,
  z = 0,
  marker = list(size = 8, color = "blue"),
  name = paste("βCP ")
)
p <- p %>% layout(
  title = list(
   # text = "Surface Plot of RΔ and |RΔ|",
    x = 0.5,                  
    y = 0.92,                 
    xanchor = "center",       
    yanchor = "top",          
    font = list(size = 18)   
  ),
  scene = list(
    xaxis = list(title = "β1"),
    yaxis = list(title = "	β2"),
    zaxis = list(title = "R-values (RΔ & |RΔ|)")
  ), 
  legend = list(
    orientation = "h",   
    x = 0.3,             
    y = -0.1             
  ))
p


# Rdelta+ and R_delta_G
beta_ols<-solve(Gplus,Zplus)
ub<-(Py-t(beta_ols)%*%Px)%*%solve(I-B)
R_delta_plus_ols<-ub%*%F_Delta_plus%*%t(ub) 
R_delta_plus_ols
R_delta_g_ols<-(t(beta_ols)- beta_cp)%*%G_delta_plus%*%t(t(beta_ols)- beta_cp)
R_delta_g_ols


p <- plot_ly()
p <- p %>% add_surface(   x = ~beta1,   y = ~beta2,   z = ~R_delta_g,  colorscale = list(c(0, "lightgreen"), c(1, "darkgreen")),   showscale = TRUE,   colorbar = list(title = "RΔG", tickfont = list(color =  "darkgreen"),                   titlefont = list(color =  'darkgreen')),   contours = list(     x = list(show = TRUE, color = "black"),     y = list(show = TRUE, color = "black"),     z = list(show = TRUE, color = "black")   ) ) 
p
p<-p %>% add_surface(   x = ~beta1,   y = ~beta2,   z = ~(R_delta_plus),   colorscale = list(c(0, "lightblue"), c(1,'darkblue' )),   showscale = TRUE,   colorbar = list(title = "RΔ+ ", tickfont = list(color ='darkblue'  ),    titlefont = list(color =  'darkblue')),   contours = list(     x = list(show = TRUE, color = "black"),     y = list(show = TRUE, color = "black"),     z = list(show = TRUE, color = "black" )   ) ) 
p <- p %>% add_markers(
  x = 1,
  y = 0,
  z = 0,
  marker = list(size = 8, color = "blue"),
  name = paste("βCP ")
)
p <- p %>% add_markers(
  x = beta_ols[1,],
  y = beta_ols[2,],
  z = R_delta_plus_ols,
  marker = list(size = 8, color = "red"),
  name = paste("βOLS2 ")
)
p <- p %>% add_markers(
  x = beta_ols[1,],
  y = beta_ols[2,],
  z = R_delta_g_ols,
  marker = list(size = 8, color = "yellow"),
  name = paste("βOLS1 ")
)
p <- p %>% layout(
  title = list(
    text = "Surface Plot of RΔG and RΔ+",
    x = 0.5,                  
    y = 0.92,                 
    xanchor = "center",       
    yanchor = "top",          
    font = list(size = 18)   
  ),
  scene = list(
    xaxis = list(title = "β1"),
    yaxis = list(title = "	β2"),
    zaxis = list(title = "R-values (RΔG & RΔ+)")
  ), 
  legend = list(
    orientation = "h",   
    x = 0.3,             
    y = -0.1             
  ))
p








###########


# Theoretical Contour 

contour(beta1, beta2, Rdelta, level=0,xlab="β1",ylab="β2", main="Contour Plot of Population RΔ")





#gamma<-seq(0,500,length=30000)
gamma<-seq(0,100,length.out=5000)
gamma<-c(gamma,100,1000,10000)

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




#########################################
###### NEW beta_gamma ALGORITHM #########
###### EMPIRICAL CASE ##################
#########################################

beta_gamma_G<-matrix(NA,nrow=2,ncol=length(gamma))

gamma<-seq(0,100,length.out=5000)
gamma<-c(gamma,100,1000,10000)
#gamma<-seq(0,500,length=30000)

R_delta_G <- numeric(length(gamma))  

for (i in seq_along(gamma)) {
  Gtemp <- Gplus + gamma[i] * G_delta_plus
  Ztemp <- Zplus + gamma[i] * G_delta_plus %*% t(beta_cp)
  
  beta_gamma_G[, i] <- solve(Gtemp, Ztemp)
  
  beta_i <- beta_gamma_G[, i, drop = FALSE]      
  diff <- beta_i - matrix(beta_cp, ncol = 1)    
  
  R_delta_G[i] <- t(diff) %*% G_delta_plus %*% diff 
}

plot(gamma, R_delta_G, type = 'p', xlab = 'γ', ylab = 'RΔG',main='RΔG vs γ',xlim=c(0,100),cex=0.7)

plot(beta_gamma_G[1, ], beta_gamma_G[2, ], type = 'p', xlim = c(0.4, 1.5), ylim = c(-1, 1.2), xlab = 'βγG1', ylab = 'βγG2', cex = 0.5, main = ' βγG with RΔ contour')
contour(beta1, beta2, Rdelta, level = 0, xlab = "β1", ylab = "β2", add = TRUE)








###### THEORETICAL BETAGAMMA FROM RDELTA+ 
gamma<-seq(0,100,length.out=5000)
gamma<-c(gamma,100,1000,10000)

Zdelta <- matrix(c(1, -1), nrow = 2, byrow = TRUE)
Gdelta <- matrix(c(1, -1, -1, 0), nrow = 2, byrow = TRUE)
Zplus <- matrix(c(3, -5), ncol = 1)
Gplus <- matrix(c(3, -3, -3, 8), nrow = 2, byrow = TRUE)


eigen_Gdelta <- eigen(Gdelta)
V <- eigen_Gdelta$vectors
Lambda <- diag(abs(eigen_Gdelta$values)) 
G_delta_plus <- V %*% Lambda %*% t(V)


F_Delta_plus <- matrix(c(1, 0, 0,
                         0, 1, 0,
                         0, 0, 0), byrow = TRUE, nrow = 3)


beta1 <- seq(0, 2, length.out = 100)
beta2 <- seq(-1, 1, length.out = 100)
beta <- expand.grid(beta1 = beta1, beta2 = beta2)


Rdelta <- matrix(0, nrow = length(beta1), ncol = length(beta2))
R_delta_plus <- matrix(0, nrow = length(beta1), ncol = length(beta2))
R_delta_g <- matrix(0, nrow = length(beta1), ncol = length(beta2))


p <- 2
B <- matrix(c(0, 0, 0,
              0, 0, -1,
              1, 0, 0),
            nrow = p + 1, byrow = TRUE)
I <- diag(p + 1)
Px <- matrix(c(1, 0, 0,
               0, 1, 0),
             nrow = p, ncol = p + 1, byrow = TRUE)
Py <- matrix(c(0, 0, 1), nrow = 1, ncol = p + 1)
beta_cp <- matrix(c(1, 0), nrow = 1)

# Loop over grid and compute values
for (i in 1:length(beta1)) {
  for (j in 1:length(beta2)) {
    beta_vec <- c(beta1[i], beta2[j])
    term1 <- 2 * sum(Zdelta * beta_vec)
    term2 <- t(beta_vec) %*% Gdelta %*% beta_vec
    Rdelta[i, j] <- 1 - term1 + term2
    
    ub <- (Py - beta_vec %*% Px) %*% solve(I - B)
    R_delta_plus[i, j] <- ub %*% F_Delta_plus %*% t(ub)
    R_delta_g[i, j] <- (beta_vec - beta_cp) %*% G_delta_plus %*% t(beta_vec - beta_cp)
  }
}




# Compute beta_gamma
u_x <- Px %*% solve(I - B)
u_y <- Py %*% solve(I - B)
# Loop over all gamma values
# Preallocate a matrix to store beta_gamma for each gamma
beta_gamma_all <- matrix(0, nrow = 2, ncol = length(gamma))

for (g in 1:length(gamma)) {
  A <-  gamma[g] * u_x %*% F_Delta_plus %*% t(u_x)
  b <- Zplus +  gamma[g] * u_x %*% F_Delta_plus %*% t(u_y)
  beta_gamma_new <- solve(Gplus + A, b)
  beta_gamma_all[,g] <- beta_gamma_new[, 1]  
}

plot(beta_gamma_all[1,],beta_gamma_all[2,],type='l',xlab='β1',ylab='β2')
contour(beta1,beta2,Rdelta,level=0,add=TRUE)



R_delta_plus_gamma <- numeric(length(gamma))

IB_inv <- solve(I - B)

for (g in 1:length(gamma)) {
  beta_gamma <- beta_gamma_all[, g, drop = FALSE]  # 2x1 column vector
  u_beta <- (Py - t(beta_gamma) %*% Px) %*% IB_inv
  R_delta_plus_gamma[g] <- u_beta %*% F_Delta_plus %*% t(u_beta)
}

plot(gamma,R_delta_plus_gamma,type='l',xlab='γ',ylab='RΔ+',main='RΔ+ vs γ',cex=0.7,xlim=c(1,100))



### PLOTSS

# TOGETHER beta_gammaG and beta_gamma_min

plot(beta_gamma_G[1, ], beta_gamma_G[2, ],
     type = 'l', xlim = c(0.4, 1.2), ylim = c(-0.6, 1.1),
     xlab = 'βγ1', ylab = 'βγ2', cex = 0.7,
    # main = 'βγG, βγ- & βγ, Population Study',
    lwd=3)
contour(beta1, beta2, Rdelta, level = 0, add = TRUE)

lines(beta_min[1, ], beta_min[2, ],
      type = 'l', cex = 0.5, col = 'darkorange', pch = 8,lwd=3)
lines(beta_gamma_all[1,],beta_gamma_all[2,],type='l',xlab='β1',ylab='β2',col='blue',lty=22,lwd=3)


legend("topleft",
       legend = c("βγG", "βγ-", "βγ"),
       col = c("black", "darkorange","blue"),
       cex=0.9,
       lty=c(1,1,22),
       lwd=c(3,3,3),
       y.intersp = 0.6,  
       x.intersp = 0.5,
       bty = "o")
# |Rdelta| e Rdelta^+ uguali perchè matrice diag 

plot(gamma,R_delta_G,type='l',ylab='RΔG,|RΔ| & RΔ+',xlab='γ',lwd=3,xlim=c(0,30),cex=0.4,ylim=c(-0.05,0.25))

lines(gamma,abs(R_delta_min), type = 'l', cex = 0.5, col = 'darkorange', lwd= 3)
lines(gamma,R_delta_plus_gamma, type = 'l', cex = 0.5, col = 'blue', lty= 22,lwd=3)
legend("topright",
       legend = c("RΔG","|RΔ|", "RΔ+"),
       col = c("black", "darkorange", "blue"),
       lty = c(1, 1, 22),
       lwd = c(3, 3, 3),
       cex = 0.9,
       y.intersp = 0.6,
       x.intersp = 0.5,
       bty = "o")

