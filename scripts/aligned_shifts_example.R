#####################################
####### EMPIRICAL STUDY #############
####### ALIGNED SHIFTS ##############
#####################################

source("./datautil.R")
source("./cd.R")
source("./ols.R")
source("./uRdeltaRHS.R")
source("./measurements.R")


set.seed(123) 
# With a11 and a12 independent
X1<- rnorm(100000)
a11<-rnorm(10000,0,1)
X1<- X1+a11
ye<- X1+ rnorm(10000)
a12<-rnorm(10000,0,3)
X2<- -ye+rnorm(10000)+a12
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
               0, 1, 0),  # Just an example
             nrow = p, ncol = p+1, byrow = TRUE)
Py <- matrix(c(0, 0, 1), nrow = 1, ncol = p+1)


data_e <- list(Xe = Xe, ye = ye)



X1<- rnorm(10000)
yo<- X1+rnorm(10000)
a22<-rnorm(10000,0,1)
X2<- -yo+rnorm(10000)+a22 
Xo<- matrix(c(X1,X2),ncol=2)

betaXo<-lm(yo~Xo-1)$coef
betaXo


data_o <- list(Xo = Xo, yo = yo)

data_combined <- c(data_e, data_o)
m <- moments(data_combined)
beta_pa<-compute_cd(m)
beta_ols<-compute_ols(m)


gamma<-seq(-1,100,length.out=500)
gamma<-c(gamma,1000,10000,100000,1000000,1000000000)
n <- length(gamma)

dims <- dim(m$Gplus)

beta_gamma <- matrix(NA, nrow = dims[1], ncol = n)




for (i in seq_along(gamma)) {
  Gtemp <- m$Gplus + (gamma[i] * m$Gdelta)
  Ztemp <- m$Zplus + (gamma[i] * m$Zdelta)
  beta_gamma[, i] <- solve(Gtemp, Ztemp)
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

# PLOT RDELTA vs. GAMMA 

plot(gamma,Rdelta_vec,type='p',xlim=c(-1,100),xlab='γ ',ylab='R∆')
abline(h=0,col='red')



###############################################

grid_n <- 200
beta1_vals <- seq(-2, 2,length.out = grid_n)
beta2_vals <- seq(-2, 2, length.out = grid_n)

grid <- expand.grid(beta1 = beta1_vals, beta2 = beta2_vals)
n_points <- nrow(grid)


Rplus <- numeric(n_points)
Rdelta <- numeric(n_points)
for (i in 1:n_points) {
  current_param <- matrix(c(grid$beta1[i], grid$beta2[i]), nrow = 2, ncol = 1)
  Rplus[i] <- in_sample_risk(data_combined, current_param)
  Rdelta[i] <- difference(data_combined, current_param)
}


# 3D PLOT, Rdelta and Rplus

library(plotly)
library(RColorBrewer)
palette<-brewer.pal(9,"Dark2")
palette2<-brewer.pal(7,"BrBG")
display.brewer.pal(7,"Accent")
# Create matrices for surfaces
Rdelta_matrix <- matrix(Rdelta, nrow = grid_n, ncol = grid_n, byrow = TRUE)
Rplus_matrix  <- matrix(Rplus, nrow = grid_n, ncol = grid_n, byrow = TRUE)
p <- plot_ly()  # Add the Rdelta surface with a yellow-to-red color scale 
p <- p %>% add_surface(   x = ~beta1_vals,   y = ~beta2_vals,   z = ~Rdelta_matrix,   colorscale = list(c(0, "yellow"), c(1, palette[2])),   showscale = TRUE,   colorbar = list(title = "RΔ ", tickfont = list(color =  palette[2]),                   titlefont = list(color =  palette[2])),   contours = list(     x = list(show = TRUE, color = "black"),     y = list(show = TRUE, color = "black"),     z = list(show = TRUE, color = "black")   ) ) 


# Add the Rplus surface with a blue-to-green color scale
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

Rplus_beta_pa<-in_sample_risk(data_combined,beta_pa) # 1.923168
p <- p %>% add_markers(
  x = beta_pa[1],
  y = beta_pa[2],
  z = Rplus_beta_pa,
  marker = list(size = 8, color = "darkblue"),
  name = paste("βCP ")
)
Rplus_beta_ols<-in_sample_risk(data_combined, beta_ols)

p <- p %>% add_markers(
  x = beta_ols[1],
  y = beta_ols[2],
  z = Rplus_beta_ols,
  marker = list(size = 8, color = "red"),
  name = paste("βOLS")
)


p <- p %>% layout(
  title = list(
    text = "Surface Plot of RΔ and R+",
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
    x = 0.2,     # Move it further right (adjust as needed)
    y = 0.2,      # Lower it (0 = bottom, 1 = top)
    len = 0.5     # Shorten the colorbar
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

# beta_ols is the global minimum for R_+:
# beta_pa is the global minimum of R_+ when R_Delta=0.












