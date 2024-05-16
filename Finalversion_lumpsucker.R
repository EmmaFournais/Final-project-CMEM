#Emma Fournais Næsheim and Rebecca Martinovic 
#Final project for computational marine ecological modelling 
#16-05-2023



library(deSolve)
library(fields)
library(ggplot2)
library(RColorBrewer)

#Parameters 
lb <- 7 #mm 
lm <- 500 #mm 
lj <- 320 #mm 
lv <- 50 #mm 
Im <- 2.79*10^-4 #g/day/mm^2 Propotionality constant (ressource ingestion of c)
Rh <- 1.5*10^-5 #g/liter Half saturation constant 
rm <- 0.0018 #/day/mm2 Proportionality constant when reached maturity
gamma <- 0.006 #/day 
u <- 0.01 #day consumer mortality
p <- 0.1 #day flow-through rate
K <- 0.0003 #g/liter max ressource biomass 
e <- 0.5 #conversion efficiency 
a <- 5000 #liter/day Attack rate 
Th <- 0.1 #day/g handling time
sigma <- 0.01 #day background mortality rate for predator when foragin on lb an
dF <- 0.01##/Day
beta <- 9.0 * 10^-6 #g/mm^-3 biomass of consumers proportionality constant 

n <- 100
deltal <- lm/n

l <- seq(deltal, lm, by = deltal)
R <- K
P <- 0.0002
C <- rep(0,n)
C[1] <- 10^-4


#-----------no fishing
derivative <- function(time, y, params) {
  C <- y[1:n]
  R <- y[n+1]
  P <- y[n+2]
  
  g <- pmax(0,gamma*(lm*R/(Rh+R)-l))
  
  I <- (Im*l^2*R)/(Rh+R)
  
  b <- (rm*l^2*R)/(Rh+R)
  
  Biomassconsumers <- function(beta,l,t,C,deltal,lv){
    ixVulnerable = l < lv
    B <- sum(beta*l[ixVulnerable]^3*C[ixVulnerable]*deltal) 
    return(B)
  }
  B <- Biomassconsumers(beta,l,t,C,deltal,lv)
  
  dp = rep(0,n)
  ixVulnerable = l < lv
  dp[ixVulnerable] <- (a*P)/(1+a*Th*B)
  
  
  # Initializing advective
  fluxC <- numeric(n+1)
  
  # Calculate advective fluxes (Ja)
  ixMature = l>lj
  ixyoung = l<lj
  fluxC[1] <-  sum(b[ixMature] * C[ixMature] * deltal) # Boundary condition
  ixall = l > lb
  
  for(i in 2:n){
    fluxC[i] <- g[i-1] * C[i-1] 
  } 
  # Set boundary fluxes to 0
  fluxC[n+1] <- 0
  
  dC_dt <- numeric(n)
 
  for (i in 1:n) {
    dC_dt[i] <-  (-(u + dp[i]) * C[i]) + ((fluxC[i] - fluxC[i+1]) / deltal) 
  }
  
  dR_dt <- p*(K-R) - (sum(I[ixall] * C[ixall] * deltal)) 
  dP_dt <- ((e*((a*B)/(1+a*Th*B))-sigma))*P 
  
  return(list(c(dC_dt,dR_dt,dP_dt)))
}

params <-
  list(
    deltal = deltal,
    u = u,
    e = e,
    a = a,
    Th = Th,
    sigma = sigma,
    beta = beta,
    p = p,
    gamma = gamma,
    lm = lm, 
    lv = lv, 
    lj = lj, 
    lb = lb, 
    K = K, 
    rm = rm 
  
  )

# Time
time <- seq(0, 10000, by = 1) # Time vector

# ODE solving
resultNo <- ode(y = c(C, R, P), time = time, func = derivative, parms = params)

C_outNo<-resultNo[,2:(n+1)]
R_outNo <-resultNo[, n+2]
P_outNo <-resultNo[, n+3]

C_outNo[C_outNo<1e-25] = 1e-25
par(mfrow=c(1,2))

image.plot(
  time,
  l,
  log(C_outNo),
  col = hcl.colors(50, "viridis"),
  main = "Consumer concentration [g/liter] over length and time - No fishing ",
  xlab = "Days",
  ylab = "Size [mm]",
  cex.lab = 1.25
)

  #----------first d(F)
derivative <- function(time, y, params) {
  C <- y[1:n]
  R <- y[n+1]
  P <- y[n+2]
  
  g <- pmax(0,gamma*(lm*R/(Rh+R)-l))
  
  I <- (Im*l^2*R)/(Rh+R)
  
  b <- (rm*l^2*R)/(Rh+R)
  
  Biomassconsumers <- function(beta,l,t,C,deltal,lv){
    ixVulnerable = l < lv
    B <- sum(beta*l[ixVulnerable]^3*C[ixVulnerable]*deltal) 
    return(B)
  }
  B <- Biomassconsumers(beta,l,t,C,deltal,lv)
  
  dp = rep(0,n)
  ixVulnerable = l < lv
  dp[ixVulnerable] <- (a*P)/(1+a*Th*B)
  
  
  # Initializing advective
  fluxC <- numeric(n+1)
  
  # Calculate advective fluxes (Ja)
  ixMature = l>lj
  ixyoung = l<lj
  fluxC[1] <-  sum(b[ixMature] * C[ixMature] * deltal) # Boundary condition
  ixall = l > lb
  dF = rep(0,n)
  dF[ixMature] <- 0.01
  dF[ixyoung] <- 0.0001
 
  
  for(i in 2:n){
    fluxC[i] <- g[i-1] * C[i-1] 
  }
  
  # Set boundary fluxes to 0
  fluxC[n+1] <- 0
  
  dC_dt <- numeric(n)
  
  for (i in 1:n) {
    dC_dt[i] <-  (-(u + dF[i] + dp[i]) * C[i]) + ((fluxC[i] - fluxC[i+1]) / deltal) # Flux[i] 
  }
  
  dR_dt <- p*(K-R) - (sum(I[ixall] * C[ixall] * deltal)) 
  dP_dt <- ((e*((a*B)/(1+a*Th*B))-sigma))*P 
  
  return(list(c(dC_dt,dR_dt,dP_dt)))
}

params <-
  list(
    deltal = deltal,
    u = u,
    e = e,
    a = a,
    Th = Th,
    sigma = sigma,
    beta = beta,
    p = p,
    gamma = gamma,
    lm = lm, 
    lv = lv, 
    lj = lj, 
    lb = lb, 
    K = K, 
    rm = rm 
  )

# Time
time <- seq(0, 10000, by = 1) # Time vector

# ODE solving
resultN <- ode(y = c(C, R, P), time = time, func = derivative, parms = params)

C_outN<-resultN[,2:(n+1)]
R_outN <-resultN[, n+2]
P_outN <-resultN[, n+3]

C_outN[C_outN<1e-25] = 1e-25
image.plot(
  time,
  l,
  log(C_outN),
  col = hcl.colors(50, "viridis"),
  main = "Consumer concentration [g/liter] over length and time - d(F) = 0.01",
  xlab = "Days",
  ylab = "Size [mm]",
  cex.lab = 1.25
)

#---------------------------------------------------

derivative <- function(time, y, params) {
  C <- y[1:n]
  R <- y[n+1]
  P <- y[n+2]
  
  g <- pmax(0,gamma*(lm*R/(Rh+R)-l))
  
  I <- (Im*l^2*R)/(Rh+R)
  
  b <- (rm*l^2*R)/(Rh+R)
  
  Biomassconsumers <- function(beta,l,t,C,deltal,lv){
    ixVulnerable = l < lv
    B <- sum(beta*l[ixVulnerable]^3*C[ixVulnerable]*deltal) 
    return(B)
  }
  B <- Biomassconsumers(beta,l,t,C,deltal,lv)
  
  dp = rep(0,n)
  ixVulnerable = l < lv
  dp[ixVulnerable] <- (a*P)/(1+a*Th*B)
  
  
  # Initializing advective
  fluxC <- numeric(n+1)
  
  # Calculate advective fluxes (Ja)
  ixMature = l>lj
  ixyoung = l<lj
  fluxC[1] <-  sum(b[ixMature] * C[ixMature] * deltal) # Boundary condition
  ixall = l > lb
  dF = rep(0,n)
  dF[ixMature] <- 0.02
  dF[ixyoung] <- 0.0001
  
  for(i in 2:n){
    fluxC[i] <- g[i-1] * C[i-1] 
  }
  
  # Set boundary fluxes to 0
  fluxC[n+1] <- 0
  
  dC_dt <- numeric(n)
  
  for (i in 1:n) {
    dC_dt[i] <-  (-(u + dF[i] + dp[i]) * C[i]) + ((fluxC[i] - fluxC[i+1]) / deltal) # Flux[i] 
  }
  
  dR_dt <- p*(K-R) - (sum(I[ixall] * C[ixall] * deltal)) 
  dP_dt <- ((e*((a*B)/(1+a*Th*B))-sigma))*P 
  
  return(list(c(dC_dt,dR_dt,dP_dt)))
}

params <-
  list(
    deltal = deltal,
    u = u,
    e = e,
    a = a,
    Th = Th,
    sigma = sigma,
    beta = beta,
    p = p,
    gamma = gamma,
    lm = lm, 
    lv = lv, 
    lj = lj, 
    lb = lb, 
    K = K, 
    rm = rm 
  )

# Time
time <- seq(0, 10000, by = 1) # Time vector

# ODE solving
result <- ode(y = c(C, R, P), time = time, func = derivative, parms = params)

C_out<-result[,2:(n+1)]
R_out <-result[, n+2]
P_out <-result[, n+3]

C_out[C_out<1e-25] = 1e-25
image.plot(time, l,log(C_out), col = hcl.colors(50, "viridis"), main = "Consumer concentration over length and time - dF = 0.01", xlab = "Days", ylab = "Lenght", cex.lab=1.25)

#---------------------------------------------
derivative <- function(time, y, params) {
  C <- y[1:n]
  R <- y[n+1]
  P <- y[n+2]
  
  g <- pmax(0,gamma*(lm*R/(Rh+R)-l))
  
  I <- (Im*l^2*R)/(Rh+R)
  
  b <- (rm*l^2*R)/(Rh+R)
  
  Biomassconsumers <- function(beta,l,t,C,deltal,lv){
    ixVulnerable = l < lv
    B <- sum(beta*l[ixVulnerable]^3*C[ixVulnerable]*deltal) 
    return(B)
  }
  B <- Biomassconsumers(beta,l,t,C,deltal,lv)
  
  dp = rep(0,n)
  ixVulnerable = l < lv
  dp[ixVulnerable] <- (a*P)/(1+a*Th*B)
  
  
  # Initializing advective
  fluxC <- numeric(n+1)
  
  # Calculate advective fluxes (Ja)
  ixMature = l>lj
  ixyoung = l<lj
  fluxC[1] <-  sum(b[ixMature] * C[ixMature] * deltal) # Boundary condition
  ixall = l > lb
  dF1 = rep(0,n)
  dF1[ixMature] <- 0.03
  dF1[ixyoung] <- 0.0001
  
  for(i in 2:n){
    fluxC[i] <- g[i-1] * C[i-1] 
  }
  
  # Set boundary fluxes to 0
  fluxC[n+1] <- 0
  
  
  dC_dt <- numeric(n)
  
  for (i in 1:n) {
    dC_dt[i] <-  (-(u + dF1[i] + dp[i]) * C[i]) + ((fluxC[i] - fluxC[i+1]) / deltal) # Flux[i] 
  }
  
  dR_dt <- p*(K-R) - (sum(I[ixall] * C[ixall] * deltal)) 
  dP_dt <- ((e*((a*B)/(1+a*Th*B))-sigma))*P 
  
  return(list(c(dC_dt,dR_dt,dP_dt)))
}

params <-
  list(
    deltal = deltal,
    u = u,
    e = e,
    a = a,
    Th = Th,
    sigma = sigma,
    beta = beta,
    p = p,
    gamma = gamma,
    lm = lm, 
    lv = lv, 
    lj = lj, 
    lb = lb, 
    K = K, 
    rm = rm 
  )


# ODE solving
result1 <- ode(y = c(C, R, P), time = time, func = derivative, parms = params)
C_out1<-result1[,2:(n+1)]
R_out1 <-result1[, n+2]
P_out1 <-result1[, n+3]

C_out1[C_out1<1e-25] = 1e-25

image.plot(
  time,
  l,
  log(C_out1),
  col = hcl.colors(50, "viridis"),
  main = "Consumer concentration [g/liter] over length and time - d(F) = 0.03",
  xlab = "Days",
  ylab = "Size [mm]",
  cex.lab = 1.25
)


#-------------------------------------------------------------------
derivative <- function(time, y, params) {
  C <- y[1:n]
  R <- y[n+1]
  P <- y[n+2]
  
  g <- pmax(0,gamma*(lm*R/(Rh+R)-l))
  
  I <- (Im*l^2*R)/(Rh+R)
  
  b <- (rm*l^2*R)/(Rh+R)
  
  Biomassconsumers <- function(beta,l,t,C,deltal,lv){
    ixVulnerable = l < lv
    B <- sum(beta*l[ixVulnerable]^3*C[ixVulnerable]*deltal) 
    return(B)
  }
  B <- Biomassconsumers(beta,l,t,C,deltal,lv)
  
  dp = rep(0,n)
  ixVulnerable = l < lv
  dp[ixVulnerable] <- (a*P)/(1+a*Th*B)
  
  
  
  # Initializing advective
  fluxC <- numeric(n+1)
  
  # Calculate advective fluxes (Ja)
  ixMature = l>lj
  ixyoung = l<lj
  fluxC[1] <-  sum(b[ixMature] * C[ixMature] * deltal) # Boundary condition
  ixall = l > lb
  dF2 = rep(0,n)
  dF2[ixMature] <- 0.04
  dF2[ixyoung] <- 0.0001
  
  
  for(i in 2:n){
    fluxC[i] <- g[i-1] * C[i-1] 
  }
  
  # Set boundary fluxes to 0
  fluxC[n+1] <- 0
  
  dC_dt <- numeric(n)
  
  for (i in 1:n) {
    dC_dt[i] <-  (-(u + dF2[i] + dp[i]) * C[i]) + ((fluxC[i] - fluxC[i+1]) / deltal) # Flux[i] 
  }
  
  dR_dt <- p*(K-R) - (sum(I[ixall] * C[ixall] * deltal)) 
  dP_dt <- ((e*((a*B)/(1+a*Th*B))-sigma))*P 
  
  return(list(c(dC_dt,dR_dt,dP_dt)))
}

params <-
  list(
    deltal = deltal,
    u = u,
    e = e,
    a = a,
    Th = Th,
    sigma = sigma,
    beta = beta,
    p = p,
    gamma = gamma,
    lm = lm, 
    lv = lv, 
    lj = lj, 
    lb = lb, 
    K = K, 
    rm = rm 
  )

# ODE solving
result2 <- ode(y = c(C, R, P), time = time, func = derivative, parms = params)
C_out2<-result2[,2:(n+1)]
R_out2 <-result2[, n+2]
P_out2 <-result2[, n+3]

C_out2[C_out2<1e-25] = 1e-25
image.plot(
  time,
  l,
  log(C_out2),
  col = hcl.colors(50, "viridis"),
  main = "Consumer concentration [g/liter] over length and time - d(F) = 0.04",
  xlab = "Days",
  ylab = "Size [mm]",
  cex.lab = 1.25
)

#--------------------------------------------------------------------------
derivative <- function(time, y, params) {
  C <- y[1:n]
  R <- y[n+1]
  P <- y[n+2]
  
  g <- pmax(0,gamma*(lm*R/(Rh+R)-l))
  
  I <- (Im*l^2*R)/(Rh+R)
  
  b <- (rm*l^2*R)/(Rh+R)
  
  Biomassconsumers <- function(beta,l,t,C,deltal,lv){
    ixVulnerable = l < lv
    B <- sum(beta*l[ixVulnerable]^3*C[ixVulnerable]*deltal) 
    return(B)
  }
  B <- Biomassconsumers(beta,l,t,C,deltal,lv)
  
  dp = rep(0,n)
  ixVulnerable = l < lv
  dp[ixVulnerable] <- (a*P)/(1+a*Th*B)
  
  
  
  # Initializing advective
  fluxC <- numeric(n+1)
  
  # Calculate advective fluxes (Ja)
  ixMature = l>lj
  ixyoung = l<lj
  fluxC[1] <-  sum(b[ixMature] * C[ixMature] * deltal) # Boundary condition
  ixall = l > lb
  dF3 = rep(0,n)
  dF3[ixMature] <- 0.05
  dF3[ixyoung] <- 0.0001
  
  for(i in 2:n){
    fluxC[i] <- g[i-1] * C[i-1]
  }
  
  # Set boundary fluxes to 0
  fluxC[n+1] <- 0
  
  dC_dt <- numeric(n)
  
  for (i in 1:n) {
    dC_dt[i] <-  (-(u + dF3[i] + dp[i]) * C[i]) + ((fluxC[i] - fluxC[i+1]) / deltal) # Flux[i] 
  }
  
  dR_dt <- p*(K-R) - (sum(I[ixall] * C[ixall] * deltal)) 
  dP_dt <- ((e*((a*B)/(1+a*Th*B))-sigma))*P 
  
  return(list(c(dC_dt,dR_dt,dP_dt)))
}

params <-
  list(
    deltal = deltal,
    u = u,
    e = e,
    a = a,
    Th = Th,
    sigma = sigma,
    beta = beta,
    p = p,
    gamma = gamma,
    lm = lm, 
    lv = lv, 
    lj = lj, 
    lb = lb, 
    K = K, 
    rm = rm 
  )


# ODE solving
result3 <- ode(y = c(C, R, P), time = time, func = derivative, parms = params)
C_out3<-result3[,2:(n+1)]
R_out3 <-result3[, n+2]
P_out3 <-result3[, n+3]

C_out3[C_out3<1e-25] = 1e-25
image.plot(time, l,log(C_out3), col = hcl.colors(50, "viridis"), main = "Consumer concentration over length and time - dF3", xlab = "time", ylab = "Lenght", cex.lab=1.25)

#-----------------------------------------------------------------
derivative <- function(time, y, params) {
  C <- y[1:n]
  R <- y[n+1]
  P <- y[n+2]
  
  g <- pmax(0,gamma*(lm*R/(Rh+R)-l))
  
  I <- (Im*l^2*R)/(Rh+R)
  
  b <- (rm*l^2*R)/(Rh+R)
  
  Biomassconsumers <- function(beta,l,t,C,deltal,lv){
    ixVulnerable = l < lv
    B <- sum(beta*l[ixVulnerable]^3*C[ixVulnerable]*deltal) 
    return(B)
  }
  B <- Biomassconsumers(beta,l,t,C,deltal,lv)
  
  dp = rep(0,n)
  ixVulnerable = l < lv
  dp[ixVulnerable] <- (a*P)/(1+a*Th*B)
  
  
  
  # Initializing advective
  fluxC <- numeric(n+1)
  
  # Calculate advective fluxes (Ja)
  ixMature = l>lj
  ixyoung = l<lj
  fluxC[1] <-  sum(b[ixMature] * C[ixMature] * deltal) # Boundary condition
  ixall = l > lb
  dF4 = rep(0,n)
  dF4[ixMature] <- 0.06
  dF4[ixyoung] <- 0.0001
  
  for(i in 2:n){
    fluxC[i] <- g[i-1] * C[i-1]
  }
  
  # Set boundary fluxes to 0
  fluxC[n+1] <- 0
  
  dC_dt <- numeric(n)
  
  
  for (i in 1:n) {
    dC_dt[i] <-  (-(u + dF4[i] + dp[i]) * C[i]) + ((fluxC[i] - fluxC[i+1]) / deltal) # Flux[i] 
  }
  
  dR_dt <- p*(K-R) - (sum(I[ixall] * C[ixall] * deltal)) 
  dP_dt <- ((e*((a*B)/(1+a*Th*B))-sigma))*P 
  
  return(list(c(dC_dt,dR_dt,dP_dt)))
}

params <-
  list(
    deltal = deltal,
    u = u,
    e = e,
    a = a,
    Th = Th,
    sigma = sigma,
    beta = beta,
    p = p,
    gamma = gamma,
    lm = lm, 
    lv = lv, 
    lj = lj, 
    lb = lb, 
    K = K, 
    rm = rm 
  )

# ODE solving
result4 <- ode(y = c(C, R, P), time = time, func = derivative, parms = params)
C_out4<-result4[,2:(n+1)]
R_out4 <-result4[, n+2]
P_out4 <-result4[, n+3]

C_out4[C_out4<1e-25] = 1e-25
image.plot(time, l,log(C_out4), col = hcl.colors(50, "viridis"), main = "Consumer concentration over length and time - dF4", xlab = "time", ylab = "Lenght", cex.lab=1.25)

#-----------------------------------------------------------------------
derivative <- function(time, y, params) {
  C <- y[1:n]
  R <- y[n+1]
  P <- y[n+2]
  
  g <- pmax(0,gamma*(lm*R/(Rh+R)-l))
  
  I <- (Im*l^2*R)/(Rh+R)
  
  b <- (rm*l^2*R)/(Rh+R)
  
  Biomassconsumers <- function(beta,l,t,C,deltal,lv){
    ixVulnerable = l < lv
    B <- sum(beta*l[ixVulnerable]^3*C[ixVulnerable]*deltal) 
    return(B)
  }
  B <- Biomassconsumers(beta,l,t,C,deltal,lv)
  
  dp = rep(0,n)
  ixVulnerable = l < lv
  dp[ixVulnerable] <- (a*P)/(1+a*Th*B)
  
  
  # Initializing advective
  fluxC <- numeric(n+1)
  
  # Calculate advective fluxes (Ja)
  ixMature = l>lj
  ixyoung = l<lj
  fluxC[1] <-  sum(b[ixMature] * C[ixMature] * deltal) # Boundary condition
  ixall = l > lb
  dF5 = rep(0,n)
  dF5[ixMature] <- 0.07
  dF5[ixyoung] <- 0.0001
  
  for(i in 2:n){
    fluxC[i] <- g[i-1] * C[i-1] 
  }
  
  # Set boundary fluxes to 0
  fluxC[n+1] <- 0
  
  dC_dt <- numeric(n)
  
  for (i in 1:n) {
    dC_dt[i] <-  (-(u + dF5[i] + dp[i]) * C[i]) + ((fluxC[i] - fluxC[i+1]) / deltal) # Flux[i] 
  }
  
  dR_dt <- p*(K-R) - (sum(I[ixall] * C[ixall] * deltal)) 
  dP_dt <- ((e*((a*B)/(1+a*Th*B))-sigma))*P 
  
  return(list(c(dC_dt,dR_dt,dP_dt)))
}

params <-
  list(
    deltal = deltal,
    u = u,
    e = e,
    a = a,
    Th = Th,
    sigma = sigma,
    beta = beta,
    p = p,
    gamma = gamma,
    lm = lm, 
    lv = lv, 
    lj = lj, 
    lb = lb, 
    K = K, 
    rm = rm 
  )

# ODE solving
result5 <- ode(y = c(C, R, P), time = time, func = derivative, parms = params)
C_out5<-result5[,2:(n+1)]
R_out5 <-result5[, n+2]
P_out5 <-result5[, n+3]

C_out5[C_out5<1e-25] = 1e-25
image.plot(time, l,log(C_out5), col = hcl.colors(50, "viridis"), main = "Consumer concentration over length and time - dF5", xlab = "time", ylab = "Lenght", cex.lab=1.25)

#-----------------------------------------------------------------------------

derivative <- function(time, y, params) {
  C <- y[1:n]
  R <- y[n+1]
  P <- y[n+2]
  
  g <- pmax(0,gamma*(lm*R/(Rh+R)-l))
  
  I <- (Im*l^2*R)/(Rh+R)
  
  b <- (rm*l^2*R)/(Rh+R)
  
  Biomassconsumers <- function(beta,l,t,C,deltal,lv){
    ixVulnerable = l < lv
    B <- sum(beta*l[ixVulnerable]^3*C[ixVulnerable]*deltal) 
    return(B)
  }
  B <- Biomassconsumers(beta,l,t,C,deltal,lv)
  
  dp = rep(0,n)
  ixVulnerable = l < lv
  dp[ixVulnerable] <- (a*P)/(1+a*Th*B)
  
  
  # Initializing advective
  fluxC <- numeric(n+1)
  
  # Calculate advective fluxes (Ja)
  ixMature = l>lj
  ixyoung = l<lj
  fluxC[1] <-  sum(b[ixMature] * C[ixMature] * deltal) # Boundary condition
  ixall = l > lb
  dF6 = rep(0,n)
  dF6[ixMature] <- 0.08
  dF6[ixyoung] <- 0.0001
  
  for(i in 2:n){
    fluxC[i] <- g[i-1] * C[i-1]
  }
  
  # Set boundary fluxes to 0
  fluxC[n+1] <- 0
  
  dC_dt <- numeric(n)
  
  for (i in 1:n) {
    dC_dt[i] <-  (-(u + dF6[i] + dp[i]) * C[i]) + ((fluxC[i] - fluxC[i+1]) / deltal) # Flux[i] 
  }
  
  dR_dt <- p*(K-R) - (sum(I[ixall] * C[ixall] * deltal)) 
  dP_dt <- ((e*((a*B)/(1+a*Th*B))-sigma))*P 
  
  return(list(c(dC_dt,dR_dt,dP_dt)))
}

params <-
  list(
    deltal = deltal,
    u = u,
    e = e,
    a = a,
    Th = Th,
    sigma = sigma,
    beta = beta,
    p = p,
    gamma = gamma,
    lm = lm, 
    lv = lv, 
    lj = lj, 
    lb = lb, 
    K = K, 
    rm = rm 
  )


# ODE solving
result6 <- ode(y = c(C, R, P), time = time, func = derivative, parms = params)
C_out6<-result6[,2:(n+1)]
R_out6 <-result6[, n+2]
P_out6 <-result6[, n+3]

C_out6[C_out6<1e-25] = 1e-25
image.plot(time, l,log(C_out6), col = hcl.colors(50, "viridis"), main = "Consumer concentration over length and time - dF = 0.8", xlab = "Days", ylab = "Lenght", cex.lab=1.25)

#---------------------------------------------------------------------------------

derivative <- function(time, y, params) {
  C <- y[1:n]
  R <- y[n+1]
  P <- y[n+2]
  
  g <- pmax(0,gamma*(lm*R/(Rh+R)-l))
  
  I <- (Im*l^2*R)/(Rh+R)
  
  b <- (rm*l^2*R)/(Rh+R)
  
  Biomassconsumers <- function(beta,l,t,C,deltal,lv){
    ixVulnerable = l < lv
    B <- sum(beta*l[ixVulnerable]^3*C[ixVulnerable]*deltal) 
    return(B)
  }
  B <- Biomassconsumers(beta,l,t,C,deltal,lv)
  
  dp = rep(0,n)
  ixVulnerable = l < lv
  dp[ixVulnerable] <- (a*P)/(1+a*Th*B)
  
  
  
  # Initializing advective
  fluxC <- numeric(n+1)
  
  # Calculate advective fluxes (Ja)
  ixMature = l>lj
  ixyoung = l<lj
  fluxC[1] <-  sum(b[ixMature] * C[ixMature] * deltal) # Boundary condition
  ixall = l > lb
  dF7 = rep(0,n)
  dF7[ixMature] <- 0.09
  dF7[ixyoung] <- 0.0001
  
  for(i in 2:n){
    fluxC[i] <- g[i-1] * C[i-1]
  }
  
  # Set boundary fluxes to 0
  fluxC[n+1] <- 0
  
  dC_dt <- numeric(n)
  
  for (i in 1:n) {
    dC_dt[i] <-  (-(u + dF7[i] + dp[i]) * C[i]) + ((fluxC[i] - fluxC[i+1]) / deltal) # Flux[i] 
  }
  
  dR_dt <- p*(K-R) - (sum(I[ixall] * C[ixall] * deltal)) 
  dP_dt <- ((e*((a*B)/(1+a*Th*B))-sigma))*P 
  
  return(list(c(dC_dt,dR_dt,dP_dt)))
}

params <-
  list(
    deltal = deltal,
    u = u,
    e = e,
    a = a,
    Th = Th,
    sigma = sigma,
    beta = beta,
    p = p,
    gamma = gamma,
    lm = lm, 
    lv = lv, 
    lj = lj, 
    lb = lb, 
    K = K, 
    rm = rm 
  )


# ODE solving
result7 <- ode(y = c(C, R, P), time = time, func = derivative, parms = params)
C_out7<-result7[,2:(n+1)]
R_out7 <-result7[, n+2]
P_out7 <-result7[, n+3]

C_out7[C_out7<1e-25] = 1e-25
image.plot(time, l,log(C_out7), col = hcl.colors(50, "viridis"), main = "Consumer concentration [g/liter] over length and time - d(F) = 0.09", xlab = "Days", ylab = "Size [mm]", cex.lab=1.25)


#Calculating the mean of the last time steps 
CmeanNo <- C_outNo[8001:10000,]
CmeanNo <- colMeans(CmeanNo)
CmeanN <- C_outN[8001:10000,]
CmeanN <- colMeans(CmeanN)
Cmean <- C_out[8001:10000,]
Cmean <- colMeans(Cmean)
Cmean1 <- C_out1[8001:10000,]
Cmean1 <- colMeans(Cmean1)
Cmean2 <- C_out2[8001:10000,]
Cmean2 <- colMeans(Cmean2)
Cmean3 <- C_out3[8001:10000,]
Cmean3 <- colMeans(Cmean3)
Cmean4 <- C_out4[8001:10000,]
Cmean4 <- colMeans(Cmean4)
Cmean5 <- C_out5[8001:10000,]
Cmean5 <- colMeans(Cmean5)
Cmean6 <- C_out6[8001:10000,]
Cmean6 <- colMeans(Cmean6)
Cmean7 <- C_out7[8001:10000,]
Cmean7 <- colMeans(Cmean7)

Blue <- colorRampPalette(brewer.pal(5, "Blues"))(length(x))
Green <- colorRampPalette(brewer.pal(5, "Greens"))(length(x))
colors <- colorRampPalette(brewer.pal(5, "Reds"))(length(x))

plot(
  l,
  CmeanNo,
  type = "l",
  ylim = c(0, 0.0000002),
  lwd = 3,
  xlab = "Size [mm]",
  ylab = "Biomass [g/liter]",
  main = "Mean distribution over the last 2000 timesteps",
  col = Blue[1],
  cex.lab = 1.2,
  cex.axis = 1.2
)
lines(l, CmeanN, col = Blue[2], lwd = 3)
lines(l, Cmean, col = Blue[3], lwd = 3)
lines(l, Cmean1, col = Blue[4], lwd = 3)  # Use colors from the 'colors' vector
lines(l, Cmean2, col = Blue[5], lwd = 3)
lines(l, Cmean3, col = Blue[6], lwd = 3)
lines(l, Cmean4, col = Blue[7], lwd = 3)
lines(l, Cmean5, col = Blue[8], lwd = 3)
lines(l, Cmean6, col = Blue[9], lwd = 3)
lines(l, Cmean7, col = Blue[10], lwd = 3)

legend(
  "topright",
  legend = c("0","0.01","0.02","0.03","0.04", "0.05", "0.06", "0.07","0.08", "0.09"),
  lty = c(1),
  lwd = 3,
  col = Blue,
  title = "d(F): [1/day]"
)

#Concentration of last 2000 timesteps for P, R and C
PlastN <- tail(P_outN[8000:10000],1)
Plast <- tail(P_out[8000:10000],1)
Plast1 <- tail(P_out1[8000:10000],1)
Plast2 <- tail(P_out2[8000:10000],1)
Plast3 <- tail(P_out3[8000:10000],1)
Plast4 <- tail(P_out4[8000:10000],1)
Plast5 <- tail(P_out5[8000:10000],1)
Plast6 <- tail(P_out6[8000:10000],1)
Plast7 <- tail(P_out7[8000:10000],1)

RlastN <- tail(R_outN[8000:10000],1)
Rlast <- tail(R_out[8000:10000],1)
Rlast1 <- tail(R_out1[8000:10000],1)
Rlast2 <- tail(R_out2[8000:10000],1)
Rlast3 <- tail(R_out3[8000:10000],1)
Rlast4 <- tail(R_out4[8000:10000],1)
Rlast5 <- tail(R_out5[8000:10000],1)
Rlast6 <- tail(R_out6[8000:10000],1)
Rlast7 <- tail(R_out7[8000:10000],1)

ClastN <- sum(CmeanN)
Clast <- sum(Cmean)
Clast1 <- sum(Cmean1)
Clast2 <- sum(Cmean2)
Clast3 <- sum(Cmean3)
Clast4 <- sum(Cmean4)
Clast5 <- sum(Cmean5)
Clast6 <- sum(Cmean6)
Clast7 <- sum(Cmean7)


x <- c(ClastN,Clast, Clast1, Clast2, Clast3, Clast4, Clast5, Clast6, Clast7)
yvalues <- c(0.01,0.02,0.03,0.04, 0.05, 0.06, 0.07,0.08, 0.09)

plot(
  yvalues,
  x,
  lwd = "8",
  col = Blue,
  type = "p",
  xlab = "Fishing intensity [1/day]",
  ylab = "Biomass [g/liter]",
  main = "Overall Concentration of Consumers at different d(F) ",
  cex.lab = 1.2,
  cex.axis = 1.2
)
legend(
  "topright",
  legend = yvalues,
  lty = c(1),
  lwd = 3,
  col = Blue,
  title = "d(F): [1/day]"
)



# Assuming dF is your vector

xP <- c(PlastN, Plast, Plast1, Plast2, Plast3, Plast4, Plast5, Plast6, Plast7)
plot(
  yvalues,
  xP,
  lwd = "8",
  col = colors,
  type = "p",
  xlab = "Fishing intensity [1/day]",
  ylab = "Biomass [individuals/liter]",
  main = "Overall concentration of Predators at different d(F)",
  cex.lab = 1.2,
  cex.axis = 1.2
)
legend(
  "topright",
  legend = yvalues,
  lty = c(1),
  lwd = 3,
  col = colors,
  title = "d(F): [1/day]"
)

xR <- c(RlastN, Rlast, Rlast1, Rlast2, Rlast3, Rlast4, Rlast5, Rlast6, Rlast7)
plot(
  yvalues,
  xR,
  lwd = "8",
  col = Green,
  type = "p",
  xlab = "Fishing intensity [1/day]",
  ylab = "Biomass [g/liter]",
  main = "Overall concentration of Ressource at different d(F)",
  cex.lab = 1.2,
  cex.axis = 1.2
)

legend(
  "topright",
  legend = yvalues,
  lty = c(1),
  lwd = 3,
  col = Green,
  title = "d(F): [1/day]"
)


