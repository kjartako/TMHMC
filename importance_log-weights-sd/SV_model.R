rm(list=ls())
source("choleskys.R")
nrep <- 1000


# SV model

y <- as.vector(read.table("SV_data.txt")$V1)
T <- length(y)
gam <- -0.021
delta <- 0.98
v <- 0.14

gg <- function(y,x){
  return(dnorm(y,mean=0,sd=exp(0.5*x),log=TRUE))
}


mu1 <- gam/(1-delta)
sigma1 <- v/sqrt(1-delta^2)


# EIS

eisMean <- function(mu,sigma,a){
  return((mu + a[1]*sigma^2)/(1-2.0*a[2]*sigma^2))
}

eisStd <- function(sigma,a){
  return(sigma/sqrt(1.0-2.0*a[2]*sigma^2))
}

eislchi <- function(mu,sigma,a){
  t1 <- sigma ^ 2
  t2 <- t1 * a[2]
  t3 <- log(2)
  t6 <- log(pi)
  t9 <- log(sigma)
  t13 <- a[1] ^ 2
  t15 <- mu ^ 2
  t23 <- 2 * t2 - 1
  t27 <- exp(-1 / t23 * (4 * t1 * a[2] * t9 + 2 * a[1] * mu + 
                           2 * t15 * a[2] + t13 * t1 + 2 * t3 * t2 + 
                           2 * t6 * t2 - t3 - t6 - 2 * t9) / 2)
  t28 <- sqrt(2)
  t30 <- sqrt(-t23)
  t33 <- sqrt(pi)
  t36 <- log(t33 / t30 * sigma * t28 * t27)
  return(t36)
}

J <- 2
R <- 6
set.seed(1)



  
  # EIS fitting process
aa <- matrix(0.0,nrow = T,ncol=2)

zz <- matrix(rnorm(R*T),R,T)
xx <- matrix(0.0,R,T)



# initial conditions
for(t in 1:T){
  aa[t,1] <- 0.5*log(y[t]^2)
  aa[t,2] <- -0.25
}
# EIS iterations
for(iter in 1:J){
  print(paste0("EIS iter #",iter))
  # sampling pass
  # t=1
  xx[,1] <- eisMean(mu1,sigma1,aa[1,]) + eisStd(sigma1,aa[1,])*zz[,1] 
  # remaining ts
  for(t in 2:T){
    xx[,t] <- eisMean(gam+delta*xx[,t-1],v,aa[t,]) + eisStd(v,aa[t,])*zz[,t]
  }
  
  # regression pass
  lchi <- rep(0.0,R)
  for(t in T:1){
    rhs <- matrix(1.0,R,3)
    rhs[,2] <- xx[,t]
    rhs[,3] <- xx[,t]^2
    
    lhs <- lchi + gg(y[t],xx[,t])
    
    fit <- lm.fit(x=rhs,y=lhs)
    aa[t,] <- fit$coefficients[2:3]
    
    if(t>1){
      lchi <- eislchi(gam+delta*xx[,t-1],v,aa[t,])
    }
  }
}


# now compute log-weights
lwts.EIS <- rep(0.0,nrep)
for(rep in 1:nrep){
  
x <- eisMean(mu1,sigma1,aa[1,]) + eisStd(sigma1,aa[1,])*rnorm(1) 
ll <- eislchi(mu1,sigma1,aa[1,]) + gg(y[1],x) - aa[1,1]*x - aa[1,2]*x^2
for(t in 2:T){
 
  tmp <- eislchi(gam+delta*x,v,aa[t,])
  x <- eisMean(gam+delta*x,v,aa[t,]) + eisStd(v,aa[t,])*rnorm(1)
  ll <- ll + tmp + gg(y[t],x) - aa[t,1]*x - aa[t,2]*x^2
}
lwts.EIS[rep] <- ll
}


# prior
lwts.prior <- rep(0.0,nrep)
for(i in 1:nrep){
  x <- mu1 + sigma1*rnorm(1)
  ll <- gg(y[1],x)
  for(t in 2:T){
    x <- gam + delta*x + v*rnorm(1)
    ll <- ll + gg(y[t],x)
  }
  lwts.prior[i] <- ll
}


#Laplace K=0
lwts.lap0 <- rep(0.0,nrep)

vSq <- v^2
deltaSq <- delta^2
GxGy <- rep(0.0,T)
L <- CIP_TriDiagChol_const1n(T, 0.5+(1.0/vSq), 0.5+(1.0+deltaSq)/vSq, -delta/vSq);
#// G_x*h_x + G_y*h_y
GxGy[1] <- log(abs(y[1]))+(gam/vSq);
GxGy[2:(T-1)] <- log(abs(y[2:(T-1)])) + gam*(1.0-delta)/vSq;
GxGy[T] <- log(abs(y[T]))+(gam/vSq);
h <- CIP_TriDiagChol_LLT_solve(L,GxGy);


for(i in 1:nrep){
  # simulate trajectory
  u <- rnorm(T)
  x <- h+CIP_TriDiagChol_LT_solve(L,u);
  # importance weight
  lp <- dnorm(x[1],mean=mu1,sd=sigma1,log=TRUE) + gg(y[1],x[1])
  for(t in 2:T){
    lp <- lp + dnorm(x[t],mean=gam+delta*x[t-1],sd=v,log=TRUE) + gg(y[t],x[t])
  }
  lp <- lp + 0.5*sum(u^2) + 0.5*T*log(2*pi) - L[2*T]
  lwts.lap0[i] <- lp
}

# Laplace K
K <- 2
x <- h

for(k in 1:K){
  
  
  Gdiag <- 0.5*y^2*exp(-x) + c(1.0/(v^2),rep((1.0+delta^2)/(v^2),T-2),1.0/(v^2))
  Godiag <- rep(-delta/(v^2),T-1)
  L <- CIP_TriDiagChol(Gdiag,Godiag) 
  ngrad <- 0.5 - 0.5*y^2*exp(-x)
  ngrad[1] <- ngrad[1] + (x[1]-mu1)/sigma1^2
  ngrad[2:T] <- ngrad[2:T] + (x[2:T] - gam - delta*x[1:(T-1)])/(v^2)
  ngrad[1:(T-1)] <- ngrad[1:(T-1)] - (delta/(v^2))*(x[2:T]-gam-delta*x[1:(T-1)])
  
  print(max(abs(ngrad)))
  
  x <- x - CIP_TriDiagChol_LLT_solve(L,ngrad)
  
}

lwts.lapK <- rep(0.0,nrep)

h <- x
for(i in 1:nrep){
  # simulate trajectory
  u <- rnorm(T)
  x <- h+CIP_TriDiagChol_LT_solve(L,u);
  # importance weight
  lp <- dnorm(x[1],mean=mu1,sd=sigma1,log=TRUE) + gg(y[1],x[1])
  for(t in 2:T){
    lp <- lp + dnorm(x[t],mean=gam+delta*x[t-1],sd=v,log=TRUE) + gg(y[t],x[t])
  }
  lp <- lp + 0.5*sum(u^2) + 0.5*T*log(2*pi) - L[2*T]
  lwts.lapK[i] <- lp
}

print("EIS")
print(sd(lwts.EIS))
print("prior")
print(sd(lwts.prior))
print("Laplace0")
print(sd(lwts.lap0))
print("LaplaceK")
print(sd(lwts.lapK))
