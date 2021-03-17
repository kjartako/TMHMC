rm(list=ls())
source("choleskys.R")
nrep <- 1000


# gam model

y <- as.vector(read.table("CEV_data.txt")$V1)
T <- length(y)
alpha <- 0.01  
beta <- 0.17
gam <- 1.18
sigx <- 0.41
sigy <- 0.0005
Delta <- 1/252


gg <- function(y,x){
  return(dnorm(y,mean=x,sd=sigy,log=TRUE))
}


mu1 <- y[1]
sigma1 <- 0.01




# EIS

eisMean <- function(mu,sigma,a){
  return((mu + a[1]*sigma^2)/(1-2.0*a[2]*sigma^2))
}

eisStd <- function(sigma,a){
  return(sigma/sqrt(1.0-2.0*a[2]*sigma^2))
}

eislchi <- function(mu,sigma,a){
  vv <- sigma^2
  m <- mu
  a1 <- a[1]
  a2 <- a[2]
  t1 <- log(vv)
  t5 <- a1 ^ 2
  t7 <- m ^ 2
  t18 <- vv ^ 2
  t22 <- sqrt(-2 * t18 * a2 + vv)
  t25 <- log(0.1e1 / vv * t22)
  t26 <- 0.1e1 / (4 * a2 * vv - 2) * (-2 * vv * a2 * t1 - 2 * a1 * m - 2 * t7 * a2 - vv * t5 + t1) - t25
  
  return(t26)
}

J <- 2
R <- 7
set.seed(1)



  
  # EIS fitting process
aa <- matrix(0.0,nrow = T,ncol=2)

zz <- matrix(rnorm(R*T),R,T)
xx <- matrix(0.0,R,T)



# initial conditions
for(t in 1:T){
  aa[t,1] <- y[t]/(sigy^2)
  aa[t,2] <- -0.5/(sigy^2)
}
# EIS iterations
for(iter in 1:J){
  print(paste0("EIS iter #",iter))
  # sampling pass
  # t=1
  xx[,1] <- eisMean(mu1,sigma1,aa[1,]) + eisStd(sigma1,aa[1,])*zz[,1] 
  # remaining ts
  for(t in 2:T){
    vv <- sigx*sqrt(Delta)*xx[,t-1]^gam
    xx[,t] <- eisMean(xx[,t-1] + Delta*(alpha-beta*xx[,t-1]),vv,aa[t,]) + eisStd(vv,aa[t,])*zz[,t]
  }
  
  plot(xx[1,],type="l")
  for(i in 2:R) lines(xx[i,],col=i)
  points(y)
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
      vv <- sigx*sqrt(Delta)*xx[,t-1]^gam
      lchi <- eislchi(xx[,t-1] + Delta*(alpha-beta*xx[,t-1]),vv,aa[t,])
      
    }
  }
}


# now compute log-weights
lwts.EIS <- rep(0.0,nrep)
for(rep in 1:nrep){
  
x <- eisMean(mu1,sigma1,aa[1,]) + eisStd(sigma1,aa[1,])*rnorm(1) 
ll <- eislchi(mu1,sigma1,aa[1,]) + gg(y[1],x) - aa[1,1]*x - aa[1,2]*x^2
for(t in 2:T){
  vv <- sigx*sqrt(Delta)*x^gam
  tmp <- eislchi(x + Delta*(alpha-beta*x),vv,aa[t,])
  x <- eisMean(x + Delta*(alpha-beta*x),vv,aa[t,]) + eisStd(vv,aa[t,])*rnorm(1)
  ll <- ll + tmp + gg(y[t],x) - aa[t,1]*x - aa[t,2]*x^2
}
lwts.EIS[rep] <- ll
}


# prior (not used in this case)



#Laplace 

sxsq <- sigx^2;
sysq <- sigy^2;
sx <- sigx;
sy <- sigy;
gamsq <- gam^2;
Delbeta <- (Delta*beta-1.0);
DelbetaSq <- (Delta*beta-1.0)^2;

#// eval point
h<-y;
hSq<-y^2;

#// Newton iterations
tmp <- rep(0.0,T-1)
tmp2 <- rep(0.0,T-1)
diag <- rep(0.0,T)
grad <- rep(0.0,T)
for(w in 1:2)
{
  tmp[1:(T-1)] <- h[2:T] + Delbeta*h[1:(T-1)] - Delta*alpha;
  tmpSq <- tmp^2
  
  for(i in 1:(T-1))
  {
    tmp2[i] <- h[i]^(2.0*gam+2.0)
  }
  tmp2 <- tmp2*sxsq*Delta;
  
  #// build Hessian
  #// diagonal
  diag[1:(T-1)] = ((DelbetaSq*hSq[1:(T-1)]-4.0*Delbeta*gam*h[1:(T-1)]*tmp[1:(T-1)]
                    + (2.0*gamsq+gam)*tmpSq[1:(T-1)])/ tmp2[1:(T-1)])-(gam / hSq[1:(T-1)])
  diag[T] <- 0.0;
  diag[1] <- diag[1] + 10000.0;
  diag[2:T] <- diag[2:T] + hSq[1:(T-1)] / tmp2[1:(T-1)];
  diag[1:T] <- diag[1:T] + 1.0/sysq;
  
  #// first off-diagonal
  odiag = (Delbeta*hSq[1:(T-1)]-2.0*gam*h[1:(T-1)] * tmp[1:(T-1)]) / tmp2[1:(T-1)]
  
  #// Cholesky
  Lnewton = CIP_TriDiagChol(diag,odiag);
  
  #// build gradient
  
  grad[1] <- -10000.0*h[1]+10000.0*y[1] + (tmpSq[1]*gam*h[1] - tmp[1]*Delbeta*hSq[1])/tmp2[1]  - (gam/h[1]) + ((y[1]-h[1])/sysq)
  
  grad[2:(T-1)] <- - ((tmp[1:(T-2)] * hSq[1:(T-2)]) / tmp2[1:(T-2)]) + (gam*tmpSq[2:(T-1)] * h[2:(T-1)] - Delbeta*tmp[2:(T-1)] * hSq[2:(T-1)]) / tmp2[2:(T-1)]  - (gam / h[2:(T-1)]) + ((y[2:(T-1)]-h[2:(T-1)])/sysq);
  
  grad[T] <- - ((tmp[T-1]*hSq[T-1]) / tmp2[T-1]) + ((y[T]-h[T])/sysq);
  
  #// update eval point
  #print(max(abs(grad)))
  h = h + CIP_TriDiagChol_LLT_solve(Lnewton,grad);	
  hSq = (h^2);
  lwts.lapK <- rep(0.0,nrep)
  for(i in 1:nrep){
    # simulate trajectory
    u <- rnorm(T)
    x <- h+CIP_TriDiagChol_LT_solve(Lnewton,u);
    # importance weight
    lp <- dnorm(x[1],mean=mu1,sd=sigma1,log=TRUE) + gg(y[1],x[1])
    for(t in 2:T){
      vv <- sigx*sqrt(Delta)*x[t-1]^gam
      mn <- x[t-1] + Delta*(alpha-beta*x[t-1])
      lp <- lp + dnorm(x[t],mean=mn,sd=vv,log=TRUE) + gg(y[t],x[t])
    }
    lp <- lp + 0.5*sum(u^2) + 0.5*T*log(2*pi) - Lnewton[2*T]
    lwts.lapK[i] <- lp
  }
  print(paste0("Laplace, K = ",w))
  print(sd(lwts.lapK))
  t <- lwts.lapK # ESS calculations
  wts <- exp(t-max(t))
  wts <- wts/sum(wts)
  print(1/sum(wts^2))
  
}


print("EIS")
print(sd(lwts.EIS))
t <- lwts.EIS # ESS calculations
wts <- exp(t-max(t))
wts <- wts/sum(wts)
print(1/sum(wts^2))
