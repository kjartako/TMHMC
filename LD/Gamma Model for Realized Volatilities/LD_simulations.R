library(mvtnorm)
library(Rcpp)
library(RcppEigen)
library(doParallel)
#Sys.setenv(PKG_CXXFLAGS = "-std=c++11")

y=read.table('../../Gamma_data.txt')[,1]
n_y=length(y)
n_z=5
n_u=1
chains = 8

M_inv = matrix(c(0.00209958921,0.00007827431,0.00119064360,-0.00259986234,0.00007827431,0.06069924229,-0.00005555461,-0.00014122489,0.00119064360,-0.00005555461,0.01155006118,-0.00392100711,-0.00259986234,-0.00014122489,-0.00392100711,0.00928497194),nrow=4)
M = solve(M_inv)
Par0=c(-2.075677,0.8926154,2.3724145,-2.9600990)
npar=length(Par0)

n_reg=2
n_burnin=500
n_sample=1000
N = n_burnin+n_sample
Ls = rep(3,N)
epsilons = rep(0.64,N)

##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
set.seed(1234)
inits=rmvnorm(chains,Par0,M_inv)

cl = makeCluster(detectCores())
registerDoParallel(cl, cores = detectCores())

res_EIS = foreach(i=1:chains,.packages=c("mvtnorm","Rcpp")) %dopar% {
  sourceCpp("HMC_EIS_Laplace.cpp") # Requires Rtools (https://cran.r-project.org/bin/windows/Rtools/)
  ptm = proc.time()
  z=matrix(rnorm(N*n_y*n_z,0,1),nrow=N)
  p_Par=rmvnorm(N,c(0,0,0,0),M)
  p_u = matrix(rnorm(N*n_u*n_y,0,1),nrow=N)
  u0 = rnorm(n_y*n_u,0,1)
  unif = runif(N)
  chainres=HMC_EIS(N,inits[i,],y,n_y,u0,z,n_u,n_z,p_Par,p_u,epsilons,Ls,unif,n_reg,M_inv)
  chainres$time = (proc.time() - ptm)[[3]]
  return(chainres)
}
stopCluster(cl)

#############

EIS_mean=matrix(NA,nrow=chains,ncol=npar)
EIS_sd=matrix(NA,nrow=chains,ncol=npar)

for (i in 1:chains)
{
  # Transform to original parameters
  res_EIS[[i]]$Par[,1] = exp(res_EIS[[i]]$Par[,1])
  res_EIS[[i]]$Par[,2] = exp(res_EIS[[i]]$Par[,2])
  res_EIS[[i]]$Par[,3] = tanh(res_EIS[[i]]$Par[,3])
  res_EIS[[i]]$Par[,4] = exp(0.5*res_EIS[[i]]$Par[,4])
  
  sample_i = res_EIS[[i]]$Par[(n_burnin+1):N,]
  EIS_mean[i,]=colMeans(sample_i)
  EIS_sd[i,]=apply(sample_i,2,sd)
}

colMeans(EIS_mean)
colMeans(EIS_sd)

##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
set.seed(1234)
inits=rmvnorm(chains,Par0,M_inv)

cl = makeCluster(detectCores())
registerDoParallel(cl, cores = detectCores())

res_Laplace = foreach(i=1:chains,.packages=c("mvtnorm","Rcpp")) %dopar% {
  sourceCpp("HMC_EIS_Laplace.cpp") # Requires Rtools (https://cran.r-project.org/bin/windows/Rtools/)
  ptm = proc.time()
  p_Par=rmvnorm(N,c(0,0,0,0),M)
  p_u = matrix(rnorm(N*n_y,0,1),nrow=N)
  u0 = rnorm(n_y,0,1)
  unif = runif(N)
  chainres=HMC_Laplace(N,inits[i,],y,n_y,u0,p_Par,p_u,epsilons,Ls,unif,M_inv,Newton_it=1)
  chainres$time = (proc.time() - ptm)[[3]]
  return(chainres)
}
stopCluster(cl)

#############

Laplace_mean=matrix(NA,nrow=chains,ncol=npar)
Laplace_sd=matrix(NA,nrow=chains,ncol=npar)

for (i in 1:chains)
{
  # Transform to original parameters
  res_Laplace[[i]]$Par[,1] = exp(res_Laplace[[i]]$Par[,1])
  res_Laplace[[i]]$Par[,2] = exp(res_Laplace[[i]]$Par[,2])
  res_Laplace[[i]]$Par[,3] = tanh(res_Laplace[[i]]$Par[,3])
  res_Laplace[[i]]$Par[,4] = exp(0.5*res_Laplace[[i]]$Par[,4])
  
  sample_i = res_Laplace[[i]]$Par[(n_burnin+1):N,]
  Laplace_mean[i,]=colMeans(sample_i)
  Laplace_sd[i,]=apply(sample_i,2,sd)
}

colMeans(Laplace_mean)
colMeans(Laplace_sd)

##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
