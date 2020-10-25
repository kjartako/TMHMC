library(mvtnorm)
library(Rcpp)
library(RcppEigen)
library(doParallel)
#Sys.setenv(PKG_CXXFLAGS = "-std=c++11")

y=read.table('../../CEV_data.txt')[,1]
n_y=length(y)
n_z=7
n_u=1
Delta=1/252

M_inv = matrix(c(7.956034e-05,1.404987e-03,-4.008700e-05,-1.879545e-04,-4.040866e-05,1.404987e-03,2.964462e-02,-6.625706e-04,
                 -3.139953e-03,-6.161980e-04,-4.008700e-05,-6.625706e-04,3.457882e-03,1.723964e-02,2.002834e-03,-1.879545e-04,
                 -3.139953e-03,1.723964e-02,8.777712e-02,8.176455e-03,-4.040866e-05,-6.161980e-04,2.002834e-03,8.176455e-03,
                 6.937637e-03),nrow=5)
M = solve(M_inv)
Par0=c(0.009904744,0.169023234,1.185682489,-1.806509839,-15.061183329)
npar=length(Par0)

n_reg=1
chains = 8
n_burnin=1000
n_sample=10000
N = n_burnin+n_sample
Ls = rep(3,N)
epsilons = rep(0.57,N)

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
  p_Par=rmvnorm(N,c(0,0,0,0,0),M)
  p_u = matrix(rnorm(N*n_u*n_y,0,1),nrow=N)
  u0 = rnorm(n_y*n_u,0,1)
  unif = runif(N)
  chainres=HMC_EIS(N,inits[i,],y,n_y,u0,z,n_u,n_z,p_Par,p_u,epsilons,Ls,unif,n_reg,M_inv,Delta)
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
  res_EIS[[i]]$Par[,4]=exp(0.5*res_EIS[[i]]$Par[,4])
  res_EIS[[i]]$Par[,5]=exp(0.5*res_EIS[[i]]$Par[,5])
  
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
  p_Par=rmvnorm(N,c(0,0,0,0,0),M)
  p_u = matrix(rnorm(N*n_y,0,1),nrow=N)
  u0 = rnorm(n_y,0,1)
  unif = runif(N)
  chainres=HMC_Laplace(N,inits[1,],y,y^2,n_y,u0,p_Par,p_u,epsilons,Ls,unif,M_inv,Delta,Newton_it=2)
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
  res_Laplace[[i]]$Par[,4]=exp(0.5*res_Laplace[[i]]$Par[,4])
  res_Laplace[[i]]$Par[,5]=exp(0.5*res_Laplace[[i]]$Par[,5])
  
  sample_i = res_Laplace[[i]]$Par[(n_burnin+1):N,]
  Laplace_mean[i,]=colMeans(sample_i)
  Laplace_sd[i,]=apply(sample_i,2,sd)
}

colMeans(Laplace_mean)
colMeans(Laplace_sd)

##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

