library(mvtnorm)
library(Rcpp)
library(RcppEigen)
library(doParallel)
#Sys.setenv(PKG_CXXFLAGS = "-std=c++11")

s0 = 0.01
df=10
alpha=20
beta=1.5

y=read.table('../../SV_data.txt')[,1]
n_y=length(y)
n_z=6
n_u=1
chains = 8

M_inv = matrix(c(0.000103540,0.001842507,-0.002242537,0.001842507,0.041489870,-0.050155357,-0.002242537,-0.050155357,0.136014191),nrow=3)
M = solve(M_inv)
Par0=c(-0.02032639,2.23705620,-3.85967962)
npar=length(Par0)

n_yeg=2
n_burnin=1000
n_sample=10000
N = n_burnin+n_sample
Ls = rep(4,N)
epsilons = rep(0.4,N)

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
  p_Par=rmvnorm(N,c(0,0,0),M)
  p_u = matrix(rnorm(N*n_u*n_y,0,1),nrow=N)
  u0 = rnorm(n_y*n_u,0,1)
  unif = runif(N)
  chainres=HMC_EIS(N,inits[i,],y,n_y,u0,z,n_u,n_z,p_Par,p_u,epsilons,Ls,unif,n_yeg,M_inv,s0,df,alpha,beta)
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
  res_EIS[[i]]$Par[,2]=tanh(res_EIS[[i]]$Par[,2])
  res_EIS[[i]]$Par[,3]=exp(0.5*res_EIS[[i]]$Par[,3])

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
  p_Par=rmvnorm(N,c(0,0,0),M)
  p_u = matrix(rnorm(N*n_y,0,1),nrow=N)
  u0 = rnorm(n_y,0,1)
  unif = runif(N)
  chainres=HMC_Laplace(N,inits[i,],y,y^2,n_y,u0,p_Par,p_u,epsilons,Ls,unif,M_inv,s0,df,alpha,beta,Newton_it=2)
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
  res_Laplace[[i]]$Par[,2]=tanh(res_Laplace[[i]]$Par[,2])
  res_Laplace[[i]]$Par[,3]=exp(0.5*res_Laplace[[i]]$Par[,3])
  
  sample_i = res_Laplace[[i]]$Par[(n_burnin+1):N,]
  Laplace_mean[i,]=colMeans(sample_i)
  Laplace_sd[i,]=apply(sample_i,2,sd)
}

colMeans(Laplace_mean)
colMeans(Laplace_sd)

##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
