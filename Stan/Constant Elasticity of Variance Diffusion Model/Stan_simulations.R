library(rstan)
options(mc.cores = parallel::detectCores())
library(mvtnorm)
library(coda)

y=read.table('../../CEV_data.txt')[,1]
n_y=length(y)
nchains=8
###################################################################################
sm = stan_model('Stan_Laplace.stan',auto_write=TRUE)
###################################################################################

Delta=1/252

set.seed(1234)
M_inv = matrix(c(7.956034e-05,1.404987e-03,-4.008700e-05,-1.879545e-04,-4.040866e-05,1.404987e-03,2.964462e-02,-6.625706e-04,
                 -3.139953e-03,-6.161980e-04,-4.008700e-05,-6.625706e-04,3.457882e-03,1.723964e-02,2.002834e-03,-1.879545e-04,
                 -3.139953e-03,1.723964e-02,8.777712e-02,8.176455e-03,-4.040866e-05,-6.161980e-04,2.002834e-03,8.176455e-03,
                 6.937637e-03),nrow=5)
Par0=c(0.009904744,0.169023234,1.185682489,-1.806509839,-15.061183329)
inits=rmvnorm(nchains,Par0,M_inv)

myinits = list()
for (i in 1:nchains)
{
  myinits=c(myinits,list(list(alpha=inits[i,1],beta=inits[i,2],gamma=inits[i,3],logSxSq=inits[i,4],logSySq=inits[i,5])))
}

###################################################################################
###################################################################################
###################################################################################
res = sampling(sm,data=list(N=n_y,y=y,Delta=Delta,Nreg=1),iter=11000, warmup = 1000, seed=123, chains=nchains,
                  init=myinits,pars=c("alpha","beta","gamma","sx","sy"))

coda=As.mcmc.list(res)
summary(coda)$statistics
