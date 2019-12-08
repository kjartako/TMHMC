library(rstan)
options(mc.cores = parallel::detectCores())
library(mvtnorm)
library(coda)

y=read.table('../../SV_data.txt')[,1]
n_y=length(y)
nchains=8
###################################################################################
sm = stan_model('Stan_Prior.stan',auto_write=TRUE)
#sm = stan_model('Stan_Laplace.stan',auto_write=TRUE)
#sm = stan_model('Stan_Fisher.stan',auto_write=TRUE)
###################################################################################

s0 = 0.01
df=10

set.seed(1234)
M_inv = matrix(c(0.000103540,0.001842507,-0.002242537,0.001842507,0.041489870,-0.050155357,-0.002242537,-0.050155357,0.136014191),nrow=3)
Par0=c(-0.02032639,2.23705620,-3.85967962)
inits=rmvnorm(nchains,Par0,M_inv)

myinits = list()
for (i in 1:nchains)
{
  #delta_p = (delta+1)/2, delta = tanh(Par0[2])
  #v_p = (s0*df)/v^2, v = exp(0.5*Par0[3])
  myinits = c(myinits,list( list( gamma=inits[i,1],delta_p=(1+tanh(inits[i,2]))/2, v_p = (s0*df)/(exp(0.5*inits[i,3])^2) )))
}

###################################################################################
###################################################################################
###################################################################################
res = sampling(sm,data=list(N=n_y,y=y),iter=2000, warmup = 1000, seed=123, chains=nchains,init=myinits,pars=c("gamma","delta","v"))  

coda=As.mcmc.list(res)
summary(coda)$statistics