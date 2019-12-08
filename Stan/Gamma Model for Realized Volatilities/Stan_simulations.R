library(rstan)
options(mc.cores = parallel::detectCores())
library(mvtnorm)
library(coda)

y=read.table('../../Gamma_data.txt')[,1]
n_y=length(y)
nchains=8
###################################################################################
sm = stan_model('Stan_Prior.stan',auto_write=TRUE)
#sm = stan_model('Stan_Laplace.stan',auto_write=TRUE)
###################################################################################

s0 = 0.01
df=10

set.seed(1234)
M_inv = matrix(c(0.00209958921,0.00007827431,0.00119064360,-0.00259986234,0.00007827431,0.06069924229,-0.00005555461,-0.00014122489,0.00119064360,-0.00005555461,0.01155006118,-0.00392100711,-0.00259986234,-0.00014122489,-0.00392100711,0.00928497194),nrow=4)
Par0=c(-2.075677,0.8926154,2.3724145,-2.9600990)
inits=rmvnorm(nchains,Par0,M_inv)

myinits = list()
for (i in 1:nchains)
{
  #delta_p = (delta+1)/2, delta = tanh(Par0[3])
  #v_p = (s0*df)/v^2, v = exp(0.5*Par0[4])
  myinits = c(myinits,list( list( tau=exp(inits[i,1]), beta=exp(inits[i,2]), delta_p=(1+tanh(inits[i,3]))/2, v_p = (s0*df)/(exp(inits[i,4])) )))
}

###################################################################################
###################################################################################
###################################################################################
res = sampling(sm,data=list(N=n_y,y=y),iter=2000, warmup = 1000, seed=123, chains=nchains,init=myinits,pars=c("tau","beta","delta","v"))
#res = sampling(sm,data=list(N=n_y,y=y),iter=2000, warmup = 1000, seed=123, chains=nchains,init=myinits,pars=c("tau","beta","delta","v"),control = list(max_treedepth = 6))

coda=As.mcmc.list(res)
summary(coda)$statistics
