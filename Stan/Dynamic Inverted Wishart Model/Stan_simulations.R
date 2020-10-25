library(rstan)
options(mc.cores = parallel::detectCores())
library(coda)

raw = t(read.table("../../Wishart_data.txt"))
################################################################################

q = 5
T = dim(raw)[1]
rr = array(dim=c(q,q,T))
k = 1

for(j in 1:q)
{
  for(i in j:q)
  {
    rr[i,j,] = raw[,k]
    rr[j,i,] = rr[i,j,]
    k = k+1
  }
}

yinv = matrix(0.0,nrow=T*q,ncol=q)
yy = matrix(0.0,nrow=T*q,ncol=q)
ldets = 0.0
for(t in 1:T)
{
  yy[((t-1)*q+1):(t*q),] = rr[,,t]
  yinv[((t-1)*q+1):(t*q),] = solve(rr[,,t])
  ldets = ldets + 2.0*sum(log(diag(chol(rr[,,t]))))
}

################################################################################
################################################################################
################################################################################
nchains=8

sm = stan_model('Stan_Prior.stan',auto_write=TRUE) 
#sm = stan_model('Stan_Laplace.stan',auto_write=TRUE) 

res = sampling(sm,data=list(T=T,yy=yy,yinv=yinv,ldets=ldets),iter=11000, warmup = 1000,
               seed=123, chains=nchains,pars=c("mu","sigmaSq","delta","hts","nu","firstx","firstz","lastz"))  

# res = sampling(sm,data=list(T=T,yy=yy,yinv=yinv,ldets=ldets,niter=0),iter=11000, warmup = 1000,
#                 seed=123, chains=nchains,pars=c("mu","sigmaSq","delta","hts","nu","firstx","firstz","lastz"))  
# res = sampling(sm,data=list(T=T,yy=yy,yinv=yinv,ldets=ldets,niter=1),iter=11000, warmup = 1000,
#                 seed=123, chains=nchains,pars=c("mu","sigmaSq","delta","hts","nu","firstx","firstz","lastz"))  
# res = sampling(sm,data=list(T=T,yy=yy,yinv=yinv,ldets=ldets,niter=2),iter=11000, warmup = 1000,
#                 seed=123, chains=nchains,pars=c("mu","sigmaSq","delta","hts","nu","firstx","firstz","lastz"))  


coda=As.mcmc.list(res)
summary(coda)$statistics
