library(rbi)
library(coda)


dta.raw <- as.vector(read.table("SV_data.txt"))
tmp <- cbind(as.vector(1:945),dta.raw)
colnames(tmp) <- c("time","value")
dta <- as.data.frame(tmp)

model <- bi_model("StochasticVolatility.bi")

bi <- libbi(model)

nsamples <- 20000
nburnin <- 10000
nrep <- 8
llist <- list()

timing <- rep(0.0,nrep)

for(rep in 1:nrep){
  bi_post <- sample(bi,target="posterior",
                  nparticles=1024,
                  nsamples=nsamples,
                  seed=rep,
                  obs=list(y=dta),
                  end_time=945,verbose=TRUE)
  post <- bi_read(bi_post)
  
  timing[rep] <- post$clock/(1e6)*(nsamples-nburnin)/nsamples # time in milliseconds
  
  gam <- post$gam_v$value
  del <- -1.0 + 2.0*post$phi_v_prime$value
  sig <- sqrt(post$sigmaSq_v$value)
  tmp <- cbind(gam[(nburnin+1):nsamples],del[(nburnin+1):nsamples],sig[(nburnin+1):nsamples])
  colnames(tmp) <- c("gam","delta","v")
  llist[[rep]] <- mcmc(data=tmp)
}


out <- as.mcmc.list(llist)

save.image("Computations")
