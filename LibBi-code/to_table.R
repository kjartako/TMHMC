library(coda)
load("Computations_new")


means <- matrix(0.0,3,8)
sds <- means
ESSs <- means
ESSspert <- means


for( i in 1:8){
  t <- summary(llist[[i]])
  means[,i] <- t$statistics[,"Mean"]
  sds[,i] <- t$statistics[,"SD"]
  ESSs[,i] <- effectiveSize(llist[[i]])
  ESSspert[,i] <- ESSs[,i]/timing[i]
}


print("time")
print(min(timing))
print(mean(timing))
print("mean, SD")
print(rowMeans(means))
print(rowMeans(sds))
print("min ESS")
print(round(min(ESSs[1,])))
print(round(min(ESSs[2,])))
print(round(min(ESSs[3,])))
print("min ESS/time")
print(min(ESSspert[1,]))
print(min(ESSspert[2,]))
print(min(ESSspert[3,]))

print("mean ESS")
print(round(rowMeans(ESSs)))
print("min ESS/time")
print(rowMeans(ESSspert))