source('~/robustSO3/Source_Code/robustFunctions.R')
library(rotations)

######
#Use bootstrap to determine critical value for H_(n)
B <- 250
Hnn <- HnnOut <- pval <- pvalOut <- rep(0,B)
kap <- 50
n <- 20
m <- 100
Sc <- genR(pi/2)

for(i in 1:B){
  
  #No outiers
  Rs <- ruars(n,rcayley,kappa=kap)
  Hnn[i] <- max(discord(Rs,type='int'))
  HnBootObs <- HnBoot(Rs,m,type='int',parametric=TRUE)
  pval[i] <- length(which(HnBootObs>=Hnn[i]))/m
  
  #One outlier
  RsOut <- ruarsCont(n,rcayley,kappa1=kap,p=1/n,Scont=Sc)
  HnnOut[i] <- max(discord(RsOut,type='int'))
  HnBootOut <- HnBoot(RsOut,m,type='int',parametric=TRUE)
  pvalOut[i] <- length(which(HnBootOut>=HnnOut[i]))/m
}

par(mfrow=c(1,2))
hist(pval,breaks=20,prob=T,main=paste("FDR: ",length(which(pval<0.05))/B))
hist(pvalOut,breaks=20,prob=T,main=paste("Power: ",length(which(pvalOut<0.05))/B))


######
#Compare dist of H_(n) to porderF.  Can't use this, the Hi statistics need to be
#independent to use this formula, but they aren't
B <- 1000
Hnn <- rep(0,B)
kap <- 50
n <- 10

for(i in 1:B){
  Rs <- ruars(n,rcayley,kappa=kap)
  Hnn[i] <- max(discord(Rs,type='int'))
}

Hnn <- sort(Hnn)
plot(ecdf(Hnn))
lines(Hnn,porderF(Hnn,n=n,k=n,df1=3,df2=3*(n-2)),col=2)

