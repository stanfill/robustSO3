source('~/robustSO3/Source_Code/robustFunctions.R')
library(rotations)

######
#Use bootstrap to determine cut-off for H_(n)
B <- 500
Hnn <- HnnOut <- pval <- pvalOut <- rep(0,B)
kap <- 50
n <- 20
m <- 100
Sc <- genR(pi/4)

for(i in 1:B){
  
  #No outiers
  Rs <- ruars(n,rcayley,kappa=kap)
  Hnn[i] <- max(discord(Rs,type='ext'))
  HnBootObs <- HnBoot(Rs,m,type='ext')
  pval[i] <- length(which(HnBootObs>Hnn[i]))/m
  
  #One outlier
  RsOut <- ruarsCont(n,rcayley,kappa1=kap,p=1/n,Scont=Sc)
  HnnOut[i] <- max(discord(RsOut,type='ext'))
  HnBootOut <- HnBoot(RsOut,m,type='ext')
  pvalOut[i] <- length(which(HnBootOut>HnnOut[i]))/m
}

hist(pval,breaks=20,prob=T,main=mean(pval))
hist(pvalOut,breaks=20,prob=T,main=mean(pvalOut))


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

