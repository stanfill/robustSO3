source('~/robustSO3/Source_Code/robustFunctions.R')
library(rotations)

######
#Compare dist of H_(n) to porderF
B <- 1000
Hnn <- rep(0,B)
kap <- 100
n <- 100

for(i in 1:B){
  Rs <- ruars(n,rcayley,kappa=kap)
  Hnn[i] <- max(discord(Rs,type='int'))
}

Hnn <- sort(Hnn)
plot(ecdf(Hnn))
lines(Hnn,porderF(Hnn,n=n,k=n,df1=3,df2=3*(n-2)),col=2)
