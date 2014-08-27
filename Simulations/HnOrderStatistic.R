source('~/robustSO3/Source_Code/robustFunctions.R')
library(rotations)
library(plyr)
library(reshape2)

######
#Use bootstrap to determine critical value for H_(n)
B <- 250
rstar <- c(0,pi/8,pi/4,pi/2)
pvalExt <- data.frame(Zero=rep(0,B),"pi/8"=0,"pi/4"=0,"pi/2"=0) 
HnInt <- HnExt <- pvalInt <- pvalExt
kap <- 5
n <- 20
m <- 200
#Sc <- genR(pi/2)

for(j in 1:length(rstar)){
  
  Sc <- genR(rstar[j])
  
  for(i in 1:B){
    
    #No outiers
    #Rs <- ruars(n,rcayley,kappa=kap)
    #Hnn[i] <- max(discord(Rs,type='int'))
    #HnBootObs <- HnBoot(Rs,m,type='int',parametric=TRUE)
    #pval[i] <- length(which(HnBootObs>=Hnn[i]))/m
    
    #One outlier
    RsOut <- ruarsCont(n,rcayley,kappa1=kap,p=1/n,Scont=Sc)
    
    #Intrinsic
    HnInt[i,j] <- max(discord(RsOut,type='int'))
    HnBootInt <- HnBoot(RsOut,m,type='int',parametric=TRUE)
    pvalInt[i,j] <- length(which(HnBootInt>=HnInt[i,j]))/m
    
    #Extrinsic
    HnExt[i,j] <- max(discord(RsOut,type='ext'))
    HnBootExt <- HnBoot(RsOut,m,type='ext',parametric=TRUE)
    pvalExt[i,j] <- length(which(HnBootExt>=HnExt[i,j]))/m
  }

}

ExtM <- melt(pvalExt,variable.name="Angle",value.name="Pval")
ExtM$Type <- "Extrinsic"

IntM <- melt(pvalInt,variable.name="Angle",value.name="Pval")
IntM$Type <- "Intrinsic"

compDF <- rbind(ExtM,IntM)
compSum <- ddply(compDF,.(Type,Angle),summarize,Power=length(which(Pval<0.05))/length(Pval))

qplot(Angle,Power,data=compSum,colour=Type,group=Type,geom='line')

###Old Plots
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

