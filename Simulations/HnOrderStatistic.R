source('~/robustSO3/Source_Code/robustFunctions.R')
Rcpp::sourceCpp('Source_Code/BootstrapCpp.cpp')
library(rotations)
library(plyr)
library(reshape2)

######
#Use bootstrap to determine critical value for H_(n), mean slippage
B <- 200
rstar <- c(0,pi/8,pi/4,pi/2,3*pi/4)
pvalExt <- data.frame(Zero=rep(0,B),"pi/8"=0,"pi/4"=0,"pi/2"=0,"3pi/4"=0) 
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
    #HnBootInt <- HnBoot(RsOut,m,type='int',parametric=TRUE)
    HnBootInt <- HnBootCpp(RsOut,m,1)
    pvalInt[i,j] <- length(which(HnBootInt>=HnInt[i,j]))/m
    
    #Extrinsic
    HnExt[i,j] <- max(discord(RsOut,type='ext'))
    #HnBootExtR <- HnBoot(RsOut,m,type='ext',parametric=TRUE)
    HnBootExt <- HnBootCpp(RsOut,m,2)
    pvalExt[i,j] <- length(which(HnBootExt>=HnExt[i,j]))/m
  }

}

ExtM <- melt(pvalExt,variable.name="Angle",value.name="Pval",measure.vars=1:ncol(pvalExt))
ExtM$Type <- "Extrinsic"

IntM <- melt(pvalInt,variable.name="Angle",value.name="Pval",measure.vars=1:ncol(pvalExt))
IntM$Type <- "Intrinsic"

compDF <- rbind(ExtM,IntM)
compSum <- ddply(compDF,.(Type,Angle),summarize,Power=length(which(Pval<0.05))/length(Pval))

qplot(Angle,Power,data=compSum,colour=Type,group=Type,geom='line',size=I(2))+
  geom_hline(yintercept=c(0,0.05),colour="gray50")+theme_bw()

######
#Use bootstrap to determine critical value for H_(n), concentration slippage
B <- 250
tau <- c(1,5,25,50)
pvalExt <- data.frame("1"=rep(0,B),"5"=0,"25"=0,"50"=0) 
HnInt <- HnExt <- pvalInt <- pvalExt
kap <- 50
n <- 20
m <- 200


for(j in 1:length(tau)){
    
  for(i in 1:B){
    
    #No outiers
    #Rs <- ruars(n,rcayley,kappa=kap)
    #Hnn[i] <- max(discord(Rs,type='int'))
    #HnBootObs <- HnBoot(Rs,m,type='int',parametric=TRUE)
    #pval[i] <- length(which(HnBootObs>=Hnn[i]))/m
    
    #One outlier
    RsOut <- ruarsCont(n,rcayley,kappa1=kap,kappa2=tau[j],p=1/n,Scont=id.SO3)
    
    #Intrinsic
    HnInt[i,j] <- max(discord(RsOut,type='int'))
    #HnBootInt <- HnBoot(RsOut,m,type='int',parametric=TRUE)
    HnBootInt <- HnBootCpp(RsOut,m,1)
    pvalInt[i,j] <- length(which(HnBootInt>=HnInt[i,j]))/m
    
    #Extrinsic
    HnExt[i,j] <- max(discord(RsOut,type='ext'))
    #HnBootExt <- HnBoot(RsOut,m,type='ext',parametric=TRUE)
    HnBootExt <- HnBootCpp(RsOut,m,2)
    pvalExt[i,j] <- length(which(HnBootExt>=HnExt[i,j]))/m
  }
  
}

ExtM <- melt(pvalExt,variable.name="Tau",value.name="Pval",measure.vars=1:ncol(pvalExt))
ExtM$Type <- "Extrinsic"

IntM <- melt(pvalInt,variable.name="Tau",value.name="Pval",measure.vars=1:ncol(pvalExt))
IntM$Type <- "Intrinsic"

compDF <- rbind(ExtM,IntM)
compSum <- ddply(compDF,.(Type,Tau),summarize,Power=length(which(Pval<0.05))/length(Pval))

qplot(Tau,Power,data=compSum,colour=Type,group=Type,geom='line',size=I(2))+geom_hline(yintercept=0)+
  theme_bw()

######
#Compare C++ and R version of HnBootstrap
Rs<-ruars(100,rcayley,kappa=50)
RHn <- HnBoot(Rs,1000,type='int',parametric=TRUE)
CppHn <- HnBootCpp(Rs,1000)

par(mfrow=c(1,2))
hist(RHn)
hist(CppHn)
layout(1)
plot(sort(RHn),sort(CppHn));abline(0,1)

library(microbenchmark)
microbenchmark(HnBoot(Rs,10,type='int',parametric=TRUE),HnBootCpp(Rs,10))
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

