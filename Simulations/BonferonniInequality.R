source('Source_Code/robustFunctions.R')
library(rotations)
library(plyr)
library(reshape2)

######
#Use Bonferonni inequality to determine critical value for H_(n), mean slippage
B <- 500
rstar <- c(0,pi/8,pi/4,pi/2,3*pi/4)
pvalExt <- data.frame(Zero=rep(0,B),"pi/8"=0,"pi/4"=0,"pi/2"=0,"3pi/4"=0) 
HnInt <- HnExt <- pvalInt <- pvalExt
kap <- 50
n <- 50
m <- 500
distribution <- 1
rangle <- rcayley
#Sc <- genR(pi/2)

for(j in 1:length(rstar)){
  
  Sc <- genR(rstar[j])
  
  for(i in 1:B){
    
    #One outlier
    RsOut <- ruarsCont(n,rangle,kappa1=kap,p=1/n,Scont=Sc)
    
    #Intrinsic
    HnInt[i,j] <- max(discord(RsOut,type='int'))
    pvalInt[i,j] <- n*pf(HnInt[i,j],3,3*(n-2),lower.tail=FALSE)
    
    #Extrinsic
    HnExt[i,j] <- max(discord(RsOut,type='ext'))
    pvalExt[i,j] <- n*pf(HnExt[i,j],3,3*(n-2),lower.tail=FALSE)
  }
  
}

ExtM <- melt(pvalExt,variable.name="Angle",value.name="Pval",measure.vars=1:ncol(pvalExt))
ExtM$Type <- "Extrinsic"

IntM <- melt(pvalInt,variable.name="Angle",value.name="Pval",measure.vars=1:ncol(pvalExt))
IntM$Type <- "Intrinsic"

compDF <- rbind(ExtM,IntM)
compSum <- ddply(compDF,.(Type,Angle),summarize,Power=length(which(Pval<=0.05))/length(Pval))

qplot(Angle,Power,data=compSum,colour=Type,group=Type,geom='line',size=I(2))+
  geom_hline(yintercept=c(0,0.05),colour="gray50")+theme_bw()
