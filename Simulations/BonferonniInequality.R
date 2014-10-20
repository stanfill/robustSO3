source('Source_Code/robustFunctions.R')
Rcpp::sourceCpp('Source_Code/BootstrapCpp.cpp')
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
  geom_hline(yintercept=c(0,0.05),colour="gray50")+theme_bw()+ylab(expression(Pr(Reject~H[0])))


######
#Compare Bonferonni cut off to using parametric bootstrap to estimate critical value
B <- 100
rstar <- c(0,pi/8,pi/4,pi/2,3*pi/4)
pvalExtBon <- data.frame(Zero=rep(0,B),"pi/8"=0,"pi/4"=0,"pi/2"=0,"3pi/4"=0) 
HnIntBon <- HnExtBon <- pvalIntBon <- pvalExtBon
HnIntBoot <- HnExtBoot <- pvalIntBoot <- pvalExtBoot <- pvalExtBon
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
    
    #Intrinsic Bonferonni
    HnIntBon[i,j] <- max(discord(RsOut,type='int'))
    pvalIntBon[i,j] <- n*pf(HnIntBon[i,j],3,3*(n-2),lower.tail=FALSE)
    
    #Extrinsic Bonferonni
    HnExtBon[i,j] <- max(discord(RsOut,type='ext'))
    pvalExtBon[i,j] <- n*pf(HnExtBon[i,j],3,3*(n-2),lower.tail=FALSE)
    
    #Intrinsic Bootstrap
    HnIntBoot[i,j] <- max(discord(RsOut,type='int'))
    HnBootInt <- HnBootCpp(RsOut,m,1,rangle)
    pvalIntBoot[i,j] <- length(which(HnBootInt>=HnIntBoot[i,j]))/m
    
    #Extrinsic Bootstrap
    HnExtBoot[i,j] <- max(discord(RsOut,type='ext'))
    HnBootExt <- HnBootCpp(RsOut,m,2,rangle)
    pvalExtBoot[i,j] <- length(which(HnBootExt>=HnExtBoot[i,j]))/m
    
  }
  
}

ExtMBon <- melt(pvalExtBon,variable.name="Angle",value.name="Pval",measure.vars=1:ncol(pvalExtBon))
ExtMBon$Type <- "Extrinsic"
ExtMBon$Method <- "Bon"

IntMBon <- melt(pvalIntBon,variable.name="Angle",value.name="Pval",measure.vars=1:ncol(pvalIntBon))
IntMBon$Type <- "Intrinsic"
IntMBon$Method <- "Bon"

ExtMBoot <- melt(pvalExtBoot,variable.name="Angle",value.name="Pval",measure.vars=1:ncol(pvalExtBoot))
ExtMBoot$Type <- "Extrinsic"
ExtMBoot$Method <- "Boot"

IntMBoot <- melt(pvalIntBoot,variable.name="Angle",value.name="Pval",measure.vars=1:ncol(pvalIntBoot))
IntMBoot$Type <- "Intrinsic"
IntMBoot$Method <- "Boot"

compDF <- rbind(ExtMBon,IntMBon,ExtMBoot,IntMBoot)
compSum <- ddply(compDF,.(Type,Method,Angle),summarize,Power=length(which(Pval<=0.05))/length(Pval))
compSum$TMethod <- paste(compSum$Type,compSum$Method)

#compSum50 <- compSum
#compSum50$Kappa <- "50"
#compSum <- rbind(compSum5,compSum50)

qplot(Angle,Power,data=compSum,colour=TMethod,group=TMethod,geom='line',size=I(1.5))+
  geom_hline(yintercept=c(0,0.05),colour="gray50")+theme_bw()+ylab(expression(Pr(Reject~H[0])))+
  scale_x_discrete(labels=expression(0,pi/8,pi/4,pi/2,3~pi/4))+scale_colour_discrete(guide=guide_legend(title=""))


qplot(Angle,Power,data=compSum,colour=Type,group=Type,geom='line',size=I(2),facets=.~Method)+
  geom_hline(yintercept=c(0,0.05),colour="gray50")+theme_bw()+ylab(expression(Pr(Reject~H[0])))

