source('Source_Code/robustFunctions.R')
Rcpp::sourceCpp('Source_Code/BootstrapCpp.cpp')
library(rotations)
library(plyr)
library(reshape2)

#########################
######
#Compare Bonferonni cut off to using parametric bootstrap to estimate critical value
rstar <- c(0,pi/8,pi/4,pi/2,3*pi/4) #discordant value misorientation angle
kap <- c(1,5,25,50)
n <- 50
B <- 250 #number of simulations to run per rstar/kappa combination

Vals <- matrix(0,nrow=length(rstar)*length(kap),ncol=B)
pvalExtBon <- cbind(data.frame(Kappa=rep(kap,each=length(rstar)),Angle=rep(rstar,length(kap))),Vals)
HnIntBon <- HnExtBon <- pvalIntBon <- pvalExtBon
HnIntBoot <- HnExtBoot <- pvalIntBoot <- pvalExtBoot <- pvalExtBon


m <- 500 #number of samples to use in parametric bootstrap
rangle <- rvmises
rownum <- 0
#Sc <- genR(pi/2)

for(k in 1:length(kap)){

  for(j in 1:length(rstar)){
    rownum <- rownum + 1
    Sc <- genR(rstar[j])
    
    for(i in 1:B){

      #One outlier
      RsOut <- ruarsCont(n,rangle,kappa1=kap[k],p=1/n,Scont=Sc)
      
      #Intrinsic Bonferonni
      HnIntBon[rownum,(i+2)] <- max(discord(RsOut,type='int'))
      pvalIntBon[rownum,(i+2)] <- n*pf(HnIntBon[rownum,(i+2)],1,1*(n-2),lower.tail=FALSE)
      
      #Extrinsic Bonferonni
      HnExtBon[rownum,(i+2)] <- max(discord(RsOut,type='ext'))
      pvalExtBon[rownum,(i+2)] <- n*pf(HnExtBon[rownum,(i+2)],1,1*(n-2),lower.tail=FALSE)
      
      #Intrinsic Bootstrap
      HnIntBoot[rownum,(i+2)] <- max(discord(RsOut,type='int'))
      HnBootInt <- HnBootCpp(RsOut,m,1,rangle)
      pvalIntBoot[rownum,(i+2)] <- length(which(HnBootInt>=HnIntBoot[rownum,(i+2)]))/m
      
      #Extrinsic Bootstrap
      HnExtBoot[rownum,(i+2)] <- max(discord(RsOut,type='ext'))
      HnBootExt <- HnBootCpp(RsOut,m,2,rangle)
      pvalExtBoot[rownum,(i+2)] <- length(which(HnBootExt>=HnExtBoot[rownum,(i+2)]))/m
      
    }    
  }
}

ExtMBon <- melt(pvalExtBon,id=c("Kappa","Angle"),value.name="Pval")
ExtMBon$Type <- "Extrinsic"
ExtMBon$Method <- "Bonferonni"

IntMBon <- melt(pvalIntBon,id=c("Kappa","Angle"),value.name="Pval")
IntMBon$Type <- "Intrinsic"
IntMBon$Method <- "Bonferonni"

ExtMBoot <- melt(pvalExtBoot,id=c("Kappa","Angle"),value.name="Pval")
ExtMBoot$Type <- "Extrinsic"
ExtMBoot$Method <- "Bootstrap"

IntMBoot <- melt(pvalIntBoot,id=c("Kappa","Angle"),value.name="Pval")
IntMBoot$Type <- "Intrinsic"
IntMBoot$Method <- "Bootstrap"

compDF <- rbind(ExtMBon,IntMBon,ExtMBoot,IntMBoot)
compSum <- ddply(compDF,.(Type,Method,Kappa,Angle),summarize,Power=length(which(Pval<=0.05))/length(Pval))
compSum$TMethod <- paste(compSum$Type,compSum$Method)

#compSum50 <- compSum
#compSum50$Kappa <- "50"
#compSum <- rbind(compSum5,compSum50)
compSum$KappaF <- factor(compSum$Kappa,labels=c("kappa==1","kappa==5","kappa==25","kappa==50"))

qplot(Angle,Power,data=compSum,colour=TMethod,group=TMethod,geom='line',size=I(1.5))+
  geom_hline(yintercept=c(0,0.05),colour="gray50")+theme_bw()+ylab(expression(Pr(Reject~H[0])))+
  scale_x_continuous(breaks=rstar,labels=expression(0,pi/8,pi/4,pi/2,3~pi/4))+
  scale_colour_discrete(name="")+facet_grid(~KappaF,labeller=label_parsed)+theme(legend.position='top')


#qplot(Angle,Power,data=compSum,colour=Type,group=Type,geom='line',size=I(2),facets=.~Method)+
  #geom_hline(yintercept=c(0,0.05),colour="gray50")+theme_bw()+ylab(expression(Pr(Reject~H[0])))