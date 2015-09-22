source('Source_Code/robustFunctions.R')
Rcpp::sourceCpp('Source_Code/BootstrapCpp.cpp')
library(rotations)
library(plyr)
library(reshape2)
library(xtable)

#########################
######
#How do the tests perform is the wrong distributional assumpition is used to 
#determine the p-value
rstar <- c(0,pi/8,pi/4,pi/2,3*pi/4) #discordant value misorientation angle
kap <- c(1,5,25,50)
n <- c(10,50)
B <- 1000 #number of simulations to run per rstar/kappa combination

Vals <- matrix(0,nrow=length(rstar)*length(kap)*length(n),ncol=B)
pvalExtBon <- cbind(data.frame(n=rep(n,each=length(kap)*length(rstar)),Kappa=rep(kap,each=length(rstar)),Angle=rstar),Vals)
HnIntBon <- HnExtBon <- pvalIntBon <- pvalExtBon
HnIntBoot <- HnExtBoot <- pvalIntBoot <- pvalExtBoot <- pvalExtBon


mEx <- 250 #number of samples to use in parametric bootstrap

#Generate the data from the cayley distribution
rangle <- rvmises

#Use the vmises distribution to generate bootstrap data and use it's limiting distribution
#in the Bonferonni approximation
rangleWrong <- rcayley
rownum <- 0

date()
for(l in 1:length(n)){
  for(k in 1:length(kap)){
    for(j in 1:length(rstar)){
      
      rownum <- rownum + 1
      Sc <- genR(rstar[j])
      
      for(i in 1:B){
        
        #One outlier
        RsOut <- ruarsCont(n[l],rangle,kappa1=kap[k],p=(1/n[l]),Scont=Sc)
        
        #Intrinsic Bonferonni
        HnIntBon[rownum,(i+3)] <- HnIntBoot[rownum,(i+3)] <- max(discord(RsOut,type='int'))
        pvalIntBon[rownum,(i+3)] <- n[l]*pf(HnIntBon[rownum,(i+3)],1,(n[l]-2),lower.tail=FALSE)
        
        #Extrinsic Bonferonni
        HnExtBon[rownum,(i+3)] <- HnExtBoot[rownum,(i+3)] <- max(discord(RsOut,type='ext'))
        pvalExtBon[rownum,(i+3)] <- n[l]*pf(HnExtBon[rownum,(i+3)],1,(n[l]-2),lower.tail=FALSE)
        
        #Intrinsic Bootstrap
        HnBootInt <- HnBootCpp(RsOut,mEx,1,rangleWrong)
        pvalIntBoot[rownum,(i+3)] <- length(which(HnBootInt>=HnIntBoot[rownum,(i+3)]))/mEx
        
        #Extrinsic Bootstrap
        HnBootExt <- HnBootCpp(RsOut,mEx,2,rangleWrong)
        pvalExtBoot[rownum,(i+3)] <- length(which(HnBootExt>=HnExtBoot[rownum,(i+3)]))/mEx
        
      }    
    }
  }
}

ExtMBon <- melt(pvalExtBon,id=c("n","Kappa","Angle"),value.name="Pval")
ExtMBon$Type <- "Extrinsic"
ExtMBon$Method <- "Bonferonni"

IntMBon <- melt(pvalIntBon,id=c("n","Kappa","Angle"),value.name="Pval")
IntMBon$Type <- "Intrinsic"
IntMBon$Method <- "Bonferonni"

ExtMBoot <- melt(pvalExtBoot,id=c("n","Kappa","Angle"),value.name="Pval")
ExtMBoot$Type <- "Extrinsic"
ExtMBoot$Method <- "Bootstrap"

IntMBoot <- melt(pvalIntBoot,id=c("n","Kappa","Angle"),value.name="Pval")
IntMBoot$Type <- "Intrinsic"
IntMBoot$Method <- "Bootstrap"

compDF <- rbind(ExtMBon,IntMBon,ExtMBoot,IntMBoot)
compSum <- ddply(compDF,.(Type,Method,n,Kappa,Angle),summarize,Power=length(which(Pval<=0.05))/length(Pval))
compSum$TMethod <- paste(compSum$Type,compSum$Method)

#compSum50 <- compSum
#compSum50$Kappa <- "50"
#compSum <- rbind(compSum5,compSum50)
compSum$KappaF <- factor(compSum$Kappa,labels=c("kappa==1","kappa==5","kappa==25","kappa==50"))
compSum$nF <- factor(compSum$n,labels=c("n==10","n==50"))

qplot(Angle,Power,data=compSum,colour=TMethod,group=TMethod,geom='line',size=I(1))+
  geom_hline(yintercept=c(0,0.05),colour="gray50")+theme_bw()+ylab(expression(Pr(Reject~H[0])))+
  scale_x_continuous(breaks=rstar,labels=expression(0,pi/8,pi/4,pi/2,3~pi/4))+
  scale_colour_discrete(name="")+facet_grid(nF~KappaF,labeller=label_parsed)+theme(legend.position='top')

#save.image("~/robustSO3/OutlierIDPaper/Results/vMisesResultsIncorrectAss_22_9_15.RData")
date()


#Make table of results
compTab <- compSum[,-c(7:9)]
t1 <- dcast(compTab,n+Kappa+Method+Type~Angle)
print(xtable(t1,digits=3),include.rownames=FALSE)
