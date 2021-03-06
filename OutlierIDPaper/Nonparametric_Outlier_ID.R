library(diptest)
source('Source_Code/robustFunctions.R')
Rcpp::sourceCpp('Source_Code/NonParaBootstrapCpp.cpp')
library(rotations)
library(dplyr)
library(reshape2)

#########################
######
#Compare Bonferonni cut off to using parametric bootstrap to estimate critical value
rstar <- c(0,pi/8,pi/4,pi/2,3*pi/4) #discordant value misorientation angle
kap <- c(1,5,25,50)
n <- c(10,50)
B <- 1000 #number of simulations to run per rstar/kappa combination

Vals <- matrix(0,nrow=length(rstar)*length(kap)*length(n),ncol=B)
nonParaPvalInt <- cbind(data.frame(Type='Intrinsic',Method='Nonparametric',
                                   n=rep(n,each=length(kap)*length(rstar)),
                                   Kappa=rep(kap,each=length(rstar)),Angle=rstar),Vals)
nonParaPvalExt <- nonParaPvalInt
nonParaPvalExt$Type <- "Extrinsic"

rangle <- rvmises
rownum <- 0
#Sc <- genR(pi/2)

for(l in 1:length(n)){
  for(k in 1:length(kap)){
    for(j in 1:length(rstar)){
      
      rownum <- rownum + 1
      Sc <- genR(rstar[j])
      
      for(i in 1:B){
        
        #One outlier
        QsOut <- ruarsCont(n[l],rangle,kappa1=kap[k],p=1/n[l],Scont=Sc,space='Q4')
        
        #Intrinsic Nonparametric
        nonParaPvalInt[rownum,(i+5)] <- max(discord(QsOut,type='int'))
        Hsstar <- HnNonParaBootCpp(QsOut,m=250,type=1)[1,]
        
        while(any(is.infinite(Hsstar))||any(is.nan(Hsstar))){
          Hsstar <- HnNonParaBootCpp(QsOut,m=250,type=1)[1,]
          print("Inf or NAN detected in Int")
        }
        nonParaPvalInt[rownum,(i+5)] <- dip.test(Hsstar)$p.value
        
        #Extrinsic Nonparametric
        nonParaPvalExt[rownum,(i+5)] <- max(discord(QsOut,type='ext'))
        HsstarExt <- HnNonParaBootCpp(QsOut,m=250,type=2)[1,]
        
        while(any(is.infinite(HsstarExt))||any(is.nan(HsstarExt))){
          HsstarExt <- HnNonParaBootCpp(QsOut,m=250,type=2)[1,]
          print("Inf or NAN detected in Ext")
        }
        
        nonParaPvalExt[rownum,(i+5)] <- dip.test(HsstarExt)$p.value
        
      }    
    }
    cat("k is ",k," of ",length(kap),"\n")
  }
  cat("l is ",l," of ",length(n),"\n")
}


nonParaPval <- rbind(nonParaPvalExt,nonParaPvalInt)
res <- melt(nonParaPval,id=c("n","Kappa","Angle","Type","Method"),value.name="Pval")

sumRes <- res%>%group_by(n,Kappa,Angle,Type,Method)%>%summarize(Power=length(which(Pval<0.05))/length(Pval))

qplot(Angle,Power,data=sumRes,facets=n~Kappa,geom='line',colour=Type,group=Type)+theme_bw()


sumRes$TMethod <- paste(sumRes$Type,sumRes$Method)
sumRes$KappaF <- factor(sumRes$Kappa,labels=c("kappa==1","kappa==5","kappa==25","kappa==50"))
sumRes$nF <- factor(sumRes$n,labels=c("n==10","n==50"))

#############
#Incorporate parametric results when correct distributional assumption is made
load("~/robustSO3/OutlierIDPaper/Results/vonMisesResultsWithNonpara_22_9_15.RData")

allRes <- rbind(compSum,sumRes)
qplot(Angle,Power,data=allRes,colour=TMethod,group=TMethod,geom='line',size=I(1))+
  geom_hline(yintercept=c(0,0.05),colour="gray50")+theme_bw()+ylab(expression(Pr(Reject~H[0])))+
  scale_x_continuous(breaks=rstar,labels=expression(0,pi/8,pi/4,pi/2,3~pi/4))+
  scale_colour_discrete(name="")+facet_grid(nF~KappaF,labeller=label_parsed)+theme(legend.position='top')

#ggsave("C:/Users/Sta36z/Dropbox/SO3_Papers/OutlierID/Figures/CayleyWithNonparametric.pdf",width=9,height=4.5)

#New plot without kappa=1 or r^*=pi/8
notallRes <- subset(allRes,Kappa>1&Angle!=pi/8)
qplot(Angle,Power,data=notallRes,colour=TMethod,group=TMethod,geom='line',size=I(1))+
  geom_hline(yintercept=c(0,0.05),colour="gray50")+theme_bw()+ylab(expression(Pr(Reject~H[0])))+
  scale_x_continuous(breaks=rstar[-2],labels=expression(0,pi/4,pi/2,3~pi/4))+
  scale_colour_discrete(name="")+facet_grid(nF~KappaF,labeller=label_parsed)+theme(legend.position='top')

#ggsave("C:/Users/Sta36z/Dropbox/SO3_Papers/OutlierID/Figures/CayleyWithNonparametricRed.pdf",width=9,height=4.5)


#save.image("~/robustSO3/OutlierIDPaper/Results/vonMisesResultsWithNonpara_22_9_15.RData")
#############
#Incorporate parametric results when incorrect distributional assumption is made
load("~/robustSO3/OutlierIDPaper/Results/vMisesResultsIncorrectAss_23_9_15.RData")

allRes <- rbind(compSum,sumRes)
qplot(Angle,Power,data=allRes,colour=TMethod,group=TMethod,geom='line',size=I(1))+
  geom_hline(yintercept=c(0,0.05),colour="gray50")+theme_bw()+ylab(expression(Pr(Reject~H[0])))+
  scale_x_continuous(breaks=rstar,labels=expression(0,pi/8,pi/4,pi/2,3~pi/4))+
  scale_colour_discrete(name="")+facet_grid(nF~KappaF,labeller=label_parsed)+theme(legend.position='top')

#New plot without kappa=1 or r^*=pi/8
notallRes <- subset(allRes,Kappa>1&Angle!=pi/8)
qplot(Angle,Power,data=notallRes,colour=TMethod,group=TMethod,geom='line',size=I(1))+
  geom_hline(yintercept=c(0,0.05),colour="gray50")+theme_bw()+ylab(expression(Pr(Reject~H[0])))+
  scale_x_continuous(breaks=rstar[-2],labels=expression(0,pi/4,pi/2,3~pi/4))+
  scale_colour_discrete(name="")+facet_grid(nF~KappaF,labeller=label_parsed)+theme(legend.position='top')

#ggsave("C:/Users/Sta36z/Dropbox/SO3_Papers/OutlierID/Figures/vMisesWrongAssWithNonparametricRed.pdf",width=9,height=4.5)


#ggsave("C:/Users/Sta36z/Dropbox/SO3_Papers/OutlierID/Figures/CayleyWrongAssWithNonparametric.pdf",width=9,height=4.5)

#####################
#Distribution of bootstrap statistic?
rstar <- 0
kap <- 5
n <- 10
rangle <- rcayley
QsOut <- ruarsCont(n,rangle,kappa1=kap,p=1/n,Scont=genR(rstar),space='Q4')

Hsstar <- HnNonParaBootCpp(QsOut,m=1000,type=1)

hist(Hsstar,breaks=50)
mean(Hsstar);var(Hsstar[1,])

qqnorm(Hsstar);qqline(Hsstar)
