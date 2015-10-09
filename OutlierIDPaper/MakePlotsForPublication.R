library(ggplot2)
library(xtable)
library(reshape2)
#########von Mises results
#Incorporate parametric results when correct distributional assumption is made
load("~/robustSO3/OutlierIDPaper/Results/vonMisesResultsWithNonpara_22_9_15.RData")

allResvM <- rbind(compSum,sumRes)
qplot(Angle,Power,data=allResvM,colour=TMethod,group=TMethod,geom='line',size=I(1))+
  geom_hline(yintercept=c(0,0.05),colour="gray50")+theme_bw()+ylab(expression(Pr(Reject~H[0])))+
  scale_x_continuous(breaks=rstar,labels=expression(0,pi/8,pi/4,pi/2,3~pi/4))+
  scale_colour_discrete(name="")+facet_grid(nF~KappaF,labeller=label_parsed)+theme(legend.position='top')

#New plot without kappa=1 or r^*=pi/8
notallResvM <- subset(allResvM,Kappa>1&Angle!=pi/8)
qplot(Angle,Power,data=notallResvM,colour=TMethod,group=TMethod,geom='line',size=I(1))+
  geom_hline(yintercept=c(0,0.05),colour="gray50")+theme_bw()+ylab(expression(Pr(Reject~H[0])))+
  scale_x_continuous(breaks=rstar[-2],labels=expression(0,pi/4,pi/2,3~pi/4))+
  scale_colour_discrete(name="")+facet_grid(nF~KappaF,labeller=label_parsed)+theme(legend.position='top')
#ggsave("C:/Users/Sta36z/Dropbox/SO3_Papers/OutlierID/Figures/vonMisesWithNonparametricRed.pdf",width=9,height=4.5)

#############
#Incorporate parametric results when incorrect distributional assumption is made
load("~/robustSO3/OutlierIDPaper/Results/vMisesResultsIncorrectAss_23_9_15.RData")

allResvM <- rbind(compSum,sumRes)
qplot(Angle,Power,data=allResvM,colour=TMethod,group=TMethod,geom='line',size=I(1))+
  geom_hline(yintercept=c(0,0.05),colour="gray50")+theme_bw()+ylab(expression(Pr(Reject~H[0])))+
  scale_x_continuous(breaks=rstar,labels=expression(0,pi/8,pi/4,pi/2,3~pi/4))+
  scale_colour_discrete(name="")+facet_grid(nF~KappaF,labeller=label_parsed)+theme(legend.position='top')

#New plot without kappa=1 or r^*=pi/8
notallResvM <- subset(allResvM,Kappa>1&Angle!=pi/8)
qplot(Angle,Power,data=notallResvM,colour=TMethod,group=TMethod,geom='line',size=I(1))+
  geom_hline(yintercept=c(0,0.05),colour="gray50")+theme_bw()+ylab(expression(Pr(Reject~H[0])))+
  scale_x_continuous(breaks=rstar[-2],labels=expression(0,pi/4,pi/2,3~pi/4))+
  scale_colour_discrete(name="")+facet_grid(nF~KappaF,labeller=label_parsed)+theme(legend.position='top')

#ggsave("C:/Users/Sta36z/Dropbox/SO3_Papers/OutlierID/Figures/vMisesWrongAssWithNonparametricRed.pdf",width=9,height=4.5)

########
rm(list=ls())
load("~/robustSO3/OutlierIDPaper/Results/vonMisesResultsWithNonpara_22_9_15.RData")
allResvM <- rbind(compSum,sumRes)
notallResvM <- subset(allResvM,Kappa>1&Angle!=pi/8)

#Make tables
notallResvM$FAngle <- factor(notallResvM$Angle,labels=c("0","pi/4","pi/2","3pi/4"))
vmTableDF <- dcast(notallResvM,n+Kappa+Type~FAngle+Method,value.var="Power")
vmTableDF$n <- as.factor(vmTableDF$n); vmTableDF$Kappa <- as.factor(vmTableDF$Kappa)
print(xtable(vmTableDF,caption="Simulation results for the von Mises distribution.",
             align=c(rep("4",3),rep("r",13))), include.rownames=FALSE, 
      floating.environment = "sidewaystable")
