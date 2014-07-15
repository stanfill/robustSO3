library(rotations,lib='../BayesCR')
#library(rotations)
library(reshape2)
source("robustFunctions.R")
Rcpp::sourceCpp('robustCpp.cpp')

#Scont<-genR(pi)
#Rs<-ruarsCont(20,rcayley,10,20,.05,id.SO3,Scont)

#plot(Rs)+aes(size=Z,alpha=Z)+scale_size(limits=c(-1,1), range=c(0.5,2.5))+theme(legend.position='none')

eps<-c(0.05,.1,.15,.2,.45)
n<-c(20,50,100)
kappa1<-40
kappa2<-20
Scont<-c(pi/2,pi)

ResDfBias<-ResDfMSE<-expand.grid(Eps=eps,n=n,Sstar=Scont)
ResDfBias$Mean<-ResDfBias$Median<-ResDfBias$tMean<-ResDfBias$wMean<-0
ResDfMSE$Mean<-ResDfMSE$Median<-ResDfMSE$tMean<-ResDfMSE$wMean<-0

B<-1000

for(i in 1:nrow(ResDfMSE)){
  
  MeanBias<-MedianBias<-tMeanBias<-wMeanBias<-0
  Scont<-genR(ResDfMSE$Sstar[i])
  
  for(j in 1:B){
    
    Rs<-ruarsCont(ResDfMSE$n[i],rcayley,kappa1,kappa2,ResDfMSE$Eps[i],id.SO3,Scont)
    Shat<-mean(Rs,type='geometric')
    Stilde<-median(Rs,type='geometric')
    
    tMean<-trimMean(Rs,.1,DistToMedian,type='geometric')$Shat
    wMean<-winzMean(Rs,.1,DistToMedian,type='geometric')$Shat
    
    #tMean<-trimMean(Rs,.1,HnFun,type='geometric')$Shat
    #wMean<-winzMean(Rs,.1,HnFun,type='geometric')$Shat
      
    MeanBias[j]<-angle(Shat)
    MedianBias[j]<-angle(Stilde)
    tMeanBias[j]<-angle(tMean)
    wMeanBias[j]<-angle(wMean)
    
  }
  
  ResDfBias$Mean[i]<-mean(MeanBias)
  ResDfBias$Median[i]<-mean(MedianBias)
  ResDfBias$tMean[i]<-mean(tMeanBias)
  ResDfBias$wMean[i]<-mean(wMeanBias)
  
  ResDfMSE$Mean[i]<-mean(MeanBias^2)
  ResDfMSE$Median[i]<-mean(MedianBias^2)
  ResDfMSE$tMean[i]<-mean(tMeanBias^2)
  ResDfMSE$wMean[i]<-mean(wMeanBias^2)
}


mResDfBias<-melt(ResDfBias,id=c("Eps","n","Sstar"))
mResDfBias$Sstar<-factor(mResDfBias$Sstar,labels=c('pi/2','pi'))
colnames(mResDfBias)[4]<-"Estimator"
mResDfBias$Estimator<-factor(mResDfBias$Estimator,labels=c("Winsozrized\nMean","Trimmed\nMean","Median","Mean"))
qplot(Eps,value,data=mResDfBias,colour=Estimator,group=Estimator,geom='line',size=I(1.25),xlab=expression(epsilon),ylab='Bias')+
  facet_grid(Sstar~n,labeller = label_parsed,scales="free_y")+coord_fixed()

mResDfMSE<-melt(ResDfMSE,id=c("Eps","n","Sstar"))
mResDfMSE$Sstar<-factor(mResDfMSE$Sstar,labels=c('pi/2','pi'))
colnames(mResDfMSE)[4]<-"Estimator"
mResDfMSE$Estimator<-factor(mResDfMSE$Estimator,labels=c("Winsozrized\nMean","Trimmed\nMean","Median","Mean"))
qplot(Eps,value,data=mResDfMSE,colour=Estimator,group=Estimator,geom='line',size=I(1.25),xlab=expression(epsilon),ylab='MSE')+
  facet_grid(Sstar~n,labeller = label_parsed,scales="free_y")+coord_fixed()

mResDfBias$Measure<-'Bias'
mResDfMSE$Measure<-'MSE'
mResDF<-rbind(mResDfBias,mResDfMSE)

write.csv(mResDF,"Results/IntrinsicRobustSimulations_1000.csv")



