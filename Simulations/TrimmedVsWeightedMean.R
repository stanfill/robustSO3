source('~/robustSO3/Source_Code/robustFunctions.R')
#Rcpp::sourceCpp('Source_Code/robustCpp.cpp')
library(rotations)
library(reshape2)
library(plyr)

B<-100
k1 <- 100
k2 <- 100
n<-c(10,50,100)
p<-c(0,.1,.2)

res<-data.frame(p=rep(p,each=(B*length(n))),n=rep(n,(B*length(p))),ProjMean=0,ProjMedian=0,Trimmed=0,Weighted=0)
S21 <- genR(pi/2)
#S21 <- id.SO3
a <- .1


for(i in 1:nrow(res)){
  Qs<-ruarsCont(n=res$n[i],rangle=rcayley,kappa1=k1,p=res$p[i],Scont=S21,
                S=id.SO3,kappa2=k2,space='SO3')  
  
  trimS <- trimMean(Qs,a,method='anneal',type='intrinsic')  
  ws <- discord(Qs,type='i')
  weightS <- weighted.mean(Qs, w=1/sqrt(ws))
  
  res$Trimmed[i]<-rot.dist(trimS$Shat,method='intrinsic')
  res$Weighted[i]<-rot.dist(weightS,method='intrinsic')
  res$ProjMean[i]<-rot.dist(mean(Qs),method='intrinsic')
  res$ProjMedian[i]<-rot.dist(median(Qs),method='intrinsic')
}


(resSum<-ddply(res,.(p,n),summarize,ProjMean=mean(ProjMean),ProjMedian=mean(ProjMedian),
      Trimmed=mean(Trimmed),Weighted=mean(Weighted)))

resSumM<-melt(resSum,id.vars=c('p',"n"),variable.name='Estimator')

qplot(p,value,data=resSumM,colour=Estimator,geom='line',lwd=I(1.5),ylab='Bias',main='Cayley')+
  facet_wrap(~n,scales="free_y")+theme_bw()+scale_x_continuous(breaks=c(0,.1,.2))

#ggsave("C:/Users/sta36z/Dropbox/SO3_Papers/OutlierIDAcc/Figures/CayleySim.pdf",height=5,width=10)

