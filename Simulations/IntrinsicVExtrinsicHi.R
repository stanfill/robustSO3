#Compare power of Intrinsic and Extrinsic Hi statistic
source('~/robustSO3/Source_Code/robustFunctions.R')
library(rotations)
library(reshape2)
library(plyr)

kappa<-c(1,2,5,10,20)
tau <- 1
n<-50
B<-500
resDF<-data.frame(kappa=rep(kappa,each=B),IID=0,ICut=0,Hi=0,EID=0,ECut=0,He=0)
rownum<-1
sc<-genR(pi/2)

for(j in 1:length(kappa)){
  for(i in 1:B){
  
    Qs<-ruarsCont(n,rcayley,kappa[j],kappa2=tau,p=1/n,S=id.SO3,Scont=id.SO3,space='Q4')
    
    Hne<-discord(Qs,type='extrinsic')
    Hni<-discord(Qs,type='intrinsic')
    
    resDF$Hi[rownum]<-max(Hni)
    resDF$IID[rownum]<-(which.max(Hni)==n)
    
    resDF$He[rownum]<-max(Hne)
    resDF$EID[rownum]<-(which.max(Hne)==n)
    
    rownum<-rownum+1
  }
}

#cut<-qf(0.95,3,3*(n-2))
cut <- qorderF(0.95,n=n,k=n,df1=3,df2=3*(n-2),ncp=0)

resDF$ICut<-resDF$Hi>cut
resDF$ECut<-resDF$He>cut
#power <- 1-pf(cut*tau/kappa,3,3*(n-2))
power <- 1-porderF(cut*tau/kappa,n=n,k=n,df1=3,df2=3*(n-2))

resDFSum<-ddply(resDF,.(kappa),summarize,Intrinsic=sum(ICut)/length(ICut),Extrinsic=sum(ECut)/length(ECut),
                IID=sum(IID)/length(IID),EID=sum(EID)/length(EID))
resDFSum$Theory <- power
resDFSum
resDF2<-melt(resDFSum,id.vars="kappa",measure.vars=c("Intrinsic","Extrinsic","Theory"),
             variable.name='Test',value.name='Power')

qplot(kappa,Power,data=resDF2,colour=Test,geom='line',size=I(2),xlab=expression(kappa))+
  theme_bw()+scale_x_continuous(breaks=kappa)#+coord_equal(20)

#ggsave("C:/Users/Sta36z/Dropbox/SO3_Papers/OutlierIDAcc/Figures/PowerPic.pdf",width=5.5,height=4)


##################
#Why are some Intrinsic Hi values negative?  Because obs close to the center.

Rs<-ruars(200,rcayley,kappa=1)
discord(Rs,type='intrinsic')
