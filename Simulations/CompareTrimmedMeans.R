source('~/robustSO3/Source_Code/robustFunctions.R')
Rcpp::sourceCpp('Source_Code/robustCpp.cpp')
library(rotations)
library(reshape2)
library(plyr)

B<-100
k1 <- 100
k2 <- 1
distj<-"Cayley"
res<-data.frame(p=rep(c(0,.1,.2),each=100),Once=0,enBloc=0,Anneal=0,Oen=0,OA=0,eA=0)
#S21 <- genR(pi/4)
S21 <- id.SO3
a <- .05

for(i in 1:nrow(res)){
  Qs<-ruarsCont(n=25,rangle=rfisher,kappa1=k1,p=res$p[i],Scont=S21,
                S=id.SO3,kappa2=k2,space='Q4')  
  
  once <- trimMean(Qs,a,method='once')
  ann <- trimMean(Qs,a,method='anneal')  
  enB <- trimMean(Qs,a,method='en')
  
  res$Once[i]<-rot.dist(once$Shat,method='intrinsic')
  res$enBloc[i]<-rot.dist(enB$Shat,method='intrinsic')
  res$Anneal[i]<-rot.dist(ann$Shat,method='intrinsic')
  
  res$Oen[i]<-as.numeric(enB$Qs==once$Qs)
  res$OA[i]<-as.numeric(ann$Qs==once$Qs)
  res$eA[i]<-as.numeric(ann$Qs==enB$Qs)
  print(i)
}

mdat<-melt(res,id.vars="p",measure.vars=c("Once","enBloc","Anneal"))

qplot(value,data=mdat,facets=variable~p)

qplot(Once,enBloc,data=res,facets=.~p)+geom_abline(intercept=0,slope=1)
qplot(Anneal,enBloc,data=res,facets=.~p)+geom_abline(intercept=0,slope=1)
qplot(Once,Anneal,data=res,facets=.~p)+geom_abline(intercept=0,slope=1)

ddply(res,.(p),summarize,Oen=sum(Oen)/length(Oen),OA=sum(OA)/length(OA),eA=sum(eA)/length(eA),
      Annea=mean(Anneal),Once=mean(Once),enB=mean(enBloc))
