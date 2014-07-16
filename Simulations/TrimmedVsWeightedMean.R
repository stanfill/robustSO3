source('~/robustSO3/Source_Code/robustFunctions.R')
Rcpp::sourceCpp('Source_Code/robustCpp.cpp')
library(rotations)
library(reshape2)
library(plyr)

B<-100
kappa <- 100
distj<-"Cayley"
res<-data.frame(p=rep(c(0,.1,.2),each=100),Trim=0,Weight=0)
S21 <- genR(pi/8)
a <- .1

for(i in 1:nrow(res)){
  Qs<-ruarsCont(n=25,rangle=rfisher,kappa1=kappa,p=res$p[i],Scont=S21,
                S=id.SO3,kappa2=kappa,space='Q4')  
  
  trimS <- trimMean(Qs,a,method='anneal')  
  ws <- HnFun(Qs)
  weightS <- weighted.mean(Qs, w=1/sqrt(ws))
  
  res$Trim[i]<-rot.dist(trimS$Shat,method='intrinsic')
  res$Weight[i]<-rot.dist(weightS,method='intrinsic')
}


ddply(res,.(p),summarize,Trim=mean(Trim),Weight=mean(Weight))

qplot(Trim,Weight,data=res,facets=.~p)+geom_abline(intercept=0,slope=1)
