source('~/robustSO3/Source_Code/robustFunctions.R')
Rcpp::sourceCpp('Source_Code/robustCpp.cpp')
library(rotations)
library(reshape2)
library(plyr)

B<-100
k1 <- 100
k2 <- 1
distj<-"Cayley"
res<-data.frame(p=rep(c(0,.1,.2),each=B),Trim=0,Weight=0)
#S21 <- genR(pi/8)
S21 <- id.SO3
a <- .1

for(i in 1:nrow(res)){
  Qs<-ruarsCont(n=50,rangle=rcayley,kappa1=k1,p=res$p[i],Scont=S21,
                S=id.SO3,kappa2=k2,space='Q4')  
  
  trimS <- trimMean(Qs,a,method='anneal')  
  ws <- HnFun(Qs)
  weightS <- weighted.mean(Qs, w=1/sqrt(ws))
  
  res$Trim[i]<-rot.dist(trimS$Shat,method='intrinsic')
  res$Weight[i]<-rot.dist(weightS,method='intrinsic')
}


ddply(res,.(p),summarize,Trim=mean(Trim),Weight=mean(Weight))

qplot(Trim,Weight,data=res,facets=.~p)+geom_abline(intercept=0,slope=1)

n<-length(ws)
x<-seq(0,max(ws),length=n)
plot(ecdf(ws))
lines(x,pf(x,3,3*(n-2)))

toCut<-which(ws>qf(.95,3,3*(n-2)))
qs<-Qs[-c(45:50),]
ws2<-HnFun(qs)

n<-length(ws2)
x<-seq(0,max(ws2),length=n)
plot(ecdf(ws2))
lines(x,pf(x,3,3*(n-2)))


Qs2<-ruars(50,rcayley,kappa=10)
Hn<-HnFun(Qs2)
n<-length(Hn)

x<-seq(0,5,length=n)
plot(ecdf(Hn),xlim=c(0,5))
lines(x,pf(x,3,3*(n-2)))
