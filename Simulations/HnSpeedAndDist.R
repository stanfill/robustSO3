#Is it faster to compute Hn using the eigenvalue method, or using the 
#full/reduced model F test method?
#So far the eigenvalue method is faster, but that may just be because it is in C++
#Not sure it is worth coding the F test version in C++

source('~/robustSO3/Source_Code/robustFunctions.R')
Rcpp::sourceCpp('Source_Code/robustCpp.cpp')
library(rotations)
library(microbenchmark)

HnApprox<-function(Qs){
  n<-nrow(Qs)
  Qhat<-mean(Qs)
  Hnia<-rep(0,n)
  SSR<-sum(rot.dist(Qs,Qhat,method='extrinsic',p=2))
  for(i in 1:n){
    
    Qsi<-Qs[-i,]
    Qhati<-mean(Qsi)
    SSRi<-sum(rot.dist(Qsi,Qhati,method='extrinsic',p=2))
    Hnia[i] <- (n-2)*(SSR - SSRi)/(SSRi)
    
  }
  return(Hnia)
}

Qs<-ruars(20,rcayley,space='Q4')
microbenchmark(HnFun(Qs),HnApprox(Qs))

Rs<-ruars(20,rcayley)
microbenchmark(HnFun(Rs),HnApprox(Rs))

#HnFun is much faster

#################################
#Distribution of Hn with/without outliers

#Without outlier (should be F with 3,3*(n-2) dfs)
Qs<-ruars(50,rcayley,kappa=50)
Hn<-HnFun(Qs)
x<-seq(0,max(Hn),length=length(Hn))
plot(ecdf(Hn))
lines(x,pf(x,3,3*(length(Hn)-2)))

#With outliers
QsW<-ruarsCont(n=25,rangle=rcayley,kappa1=100,p=res$p[i],Scont=id.SO3,
              S=id.SO3,kappa2=1,space='Q4') 
HnW<-HnFun(QsW)
xW<-seq(0,max(HnW),length=length(HnW))
plot(ecdf(HnW))
lines(xW,pf(xW,3,3*(length(HnW)-2)))

#######
#Distribution of Hn en blocs of size t with/without outliers
t <- 1

#Without outlier (should be F with 3,3*(n-2) dfs)
Qs<-ruars(50,rcayley,kappa=100)
Hn<-HnBlocCpp(Qs,t)
x<-seq(0,max(Hn$Hn),length=length(Hn$Hn))
plot(ecdf(Hn$Hn))
lines(x,pf(x,3*t,3*(length(Hn$Hn)-1-t)),col=2)

#With outliers
QsW<-ruarsCont(n=25,rangle=rcayley,kappa1=100,p=res$p[i],Scont=id.SO3,
               S=id.SO3,kappa2=1,space='Q4') 
HnW<-HnBloc(QsW,t)
xW<-seq(0,max(HnW),length=length(HnW))
plot(ecdf(HnW))
lines(xW,pf(xW,3,3*(length(HnW)-1-t)))
