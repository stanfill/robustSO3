#First to show that the SSE/8 of wrt Euclidean distance for a sample of n rotations
#scaled by inverse sqrt of variance matrix has a central chi square 
#3n distribution as kappa increases.

library(rotations)
B<-1000
n<-50
kappa<-50
SSE<-SSEHat<-rep(0,B)

rs<-rcayley(1000,kappa=kappa)
a<-mean(cos(rs/2)^2)
b<-mean(sin(rs/2)^2)/3
InvSig<-diag(c(1/sqrt(a),rep(1/sqrt(b),3)))
A<-diag(c(0,rep(1,3)))

for(i in 1:B){
  qs<-ruars(n,rcayley,kappa=kappa,space='Q4')
  
  shat<-mean(qs)
  Ahati<-diag(1,4)-t(shat)%*%shat
    
  rshat<-rot.dist(qs,shat,method='intrinsic')
  ahat<-mean(cos(rshat/2)^2)
  bhat<-mean(sin(rshat/2)^2)/3
  InvSigHat<-diag(c(1/sqrt(ahat),rep(1/sqrt(bhat),3)))  
  
  for(j in 1:n){
    qsj<-matrix(qs[j,],4,1)
    SSE[i]<-SSE[i]+t(qsj)%*%t(InvSig)%*%A%*%InvSig%*%qsj
    SSEHat[i]<-SSEHat[i]+t(qsj)%*%t(InvSig)%*%Ahati%*%InvSig%*%qsj
  }
}

x<-seq(0,max(SSE),length=B)
plot(ecdf(SSE))
lines(x,pchisq(x,3*(n-1)),col=2)

plot(ecdf(SSEHat))
lines(x,pchisq(x,3*(n-1)),col=2)
lines(x,pchisq(x,3*n),col=3)

####
#Numerator of Hn is the difference in the squared SSEs when the jth observations
# is deleted minus the SSE with the full sample

#Compare the distribution of the difference in SSEs between the full model and when the 1st
#observation is deleted against the chisq distribution with 3 dfs.
#Then compare full/reduced model F test to the F distribution
B<-1000
n<-50
kappa<-25
SSEFull<-SSERed<-SSEDiff<-rep(0,B)

rs<-rfisher(1000,kappa=kappa)
a<-mean(cos(rs/2)^2)
b<-mean(sin(rs/2)^2)/3
InvSig<-diag(c(1/sqrt(a),rep(1/sqrt(b),3)))
A<-diag(c(0,rep(1,3)))

for(i in 1:B){
  qs<-ruars(n,rfisher,kappa=kappa,space='Q4')
  
  shat<-mean(qs)
  Ahati<-diag(1,4)-t(shat)%*%shat
  
  shat1<-mean(qs[-1,])
  Ahatij<-diag(1,4)-t(shat1)%*%shat1
  
  qsj<-matrix(qs[1,],4,1)
  SSEFull[i]<-SSEFull[i]+t(qsj)%*%InvSig%*%Ahati%*%InvSig%*%qsj 
  
  for(j in 2:n){
    qsj<-matrix(qs[j,],4,1)
    SSEFull[i]<-SSEFull[i]+t(qsj)%*%InvSig%*%Ahati%*%InvSig%*%qsj
    SSERed[i]<-SSERed[i]+t(qsj)%*%InvSig%*%Ahatij%*%InvSig%*%qsj
  }
}

SSEDiff<-SSEFull-SSERed
Ftest <- SSEDiff/(SSERed/(n-2))

x<-seq(0,max(SSEDiff),length=B)
plot(ecdf(SSEDiff))
lines(x,pchisq(x,3*(1)),col=2)

x<-seq(0,max(Ftest),length=B)
plot(ecdf(Ftest))
lines(x,pf(x,3,3*(n-2)),col=2)

###########
#How does this calculation compare to Hn?
Rcpp::sourceCpp('Source_Code/robustCpp.cpp')

kap<-1000
qs<-ruars(n,rcayley,kappa=kap,space='Q4')

shat<-mean(qs)

rs<-rot.dist(qs,shat,method='intrinsic')
rs<-rcayley(1000,kappa=kap)
a<-mean(cos(rs/2)^2)
b<-mean(sin(rs/2)^2)/3
InvSig<-diag(c(1/sqrt(a),rep(1/sqrt(b),3)))


Ahati<-diag(1,4)-t(shat)%*%shat
SSEFull<-0
SSERed<-rep(0,n)

for(j in 1:n){
  qsj<-matrix(qs[j,],4,1)
  SSEFull<-SSEFull+t(qsj)%*%t(InvSig)%*%Ahati%*%InvSig%*%qsj
  
  shat1<-mean(qs[-j,])
  Ahatij<-diag(1,4)-t(shat1)%*%shat1
  InvSigj<-InvSig
  #rsj<-rot.dist(qs[-j,],shat1,method='intrinsic')
  #aj<-mean(cos(rsj/2)^2)
  #bj<-mean(sin(rsj/2)^2)/3
  #InvSigj<-diag(c(1/sqrt(aj),rep(1/sqrt(bj),3)))
  
  cnt<-(1:n)[-j]
  
  for(i in cnt){
    SSERed[i]<-SSERed[i]+t(qsj)%*%t(InvSigj)%*%Ahatij%*%InvSigj%*%qsj
  }
}

Hnish<-(n-2)*(SSEFull-SSERed)/SSERed
Hn<-HnCpp(qs)
plot(Hnish,Hn,asp=1);abline(0,1)

x<-seq(0,max(c(Hnish,Hn)),length=100)
plot(ecdf(Hnish))
lines(x,pf(x,3,3*(n-2)))
lines(ecdf(Hn),col=2)



######
#Is the projection matrix even idempotent?
qs<-ruars(20,rcayley,kappa=10,space='Q4')
shat<-matrix(mean(qs),4,1)
shatj<-matrix(mean(qs[-5,]),4,1)

#No, symmetric but not idempotent
Aj<-shatj%*%t(shatj)-shat%*%t(shat)
round(Aj%*%Aj-Aj,5)

#Yes, symmetric and idempotent
denAj<-diag(1,4)-shat%*%t(shat)
round(denAj%*%denAj-denAj,10)


###
#Compare scaled data and original data
Qs<-ruars(20,rcayley,space='Q4',kappa=50)
scaleQs<-Qs
ahat<-mean(cos(rshat/2)^2)
bhat<-mean(sin(rshat/2)^2)/3
InvSigHat<-diag(c(1/sqrt(ahat),rep(1/sqrt(bhat),3))) 

for(i in 1:nrow(Qs)){
  Qsi<-matrix(Qs[i,],4,1)
  scaleQs[i,]<-InvSigHat%*%Qsi
  #scaleQs[i,]<-scaleQs[i,]/sqrt(sum(scaleQs[i,]^2))
}
is.Q4(scaleQs)

plot(scaleQs,col=1)
plot(Qs,col=1)
