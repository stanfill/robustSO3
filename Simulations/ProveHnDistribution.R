#First to show that the SSE/8 of wrt Euclidean distance for a sample of n rotations
#scaled by inverse sqrt of variance matrix has a central chi square 
#3n distribution as kappa increases

library(rotations)
B<-1000
n<-50
kappa<-25
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
  
  #rshat<-rot.dist(qs,shat,method='intrinsic')
  
  #ahat<-mean(cos(rshat/2)^2)
  #bhat<-mean(sin(rshat/2)^2)/3
  #InvSigHat<-diag(c(1/sqrt(ahat),rep(1/sqrt(bhat),3)))  
  
  for(j in 1:n){
    qsj<-matrix(qs[j,],4,1)
    SSE[i]<-SSE[i]+t(qsj)%*%InvSig%*%A%*%InvSig%*%qsj
    SSEHat[i]<-SSEHat[i]+t(qsj)%*%InvSig%*%Ahati%*%InvSig%*%qsj
  }
}

x<-seq(0,max(SSE),length=B)
plot(ecdf(SSE))
lines(x,pchisq(x,3*n),col=2)

plot(ecdf(SSEHat))
lines(x,pchisq(x,3*(n-1)),col=2)
lines(x,pchisq(x,3*n),col=3)


#Numerator of Hn is the difference in the squared SSEs when the jth observations
# is deleted minus the SSE with the full sample

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
