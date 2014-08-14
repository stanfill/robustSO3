#Show the Hi statistic has the F distribution
library(rotations)
n<-100
kap<-50
Qs<-ruars(n,rvmises,kappa=kap,space='Q4')
That<-t(Qs)%*%Qs
lamHat<-eigen(That)$values[1]

lamHati<-rep(0,n)
Hnum<-rep(0,n)
Hdenom<-rep(0,n)
HdenomApprox<-rep(0,n)
rsI2<-matrix(0,n,n-1)


for(i in 1:n){  
  Qi<-Qs[i,]
  Thati<-That-Qi%*%t(Qi)
  lamHati[i]<-eigen(Thati)$values[1]
  
  Qsi<-Qs[-i,]
  shati<-mean(Qsi)
  rsI2[i,]<-rot.dist(Qsi,shati,method='intrinsic',p=2)
  HdenomApprox[i]<-sum(rsI2[i,])/4
  
}

Hdenom<-(n-1-lamHati)
plot(Hdenom,HdenomApprox);abline(0,1)
Hdenom-HdenomApprox
Hnum<-1+lamHati-lamHat

#Does 2kr^2 have an approx chi^2_3 dist? Yes if it is Fisher, no if its Cayley.
stat<-sort(kap*(1-cos(sqrt(rsI2[1,]))^2))
#stat<-sort(kap*rsI2[1,])
plot(ecdf(stat))
lines(stat,pchisq(stat,1))

#Shape of the distribution as a function of kappa
kappa<-c(10,50,100)
rs<-seq(-pi/4,pi/4,length=100)
dfComp<-data.frame(r=rs,kap=rep(kappa,each=100))
dfComp$df<-dfisher(dfComp$r,kappa=dfComp$kap)
dfComp$kap<-as.factor(dfComp$kap)
qplot(r,df,data=dfComp,colour=kap,group=kap,geom='line')


#######################
#Proposition 3.4 From Leon et al 
#######################

library(rotations)
sfun<-function(Rs,S=id.SO3){
  S<-matrix(S,3,3)
  n<-nrow(Rs)
  s<-matrix(0,n,3)
  for(i in 1:n){
    Ri<-matrix(Rs[i,],3,3)
    inner<-S-2*solve(S+Ri)
    s[i,]<-c(-1,1,-1)*rev(inner[upper.tri(inner)])
  }
  return(s)
}


Rs<-ruars(50,rcayley,kappa=50)
s<-sfun(Rs-mean(Rs))
qqnorm(s[,3]);qqline(s[,3])
ss<-rowSums(s^2)
rs<-rot.dist(Rs,mean(Rs),method='intrinsic',p=2)

plot(ss,rs/4);abline(0,1)

########################
#LM approach again
########################
library(rotations)

Qs<-ruars(20,rcayley,space='Q4')
Qsi<-Qs[-1,]
Si<-matrix(mean(Qsi),nrow=1)
S<-matrix(mean(Qs),nrow=1)
A<-t(Si)%*%Si-t(S)%*%S

Qsq<-Qsi%*%A%*%t(Qsi)

#
rs<-rvmises(100,kappa=100)
1-mean(cos(rs))^2

#######################
#Old method
#######################
#First to show that the SSE/8 of wrt Euclidean distance for a sample of n rotations
#scaled by inverse sqrt of variance matrix has a central chi square 
#3n distribution as kappa increases.

library(rotations)
B<-1000
n<-50
kappa<-100
SSE<-SSEHat<-rep(0,B)

rs<-rcayley(1000,kappa=kappa)
a<-mean(cos(rs/2)^2)-mean(cos(rs/2))^2
#a<-var(cos(rs/2))
b<-mean(sin(rs/2)^2)/3

InvSig<-diag(c(1/sqrt(a),rep(1/sqrt(b),3)))
A<-diag(c(0,rep(1,3)))

for(i in 1:B){
  qsOrig<-ruars(n,rcayley,kappa=kappa,space='Q4')
  
  shat<-mean(qsOrig)
  qs<-qsOrig#-shat
  Ahati<-diag(1,4)-t(shat)%*%shat
    
  #rshat<-rot.dist(qs,id.Q4,method='intrinsic')
  #ahat<-mean(cos(rshat/2)^2)-mean(cos(rshat/2))^2
  #ahat<-var(cos(rshat/2))
  #bhat<-mean(sin(rshat/2)^2)/3
  #bhat<-var(sin(rshat/2))/3
  #InvSigHat<-diag(c(1/sqrt(ahat),rep(1/sqrt(bhat),3)))  
  
  for(j in 1:n){
    qsj<-matrix(qs[j,],4,1)
    SSE[i]<-SSE[i]+t(qsj)%*%t(InvSig)%*%A%*%InvSig%*%qsj
    conti<-t(qsj)%*%t(InvSig)%*%Ahati%*%InvSig%*%qsj
    SSEHat[i]<-SSEHat[i]+conti
  }
}

x<-seq(0,max(SSE),length=B)
plot(ecdf(SSE))
lines(x,pchisq(x,3*(n-1)),col=2)

x<-seq(0,max(SSEHat),length=B)
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
n<-50
qs<-ruars(n,rcayley,kappa=kap,space='Q4')

shat<-mean(qs)

rs<-rot.dist(qs,shat,method='intrinsic')
rs<-rcayley(1000,kappa=kap)
a<-mean(cos(rs/2)^2)-mean(cos(rs/2))^2
#b<-mean(sin(rs/2)^2)/3
InvSig<-diag(c(1/sqrt(a),rep(1/sqrt(b),3)))


Ahati<-diag(1,4)-t(shat)%*%shat
SSEFull<-0
SSERed<-rep(0,n)

for(j in 1:n){
  qsj<-matrix(qs[j,],4,1)
  SSEFull<-SSEFull+t(qsj)%*%t(InvSig)%*%Ahati%*%InvSig%*%qsj
  
  shat1<-mean(qs[-j,])
  Ahatij<-diag(1,4)-t(shat1)%*%shat1
  #InvSigj<-InvSig
  rsj<-rot.dist(qs[-j,],shat1,method='intrinsic')
  aj<-mean(cos(rsj/2)^2)-mean(cos(rsj/2))^2
  bj<-mean(sin(rsj/2)^2)/3
  InvSigj<-diag(c(1/sqrt(aj),rep(1/sqrt(bj),3)))
  
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
