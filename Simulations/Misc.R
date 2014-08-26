library(rotations)
rstar<-pi/2
Sstar<-genR(rstar)
kap<-500  

rs<-(rvmises(1000,kappa=kap))

rs2<-(mis.angle(genR(rs,S=Sstar)))

layout(1)
plot((rs),rs2)
abline(h=rstar,v=0);abline(mean(rs2),1);abline(mean(rs2),-1)

plot(kap*rs^2,kap*rs2^2)
abline(h=kap*rstar^2)

par(mfrow=c(1,2))
hist(rs,main=c(1/kap,var(rs)))
hist(rs2,main=c(1/(3*kap),var(rs2)))

rsStat<-sort(kap*rs^2)
plot(ecdf(rsStat));lines(rsStat,pchisq(rsStat,1),col=2)

sc <- .2
rs2Stat <- sort(kap*(rs2^2)/sc)
plot(ecdf(rs2Stat));lines(rs2Stat,pchisq(rs2Stat,1,ncp=kap*rstar^2/sc),col=2)
#
x<-rnorm(100)
xstar<-x+5
plot(x,xstar)
#
rseq <- seq(-1,1,length=100)
plot(rseq,acos(rseq),type='l')
abline(h=0,v=0)


####
kap<-5
rstar<-seq(pi/10,pi,length=20)
ratio <- rep(0,length(rstar))

for(i in 1:length(rstar)){
  Sstar <- genR(rstar[i])
  rs<-(rvmises(10000,kappa=kap))
  rs2<-mis.angle(genR(rs,S=Sstar))
  varrs <- apply(matrix(rs,nrow=100),2,var)
  varrs2 <- apply(matrix(rs2,nrow=100),2,var)
  ratio[i] <- mean(varrs2)/mean(varrs)
}
plot(rstar,ratio,main=mean(ratio))


############################
#Look at shape of dist of dorderF
rs <- seq(0,10,length=100)
plot(rs,dorderF(rs,n=20,df1=3,df2=3*(18)),type='l')
lines(rs,df(rs,3,3*18),type='l',col=2)

plot(rs,porderF(rs,n=20,df1=3,df2=3*(18)),type='l')
lines(rs,pf(rs,3,3*18),col=2)

#Make sure porderF is correct
maxF <- rep(0,500)
for(i in 1:length(maxF)){
  maxF[i] <- max(rf(10,5,10))  
}
maxF <- sort(maxF)
plot(ecdf(maxF))  
lines(maxF,porderF(maxF,n=10,df1=5,df2=10),col=2)

