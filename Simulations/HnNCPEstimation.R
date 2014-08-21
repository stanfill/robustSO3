library(rotations)
source('~/robustSO3/Source_Code/robustFunctions.R')

#von Mises ncp
###No outliers therefore ncp=0
n <- 20
kap <- 25
B <- 1000
testStat <- rep(0,B)
for(i in 1:B){
  Rs <- ruars(n,rvmises,kappa=kap)
  rs <- matrix(rot.dist(Rs,id.SO3,method='intrinsic'),ncol=1)
  testStat[i] <- kap*t(rs)%*%rs
}
testStat <- sort(testStat)
plot(ecdf(testStat))
lines(testStat,pchisq(testStat,n),col=2)

########################
###One outlier giving ncp=k*r^2 where r is expected value of outlier
n <- 50
kap <- 25
B <- 1000
testStat <- rsOut <- rep(0,B)
rstar <- pi/8
Sstar <- genR(rstar)
ncp <- kap*(rstar^2)

for(i in 1:B){
  Rs <- ruarsCont(n,rvmises,kappa1=kap,p=1/n,Scont=Sstar)
  rs <- matrix(rot.dist(Rs,id.SO3,method='intrinsic'),ncol=1)
  testStat[i] <- kap*t(rs)%*%rs
}
testStat <- sort(testStat)
plot(ecdf(testStat))
lines(testStat,pchisq(testStat,n,ncp=ncp),col=2)

########################
###Hi statistic under Ho when its true
n <- 20
kap <- 25
B <- 1000
Hi <- rep(0,B)
for(i in 1:B){
  Rs <- ruars(n,rvmises,kappa=kap)
  Hi[i] <- discord(Rs, type='intrinsic',obs=n)
}
Hi <- sort(Hi)
plot(ecdf(Hi))
lines(Hi,pf(Hi,1,n-2),col=2)

########################
###Hi statistic under Ha when Ho should be rejected, diff mean alternative
n <- 50
kap <- 5
B <- 1000
Hi <- denom <- num <- rep(0,B)
rstar <- pi/8
Sstar <- genR(rstar)
ncpstar <- kap*rstar^2

for(i in 1:B){
  
  rs <- rvmises(n,kappa=kap)
  
  Rs <- genR(rs[-n])
  #denom[i] <- sum(rot.dist(Rs,id.SO3,method='intrinsic')^2)
  
  Outlier <- genR(rs[n],S=Sstar)
  #num[i] <- rot.dist(Outlier,id.SO3,method='intrinsic')^2
  
  #Hi[i] <- (num[i])/(denom[i]/(n-1))
  Rs <- as.SO3(rbind(Rs,Outlier))
  Hi[i] <- discord(Rs,type='i',obs=n)
  
}

#Compare numerator to non-central chi-square
num <- sort(num)
plot(ecdf(kap*num*3))
lines(3*kap*num,pchisq(3*kap*num,1,ncp=3*ncpstar),col=2)

#Compare their ratio to non-central F
scaler <- 1.5
Hi <- sort(Hi)
plot(ecdf(scaler*Hi))
lines(scaler*Hi,pf(scaler*Hi,1,n-1,ncp=scaler*ncpstar),col=2)

######
#Why doesn't numerator match theory? Too concentrated!  But why?!
#Because angles of rotation don't add like that.  If mis.angle(R1)=pi/2
#then mis.angle(R1,genR(pi/2))!=pi/2+pi/2=pi.  SO, just because the
#outlier is centered around R=genR(pi/2) doesn't mean E(R)=pi/2 where
n <- 1000
kap <- 100
rstar <- pi/4
rs <- sort(rvmises(n, kap) + rstar)
cStat <- kap*rs^2

plot(ecdf(cStat))
lines(cStat,pchisq(cStat,1,ncp=(kap*rstar^2)),col=2)

#
Rs<-genR(rs)
rsAbs <- rot.dist(Rs,id.SO3,method='intrinsic')
cStatSO3 <- kap*rsAbs^2
plot(ecdf(cStatSO3))
lines(cStatSO3,pchisq(cStatSO3,1,ncp=(kap*rstar^2)),col=2)

########################
###Hi statistic under Ha when Ho should be rejected, diff concentration alternative
###This is correct, when the outlier has same central orientation but different concentration tau
###then (tau/kappa)Hi~ F_1,n-1 central
n <- 25
kap <- 100
tau <- 5
B <- 1000
Hi <- rep(0,B)

for(i in 1:B){
  
  rs <- rvmises(n-1,kap)
  out <- rvmises(1,tau)
  Rs <- genR(c(rs,out))
  #rs2 <- mis.angle(Rs)^2
  #Hi[i] <- (rs2[n])/(sum(rs2[-n])/(n-1))
  Hi[i] <- discord(Rs, type='intrinsic', obs=n)
  
}

Hi <- sort(Hi)
scaleHi <- tau*Hi/kap
plot(ecdf(scaleHi))
lines(scaleHi,pf(scaleHi,1,n-2),col=2)


########################
###See if non-central chi^2 and F statistics work like I think they do
B<-1000
n<-10
mu<-pi/2
sigma<-100
chistat <- Fstat <- rep(0,B)

for(i in 1:B){  
  X1 <- rnorm(n,mu,sigma)
  chistat[i] <- sum((X1/sigma)^2)
  X2 <- rnorm(n)
  Fstat[i] <- chistat[i]/sum(X2^2)
}

chistat <- sort(chistat)
plot(ecdf(chistat))
lines(chistat,pchisq(chistat,n,ncp=n*(mu/sigma)^2),col=2)

Fstat <- sort(Fstat)
plot(ecdf(Fstat))
lines(Fstat,pf(Fstat,n,n,ncp=n*(mu/sigma)^2),col=2)
