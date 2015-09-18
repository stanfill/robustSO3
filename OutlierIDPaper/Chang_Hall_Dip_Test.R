library(rotations)
library(ks)
library(diptest)
source('Source_Code/robustFunctions.R')
Rcpp::sourceCpp('Source_Code/NonParaBootstrapCpp.cpp')

inv_gamma_betaDist <- function(b,dhat){
  return((dhat-2^(4*b-1)*(b-1)*beta(b,b)^2)^2)
}


inv_gamma_scaledtDist <- function(b,dhat){
  return((dhat-2*b*beta(b-.5,.5)^2)^2)
}

#######

#Generate data with/without outlier
rstar <- pi/2
kap <- 50
n <- 10
rangle <- rcayley
QsOut <- ruarsCont(n,rangle,kappa1=kap,p=1/n,Scont=genR(rstar),space='Q4')

#Bootstrap data to get distribution of Hn
Hsstar <- HnNonParaBootCpp(QsOut,m=1000,type=1)

#####################
#Cheng and Hall calibrated dip test

#1 - compute dip test statistic based on distribution of Hsstar
delta <- dip(Hsstar)

#2 - Estimate density of Hsstar using Gaussian kernel estimate and find largest model
hOpt <- hns(Hsstar, deriv.order=0) #optimal bandwidth
fhat <- kdde(Hsstar[1,], h=hOpt, deriv.order=0)
fhatMode <- which.max(fhat$estimate)

#3 - Estimate second derivative of density of Hsstar using Gaussian kernel estimate
h2Opt <- hns(Hsstar,deriv.order=2)
fhat2nd <- kdde(Hsstar[1,], h=h2Opt, deriv.order=2, eval.points=fhat$eval.points)

#4 - Compute dhat
dhat <- abs(fhat2nd$estimate[fhatMode])/(fhat$estimate[fhatMode]^3)

#5 - Estimate beta from gamma^{-1}(dhat) and generate Xstar
if(dhat<2*pi){
  #Sample from beta distribution
  betaRes <- optimize(f=inv_gamma_betaDist,interval=c(1.001,254),dhat=dhat)
  betaHat <- betaRes$minimum
  Xstar <- matrix(rbeta(500*length(Hsstar),betaHat,betaHat),nrow=500)
}else{
  #Sample from rescaled t-distribution
  betaRes <- optimize(f=inv_gamma_scaledtDist,interval=c(0.05,100),dhat=dhat)
  betaHat <- betaRes$minimum
  eta <- 2*betaHat-1
  Xstar <- matrix((1/sqrt(eta))*rt(500*length(Hsstar),df=eta),nrow=500)
}

#6 - Generate bootstrap distribution of delta based on Xstar
deltaStar <- apply(Xstar,1,dip)

#7 - Reject H0?
delta>quantile(deltaStar,(1-0.05))

#####################
#Compress into a function


calib_dip.test <- function(Hsstar,m=500,type=1,alpha=0.05){
  ##Hsstar - quaternion observations
  ##m - number of bootstrap samples of the data to take
  ##type - type of test statistic to use: 1=intrinisic, 2=extrinsic
  ##alpha - level of the test
  ###Returns 0/1 corresponding to "fail to reject"/reject

  #1 - compute dip test statistic based on distribution of Hsstar
  delta <- dip(Hsstar)
  
  #2 - Estimate density of Hsstar using Gaussian kernel estimate and find largest model
  hOpt <- hns(Hsstar, deriv.order=0) #optimal bandwidth
  fhat <- kdde(Hsstar[1,], h=hOpt, deriv.order=0)
  fhatMode <- which.max(fhat$estimate)
  
  #3 - Estimate second derivative of density of Hsstar using Gaussian kernel estimate
  h2Opt <- hns(Hsstar,deriv.order=2)
  fhat2nd <- kdde(Hsstar[1,], h=h2Opt, deriv.order=2, eval.points=fhat$eval.points)
  
  #4 - Compute dhat
  dhat <- abs(fhat2nd$estimate[fhatMode])/(fhat$estimate[fhatMode]^3)
  
  #5 - Estimate beta from gamma^{-1}(dhat) and generate from correct dist
  if(dhat<2*pi){
    #Sample from beta distribution
    betaRes <- optimize(f=inv_gamma_betaDist,interval=c(1.001,254),dhat=dhat)
    betaHat <- betaRes$minimum
    Xstar <- matrix(rbeta(500*length(Hsstar),betaHat,betaHat),nrow=500)
  }else{
    #Sample from rescaled t-distribution
    betaRes <- optimize(f=inv_gamma_scaledtDist,interval=c(0.5001,1000),dhat=dhat)
    betaHat <- betaRes$minimum
    eta <- 2*betaHat-1
    Xstar <- matrix((1/sqrt(eta))*rt(500*length(Hsstar),df=eta),nrow=500)
  }
  
  
  #6 - Generate bootstrap distribution of delta based on Xstar
  deltaStar <- apply(Xstar,1,dip)
  
  #7 - Reject H0?
  return(delta>quantile(deltaStar,(1-alpha)))
  
}

#################
#Compare calibrated and uncalibrated dip test

rstar <- c(0,pi/8,pi/4,pi/2)
kap <- 50
n <- 50
rangle <- rcayley
B <- 250
calib <- uncalib <- matrix(0,B,4)
alpha <- 0.1

for(j in 1:4){
  for(i in 1:B){
    
    #Generate data with/without outlier
    QsOut <- ruarsCont(n,rangle,kappa1=kap,p=1/n,Scont=genR(rstar[j]),space='Q4')
    
    #Bootstrap data to get distribution of Hn
    Hsstar <- HnNonParaBootCpp(QsOut,m=500,type=1)
    
    calib[i,j] <- calib_dip.test(Hsstar,alpha=alpha)
    uncalib[i,j] <- dip.test(Hsstar)$p.value<alpha
    
  }
}
#Estimated P(Reject H0)
colSums(calib)/B
colSums(uncalib)/B









