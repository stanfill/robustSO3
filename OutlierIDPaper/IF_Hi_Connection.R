library(rotations)

influence <- function(Rs,const,type='intrinsic'){
  #Compute the influence for each observation
  #c = E[1+2cos(r)]
  #d = E{[1+3cos(r)]/[12*sqrt(1-cos(r))]}
  
  rs <- mis.angle(Rs)
  if(type=='extrinsic'){
    
    ifval <- 3*sin(rs)/const
    
  }else if(type=='intrinsic'){
    
    top <- sin(rs)
    bottom <- 2*const*sqrt(1-cos(rs))
    ifval <- top/bottom
    
  }else{
    ifval <- rep(0,length(rs))
  }
  
  return(ifval)
}


Rs<-ruars(200,rcayley,kappa=1)
rsHat <- mis.angle(Rs)
#Estimate chat
chat <- mean(1+2*cos(rsHat))

HiE <- discord(Rs,type='extrinsic')
ifE <- influence(Rs,chat,type='extrinsic')
plot(HiE,ifE)

#Estimate dhat
dhat <- 1+3*cos(rsHat)
dhat <- mean(dhat/12*sqrt(1-cos(rsHat)))

HiI <- discord(Rs,type='intrinsic')
ifI <- influence(Rs,dhat,type='intrinsic')
plot(HiI,ifI)
