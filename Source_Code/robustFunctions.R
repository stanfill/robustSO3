ruarsCont<-function(n,rangle,kappa1,kappa2=kappa1,p,S=id.SO3,Scont,space='SO3'){
  
  #n   		- sample size
  #rangle - angular distribution from which to simulate
  #kappa1	- concentration parameter for F data
  #kappa2 - concentration for contaminated data
  #p			- percent of n that will be contaminated
  #S			- central direction of normal data
  #Scont	-	central direction of contaminated data
  #space  - SO3 (default) or quaternions("Q4")
    
  nCont<-floor(p*n)

  nNorm <- n-nCont

  rsNorm<-rangle(nNorm,kappa=kappa1)
  RsNorm<-genR(rsNorm,S)	#Simulate from the normal distribution
  
  if(nCont>0){
    rsCont<-rangle(nCont,kappa=kappa2)
    RsCont<-genR(rsCont,Scont) #Simulated from the contaminated distribution
    Rs<-as.SO3(rbind(RsNorm,RsCont))
  }else{
    Rs<-RsNorm
  }
  
  if(space=='Q4')
    Rs<-as.Q4(Rs)
  
  return(Rs)
  
}

#HnFun<-function(Qs,full=TRUE){
  #Compute the statistic proposed by FLW(?) that is a function of the largest eigenvalue
  #when observation i was removed
  #Written for quaternions, so if SO3 is given, make them quaternions
  
  #if(class(Qs)=="SO3")
  #  Qs<-as.Q4(Qs)
  
  #if(full){
  #  Hn<- as.vector(HnCpp(Qs))
  #}else{
  #  Hn<- as.vector(HnCpp(Qs))
  #}
  #return(Hn)
#}

HnBloc<-function(Qs,t){
  #Compute the Hn statistic when each possible set of t observations is deleted
  #warning("This function does not compute the right statistic, see Figue... and Gome (2005)")
  n<-nrow(Qs)
  groups <- combn(n,t)
  
  Qhat<-mean(Qs)
  Hnia<-rep(0,n)
  SSR<-sum(rot.dist(Qs,Qhat,method='extrinsic',p=2))
  
  for(i in 1:ncol(groups)){
    Qsi <- Qs[-groups[,i],]
    Qhati<-mean(Qsi)
    SSRi<-sum(rot.dist(Qsi,Qhati,method='extrinsic',p=2))
    Hnia[i] <- ((SSR - SSRi)/t)/((SSRi)/(n-t-1))
  }
  
  return(list(groups=groups,Hn=Hnia))
  
}

HnBlocCpp<-function(Qs,t){
  #Compute the Hn statistic when each possible set of t observations is deleted
  Qs<-as.Q4(Qs)
  n<-nrow(Qs)
  groups <- combn(n,t)
  
  Hnia <- HnCppBloc(Qs,groups)
  
  return(list(groups=groups,Hn=Hnia))
  
}

trimMean<-function(Qs,a,method='once',type='intrinsic',...){
  #Trim the most extreme a% based on the HnFun results
  #Qs - the sample
  #a - percent of sample to remove
  #method - method to determine which observations are outlying.  Options are:
  #                 'once' = compute Hn statics once and remove the out
  #                 'en bloc'  = compute Hn statistics for all blocks of size n*a
  #                 'anneal' = compute Hn, remove most outlying, recompute Hn remove most outlying...
  #... - additional arguements passed to mean function
  if(class(Qs)!="Q4")
    Qs<-as.Q4(Qs)
  
  
  n<-nrow(Qs)
  nCut<-floor(min(max(0,n*a),n)) #remove at least 0, at most n
  
  #Written for quaternions so change to quaternions if given matrices
  if(class(Qs)=="SO3")
    Qs<-Q4(Qs)
  
  if(nCut==0){
    #if(only){
    #  return(mean(Qs,...))
    #}else{
      return(list(Qs=Qs,Shat=mean(Qs,...)))
    #}
  }
  
  method <- try(match.arg(method,c("once","en bloc","anneal")))
  
  if(class(method)=='try-error')
    stop("Incorrect 'method' choice.")
  
  if(method=="anneal"){
    
    for(i in 1:nCut){
      Hn<-discord(Qs,type)
      Qs<-Qs[-which.max(Hn),]
    }
    #if(only){
    #  return(mean(Qs,...))
    #}else{
      return(list(Qs=Qs,Shat=mean(Qs,...)))
    #}
    
  }else if(method=='once'){
    Hn<-discord(Qs,type)
    toCut<-which(order(Hn)>(n-nCut))
    tQs<-Qs[-toCut,]
    #if(only){
    #  return(mean(tQs,...))
    #}else{
      return(list(Qs=tQs,Shat=mean(tQs,...)))
    #}
  }else{
    
    HnB <- HnBlocCpp(Qs,nCut)
    iToCut <- which.max(HnB$Hn)
    tQs <- Qs[-HnB$groups[,iToCut],]
    #if(only){
    # return(mean(tQs,...))
    #}else{
      return(list(Qs=tQs,Shat=mean(tQs,...)))
    #}
    
  }
  
}

bootSE<-function(samp,est,m,origEst,weight=F,...){
  
  #samp - the sample
  #est - estimator to estimate the SE of
  #m - number of replicates to use to estimate SE
  #origEst - the sample estimate
  #weight - if its the weighted mean then weights need to be updated
  
  n<-nrow(samp)
  SEhat<-0
  
  if(weight){
    
    for(i in 1:m){
      sampi<-samp[sample(1:n,replace=T),]
      Hns<-HnFun(sampi)
      esti<-est(sampi,w=(1/Hns),...)
      SEhat<-SEhat+rot.dist(origEst,esti)^2
    }
    
  }else{
    for(i in 1:m){
      sampi<-samp[sample(1:n,replace=T),]
      esti<-est(sampi,...)
      SEhat<-SEhat+rot.dist(origEst,esti)^2
    }
  }
  
  return(SEhat/n)
  
}

trimMeanOld<-function(Qs,a,discordFun,anneal=F,...){
  #Trim the most extreme a% based on the HnFun results
  #Qs - the sample
  #a - percent of sample to remove
  #discordFun - function to identify extreme observations, larger value more extreme obs
  #anneal - T/F, remove all at once (F) or one at a time (T)
  #... - additional arguements passed to mean function
  n<-nrow(Qs)
  nCut<-floor(min(max(0,n*a),n)) #remove at least 0, at most n
  
  #Written for quaternions so change to quaternions if given matrices
  if(class(Qs)=="SO3")
    Qs<-Q4(Qs)
  
  if(nCut==0){
    return(list(Qs=Qs,Shat=mean(Qs)))
  }
  
  if(anneal){
    
    for(i in 1:nCut){
      Hn<-discordFun(Qs)
      Qs<-as.Q4(Qs[-which.max(Hn),])
    }
    return(list(Qs=Qs,Shat=mean(Qs,...)))
    
  }else{
    Hn<-discordFun(Qs)
    toCut<-which(order(Hn)>(n-nCut))
    tQs<-as.Q4(Qs[-toCut,])
    return(list(Qs=tQs,Shat=mean(tQs,...)))
  }
}


# require(hypergeo)
# bipolarWatson<-function(r,kappa,Haar=T){
#   
#   kernal <- exp(kappa*(cos(r)+1)/2)
# 
#   const <- 2*pi*genhypergeo(z=kappa,U=0.5,L=2)
#   
#   if(Haar){
#     return(kernal/const)
#   }
#   
#   return(kernal*(1-cos(r))/(const))
#   
# }
#integrate(bipolarWatson,-pi,pi,kappa=1,Haar=F)


##################
#To test for outliers, want to compare the largest Hi value
#to the distribution for the largest value, not to the F

dorderF <- function(x,n,k=n,df1,df2,ncp=0){
  #x - the value at which to evaluate the pdf
  #n - the sample size
  #k - the kth observation, k in {1,...,n}
  #df1 - numerator df
  #df2 - denominator df
  #ncp - non-centrality parameter
  
  if(k<1 || k>n){
    stop("k is not in the range 1,...,n")
  }
  
  if(k==n){
    
    p1 <- n
    p3 <- 1
    
  }else{
    
    p1 <- factorial(n)/(factorial(k-1)*factorial(n-k))
    p3 <- (1-pf(x,df1,df2,ncp))^(n-k)
    
  }
  
  p2 <- pf(x,df1,df2,ncp)^(k-1)
  p4 <- df(x,df1,df2,ncp)
  return(p1*p2*p3*p4)
}

porderF <- function(x,n,k=n,df1,df2,ncp=0,lower.tail=TRUE){
  #Integrate dorderF from 0 to x to estimate F(x)=P(X<=x)
  
  lt <- rep(0,length(x))
  
  for(i in 1:length(x)){
    lt[i] <- integrate(dorderF,0,x[i],n=n,k=k,df1=df1,df2=df2,ncp=ncp)$value
  }
  if(lower.tail){
    return(lt)
  }else{
    return(1-lt)
  }
  
}

helpQOF <- function(x,q,n,k=n,df1,df2,ncp=0){
  
  Fx <- porderF(x,n,k,df1,df2,ncp)
  return((Fx-q)^2)
  
}

qorderF <- function(q,n,k=n,df1,df2,ncp=0){
  #Use optim to find the qth quantile of order F
  
  sol <- optimize(f=helpQOF,interval=c(0,100),q=q,n=n,k=k,df1=df1,df2=df2,ncp=ncp)
  return(sol$minimum)
  
}
