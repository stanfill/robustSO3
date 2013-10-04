ruarsCont<-function(n,rangle,kappa,p,S=id.SO3,Scont){
  
  #n   		- sample size
  #rangle - angular distribution from which to simulate
  #kappa	- concentration parameter
  #p			- percent of n that will be contaminated
  #S			- central direction of normal data
  #Scont	-	central direction of contaminated data
  
  rs<-rangle(n,kappa=kappa)
  
  nCont<-floor(p*n)
  nNorm<-n-nCont
  
  RsCont<-genR(rs[1:nCont],Scont) #Simulated from the contaminated distribution
  RsNorm<-genR(rs[-c(1:nCont)],S)	#Simulate from the normal distribution
  
  Rs<-rbind(RsNorm,RsCont)
  return(as.SO3(Rs))
  
}

HnFun<-function(Qs,full=TRUE){
  #Compute the statistic proposed by FLW(?) that is a function of the largest eigenvalue
  #when observation i was removed
  
  if(ncol(Qs)==9){
    Qs<-Q4(Qs)
  }
  
  Tn<-t(Qs)%*%Qs
  n<-nrow(Qs)
  tau4n<-svd(Tn)$d[1]
  tau4n1<-rep(0,n)
  
  for(i in 1:n){
    Qsi<-Qs[-i,]
    Tn1<-t(Qsi)%*%Qsi
    tau4n1[i]<-svd(Tn1)$d[1]
  }
  
  if(full){
    Hn<- (n-2)*(1+tau4n1-tau4n)/(n-1-tau4n1)
  }else{
    Hn<- 1+tau4n1-tau4n
  }
  return(Hn)
}

MeanMove<-function(Qs,...){
  #Compute geodesic distance between full sample mean and mean when obs i is removed
  n<-nrow(Qs)
  Shat<-mean(Qs)
  ds<-rep(0,n)
  for(i in 1:n){
    Qsi<-as.Q4(Qs[-i,])
    Shati<-mean(Qsi)
    ds[i]<-dist(Shat,Shati,...)
  }
  return(ds)
}

trimMean<-function(Qs,a,discordFun,anneal=F){
  #Trim the most extreme a% based on the HnFun results
  #Qs - the sample
  #a - percent of sample to remove
  #discordFun - function to identify extreme observations, larger value more extreme obs
  #anneal - T/F, remove all at once (F) or one at a time (T)
  n<-nrow(Qs)
  nCut<-floor(min(max(0,n*a),n)) #remove at least 0, at most n
  
  if(nCut==0){
    return(list(Qs=Qs,Shat=mean(Qs)))
  }
  
  if(anneal){
    
    for(i in 1:nCut){
      Hn<-discordFun(Qs)
      Qs<-as.Q4(Qs[-which.max(Hn),])
    }
    return(list(Qs=Qs,Shat=mean(Qs)))
    
  }else{
    Hn<-discordFun(Qs)
    toCut<-which(order(Hn)>(n-nCut))
    tQs<-as.Q4(Qs[-toCut,])
    return(list(Qs=tQs,Shat=mean(tQs)))
  }
}

winzMean<-function(Rs,a,discordFun,anneal=F){
  #Project the most extreme a% observations closer
  #Rs - the sample
  #a - percent of sample to remove
  #discordFun - function to identify extreme observations, larger value more extreme obs
  #anneal - T/F, remove all at once (F) or one at a time (T)
  n<-nrow(Rs)
  nCut<-floor(min(max(0,n*a),n)) #project at least 0, at most n
  
  if(nCut==0){
    return(mean(Qs))
  }
  Shat<-mean(Rs)
  
  if(anneal){
    return(Shat)
  }else{
    
    crs<-dist(Rs,Shat,method='intrinsic')
    huberC<-sort(crs)[n-nCut] #There should be nCut observations further away than
                              #this value
    
    toWinz<-which(crs>huberC)
    badRs<-as.SO3(Rs[toWinz,]) #These are the observations to project closer to center
    us<-axis(badRs)
    winzRs<-SO3(us,rep(huberC,length(toWinz)))
    winzRs<-center(winzRs,t(Shat))
    wRs<-Rs
    wRs[toWinz,]<-winzRs
    
    return(list(Rs=wRs,Shat=mean(wRs)))
  }
  
}