ruarsCont<-function(n,rangle,kappa,p,S=id.SO3,Scont,space='SO3'){
  
  #n   		- sample size
  #rangle - angular distribution from which to simulate
  #kappa	- concentration parameter
  #p			- percent of n that will be contaminated
  #S			- central direction of normal data
  #Scont	-	central direction of contaminated data
  #space  - SO3 (default) or quaternions("Q4")
  
  rs<-rangle(n,kappa=kappa)
  
  nCont<-floor(p*n)
  nNorm<-n-nCont
  
  RsCont<-genR(rs[1:nCont],Scont) #Simulated from the contaminated distribution
  RsNorm<-genR(rs[-c(1:nCont)],S)	#Simulate from the normal distribution
  
  Rs<-as.SO3(rbind(RsNorm,RsCont))
  
  if(space=='Q4')
    Rs<-Q4(Rs)
  
  return(Rs)
  
}

HnFun<-function(Qs,full=TRUE){
  #Compute the statistic proposed by FLW(?) that is a function of the largest eigenvalue
  #when observation i was removed
  #Written for quaternions, so if SO3 is given, make them quaternions
  
  if(class(Qs)=="SO3")
    Qs<-Q4(Qs)
  
  if(full){
    Hn<- as.vector(HnCpp(Qs))
  }else{
    Hn<- as.vector(HnCpp(Qs))
  }
  return(Hn)
}

MeanMove<-function(Qs,...){
  #Compute geodesic distance between full sample mean and mean when obs i is removed
  
  #Written for quaternions so change to quaternions if given matrices
  if(class(Qs)=="SO3")
    Qs<-Q4(Qs)
  
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
  
  #Written for rotations, so if given quaternions make them rotations
  if(class(Rs)=="Q4")
    Rs<-SO3(Rs)
  
  if(nCut==0){
    return(mean(Rs))
  }
  Shat<-mean(Rs)
  
  if(anneal){
    return(Shat)
  }else{
    
    crs<-discordFun(Rs)
    huberC<-sort(crs)[n-nCut] #There should be nCut observations further away than
                              #this value
    
    toWinz<-which(crs>huberC)
    badRs<-as.SO3(Rs[toWinz,]) #These are the observations to project closer to center
    badRs<-center(badRs,Shat)
    us<-axis(badRs)
    winzRs<-SO3(us,rep(huberC,length(toWinz)))
    winzRs<-center(winzRs,t(Shat))
    wRs<-Rs
    wRs[toWinz,]<-winzRs
    
    return(list(Rs=wRs,Shat=mean(wRs)))
  }
  
}

HuberMean<-function(RS,c){
  #Find the multidimensional Huber estimator based on 
  #projected mean
  #Rs - the sample
  #c  - the value the influence function should not exceed
  
  if(class(Rs)=="Q4")
    Rs<-SO3(Rs)
  iters<-0
  
  while(iters<100){
  
    shat<-mean(Rs)
    rs<-dist(Rs,shat,method='intrinsic')          #Estimate r_i based on Shat
    #dhat<-mean(1+2*cos(rs))/3             #Estimate dhat
    #ifi<-sin(rs)/dhat                    #Evaluate infulence function
  
    toofar<-which(rs>(c+10e-5))
    numFar<-length(toofar)
  
    if(numFar==0){
      break
    }else{
      #print(toofar)
      #print(rs[toofar])
      toMove<-as.SO3(Rs[toofar,])
      toMove<-center(toMove,shat)
      usStar<-axis(toMove)
      Moved<-center(SO3(usStar,rep(c,numFar)),t(shat))
      Rs[toofar,]<-Moved
    }
    iters<-iters+1
    #Rs<-center(Rs,t(shat)) #add the center back
  }
  #print(iters)
  return(list(Rs=Rs,Shat=mean(Rs)))
}