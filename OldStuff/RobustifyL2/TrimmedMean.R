HnFun<-function(Qs){
  #Compute the statistic proposed by FLW(?) that is a function of the largest eigenvalue
  #when observation i was removed
  
  Tn<-t(Qs)%*%Qs
  n<-nrow(Qs)
  tau4n<-svd(Tn)$d[1]
  tau4n1<-rep(0,n)
  
  for(i in 1:n){
    Qsi<-Qs[-i,]
    Tn1<-t(Qsi)%*%Qsi
    tau4n1[i]<-svd(Tn1)$d[1]
  }
  Hn<-(n-2)*(1+tau4n1-tau4n)/(n-1-tau4n1)
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

Qs<-ruars(50,rvmises,kappa=5,space='Q4')
#plot(Rs)
#Qs<-Q4(Rs)
hns<-HnFun(Qs)
ds<-MeanMove(Qs,type='intrinsic',p=2)
plot(hns,ds,pch=19)

system.time(HnFun(Qs))  #This is faster
system.time(MeanMove(Qs))

###################
#Compute the trimmed mean based on HnFun

trimMean<-function(Qs,a,discordFun,anneal=F){
  #Trim the most extreme a% based on the HnFun results
  #Qs - the sample
  #a - percent of sample to remove
  #discordFun - function to identify extreme observations, larger value more extreme obs
  #anneal - T/F, remove all at once (F) or one at a time (T)
  
  n<-nrow(Qs)
  nCut<-floor(min(max(0,n*a),n/2)) #remove atleast 0, atmost n/2
  
  if(nCut==0){
    return(mean(Qs))
  }
  
  if(anneal){
    
    for(i in 1:nCut){
      Hn<-discordFun(Qs)
      Qs<-as.Q4(Qs[-which.max(Hn),])
    }
    return(mean(Qs))
    
  }else{
    Hn<-discordFun(Qs)
    toCut<-which(order(Hn)>(n-nCut))
    return(mean(as.Q4(Qs[-toCut,])))
  }
}


#Compare trimmed means based on Hn statistic, annealed and not
QsFish<-ruars(50,rfisher,kappa=10,space='Q4')
dist(mean(QsFish))
dist(trimMean(QsFish,a=.05,discordFun=HnFun))
dist(trimMean(QsFish,a=.05,discordFun=HnFun,anneal=T))
plot(SO3(QsFish))

QsCay<-ruars(50,rcayley,kappa=10,space='Q4')
dist(mean(QsCay))
dist(trimMean(QsCay,a=.05,discordFun=HnFun))
dist(trimMean(QsCay,a=.05,discordFun=HnFun,anneal=T))
plot(SO3(QsCay))

QsMis<-ruars(50,rvmises,kappa=5,space='Q4')
dist(mean(QsMis))
dist(trimMean(QsMis,a=.05,discordFun=HnFun))
dist(trimMean(QsMis,a=.05,discordFun=HnFun,anneal=T))
plot(SO3(QsMis))

#Compare trimmed means based on distance from mean
meanFN<-function(x){return(dist(x,mean(x)))}

QsFish<-ruars(50,rfisher,kappa=10,space='Q4')
dist(mean(QsFish))
dist(trimMean(QsFish,a=.05,discordFun=meanFN))
dist(trimMean(QsFish,a=.05,discordFun=meanFN,anneal=T))
plot(SO3(QsFish))

QsCay<-ruars(50,rcayley,kappa=10,space='Q4')
dist(mean(QsCay))
dist(trimMean(QsCay,a=.05,discordFun=meanFN))
dist(trimMean(QsCay,a=.05,discordFun=meanFN,anneal=T))
plot(SO3(QsCay))

QsMis<-ruars(50,rvmises,kappa=5,space='Q4')
dist(mean(QsMis))
dist(trimMean(QsMis,a=.05,discordFun=meanFN))
dist(trimMean(QsMis,a=.05,discordFun=meanFN,anneal=T))
plot(SO3(QsMis))

