#this code allows for simulation from contaminated distributions

library(rotations2)

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

rangle<-rcayley
n<-50
kappa<-100
p<-.25
S<-id.SO3
Scont<-genR(pi/2)

Rs<-ruarsCont(n,rangle,kappa,p,S,Scont)

plot(Rs,center=id.SO3,col=2)
plot(Rs,center=median(Rs),show_estimates=c("proj.mean","proj.median"),col=2)

#Chang method
plot(Rs,center=median(Rs),show_estimates=c("proj.mean","proj.median"),col=1,
		 mean_regions='chang',median_regions='chang',alp=.05)

plot(Rs,center=median(Rs),show_estimates=c("proj.mean","proj.median"),col=2,
		 mean_regions='chang',median_regions='chang',alp=.1)

plot(Rs,center=median(Rs),show_estimates=c("proj.mean","proj.median"),col=3,
		 mean_regions='chang',median_regions='chang',alp=.1)

#Zhang method
plot(Rs,center=median(Rs),show_estimates=c("proj.mean","proj.median"),col=2,
		 mean_regions='zhang',median_regions='zhang',alp=.1)
