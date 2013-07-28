library(rotations2)
library(scatterplot3d)

ruarsCont<-function(n,rangle,kappa,p,S=id.SO3,Scont){
	
	rs<-rangle(n,kappa=kappa)
	
	nCont<-floor(p*n)
	nNorm<-n-nCont
	
	RsCont<-genR(rs[1:nCont],Scont) #Simulated from the contaminated distribution
	RsNorm<-genR(rs[-c(1:nCont)],S)	#Simulate from the normal distribution
	
	Rs<-rbind(RsNorm,RsCont)
	return(as.SO3(Rs))
	
}

#Simulate rotations
n<-25
S2<-genR(pi/4)
Rs<-ruarsCont(n,rcayley,kappa=100,.25,Scont=S2)
Shat<-mean(Rs)
Stilde<-median(Rs)

#Transform them into vectors in R3, angle*axis 
rs<-angle(Rs)
us<-axis2(Rs)

#Add the angles and axis of the estiamtors
us<-rbind(us,axis2(Shat))
us<-rbind(us,axis2(Stilde))
rs<-c(rs,angle(Shat),angle(Stilde))
ws<-rs*us
ws<-rbind(ws,c(0,0,0))


#Red is mean, green is median, blue is truth
scatterplot3d(ws[,1],ws[,2],ws[,3],pch=19,color=c(rep(1,n),2,3,4))

scatterplot3d(ws[,1],ws[,2],ws[,3],pch=19,color=c(rep(1,n),2,3,4),type='h')