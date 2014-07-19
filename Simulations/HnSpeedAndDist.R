#Is it faster to compute Hn using the eigenvalue method, or using the 
#full/reduced model F test method?
#So far the eigenvalue method is faster, but that may just be because it is in C++
#Not sure it is worth coding the F test version in C++

source('~/robustSO3/Source_Code/robustFunctions.R')
Rcpp::sourceCpp('Source_Code/robustCpp.cpp')
library(rotations)
library(microbenchmark)

HnApprox<-function(Qs){
  n<-nrow(Qs)
  Qhat<-mean(Qs)
  Hnia<-rep(0,n)
  SSR<-sum(rot.dist(Qs,Qhat,method='extrinsic',p=2))
  for(i in 1:n){
    
    Qsi<-Qs[-i,]
    Qhati<-mean(Qsi)
    SSRi<-sum(rot.dist(Qsi,Qhati,method='extrinsic',p=2))
    Hnia[i] <- (n-2)*(SSR - SSRi)/(SSRi)
    
  }
  return(Hnia)
}

Qs<-ruars(20,rcayley,space='Q4')
microbenchmark(HnFun(Qs),HnApprox(Qs))

Rs<-ruars(20,rcayley)
microbenchmark(HnFun(Rs),HnApprox(Rs))

#HnFun is much faster

#################################
#Distribution of Hn with/without outliers

#Without outlier (should be F with 3,3*(n-2) dfs)
Qs<-ruars(50,rcayley,kappa=50)
Hn<-HnFun(Qs)
x<-seq(0,max(Hn),length=length(Hn))
plot(ecdf(Hn))
lines(x,pf(x,3,3*(length(Hn)-2)))

#With outliers
QsW<-ruarsCont(n=25,rangle=rcayley,kappa1=100,p=res$p[i],Scont=id.SO3,
              S=id.SO3,kappa2=1,space='Q4') 
HnW<-HnFun(QsW)
xW<-seq(0,max(HnW),length=length(HnW))
plot(ecdf(HnW))
lines(xW,pf(xW,3,3*(length(HnW)-2)))

#######
#Distribution of Hn en blocs of size t with/without outliers
t <- 1

#Without outlier (should be F with 3,3*(n-2) dfs)
Qs<-ruars(50,rcayley,kappa=100)
Hn<-HnBlocCpp(Qs,t)
x<-seq(0,max(Hn$Hn),length=length(Hn$Hn))
plot(ecdf(Hn$Hn))
lines(x,pf(x,3*t,3*(length(Hn$Hn)-1-t)),col=2)

#With outliers
QsW<-ruarsCont(n=25,rangle=rcayley,kappa1=100,p=res$p[i],Scont=id.SO3,
               S=id.SO3,kappa2=1,space='Q4') 
HnW<-HnBloc(QsW,t)
xW<-seq(0,max(HnW),length=length(HnW))
plot(ecdf(HnW))
lines(xW,pf(xW,3,3*(length(HnW)-1-t)))

###############
#Pretty plots of Hn ECDF for various kappa
library(rotations)
source('~/robustSO3/Source_Code/robustFunctions.R')
Rcpp::sourceCpp('Source_Code/robustCpp.cpp')

kappa<-c(1,2,5,10,Inf)
x<-seq(0,4,length=1000)
resDFCay<-data.frame(Hn=x,kappa=rep(kappa,each=length(x)),ECDF=0)
resDFMises<-data.frame(Hn=x,kappa=rep(kappa,each=length(x)),ECDF=0)


for(i in 1:length(kappa)){

  row <- (i-1)*1000+1
  
  if(kappa[i]<Inf){
    Qs<-ruars(1000,rcayley,kappa=kappa[i])
    Hn<-HnFun(Qs)
    ecdfHN<-ecdf(Hn)
    resDFCay$ECDF[c(row:(row+999))]<-ecdfHN(x)
    
    Qs<-ruars(1000,rvmises,kappa=kappa[i])
    Hn<-HnFun(Qs)
    ecdfHN<-ecdf(Hn)
    resDFMises$ECDF[c(row:(row+999))]<-ecdfHN(x)
    
    
  }else{
    resDFCay$ECDF[c(row:(row+999))] <- pf(x,3,3*1000)
    
    resDFMises$ECDF[c(row:(row+999))] <- pf(x,1,1000)
  }
}

resDFCay$kappa<-as.factor(resDFCay$kappa)
resDFMises$kappa<-as.factor(resDFMises$kappa)

resDFCay$kappa<-factor(resDFCay$kappa,levels=rev(levels(resDFCay$kappa)))
resDFMises$kappa<-factor(resDFMises$kappa,levels=rev(levels(resDFMises$kappa)))

qplot(Hn,ECDF,data=resDFCay,group=kappa,colour=kappa,xlab=expression(H[n]),
  ylab=expression(F[n](x)),geom='line',lwd=I(2),main='Cayley')+theme_bw()+
  scale_color_grey()+coord_fixed(4)

qplot(Hn,ECDF,data=resDFMises,group=kappa,colour=kappa,xlab=expression(H[n]),
  ylab=expression(F[n](x)),geom='line',lwd=I(2),main='von Mises')+theme_bw()+
  scale_color_grey()+coord_fixed(4)

