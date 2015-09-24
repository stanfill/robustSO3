#This code creates the plots of the Hn distribution with the theoretical 
#F limiting distribution
library(rotations)
kappa<-c(.2,.5,1,5,Inf)
n<-2500
HnSeq <- seq(0,4,length=n)
Caylabs <- expression(F["3,3(n-2)"],kappa=="5.0",kappa=="1.0",kappa==0.5,kappa==0.2)
vMiseslabs <- expression(F["1,1(n-2)"],kappa=="5.0",kappa=="1.0",kappa==0.5,kappa==0.2)


#Intrinsic test statistic
resDFCayInt<-data.frame(Hn=HnSeq,kappa=rep(kappa,each=n),ECDF=0)
resDFMisesInt<-data.frame(Hn=HnSeq,kappa=rep(kappa,each=n),ECDF=0)

for(i in 1:length(kappa)){
  
  row <- (i-1)*n+1
  
  if(kappa[i]<Inf){
    Qs<-ruars(n,rcayley,kappa=kappa[i])
    Hn<-discord(Qs,type='i')
    ecdfHN<-ecdf(Hn)
    resDFCayInt$ECDF[c(row:(row+(n-1)))]<-ecdfHN(HnSeq)

    Qs<-ruars(n,rvmises,kappa=kappa[i])
    Hn<-discord(Qs,type='i')
    ecdfHN<-ecdf(Hn)
    resDFMisesInt$ECDF[c(row:(row+(n-1)))]<-ecdfHN(HnSeq)
    
  }else{
    resDFCayInt$ECDF[c(row:(row+(n-1)))] <- pf(HnSeq,3,3*(n-2))
    
    resDFMisesInt$ECDF[c(row:(row+(n-1)))] <- pf(HnSeq,1,(n-2))
  }
}

resDFCayInt$kappa<-factor(resDFCayInt$kappa,levels=rev(kappa))
resDFMisesInt$kappa<-factor(resDFMisesInt$kappa,levels=rev(kappa))

qplot(Hn,ECDF,data=resDFCayInt,group=kappa,colour=kappa,geom='line',lwd=I(1.2))+
  ylab(expression(F[n](x)))+xlab(expression(x))+xlim(c(0,4))+theme_bw()+
  coord_fixed(4)+scale_color_grey(labels=Caylabs,name="")
#ggsave("C:/Users/Sta36z/Dropbox/SO3_Papers/OutlierID/Figures/CayleyHnIntrinsicECDF.pdf",height=5,width=5)

qplot(Hn,ECDF,data=resDFMisesInt,group=kappa,colour=kappa,geom='line',lwd=I(1.2))+
  ylab(expression(F[n](x)))+xlab(expression(x))+xlim(c(0,4))+theme_bw()+
  coord_fixed(4)+scale_color_grey(labels=vMiseslabs,name="")
#ggsave("C:/Users/Sta36z/Dropbox/SO3_Papers/OutlierID/Figures/vonMisesIntrinsicHnECDF.pdf",height=5,width=5)

#########################
#Extrinsic test statistic
resDFCayExt<-data.frame(Hn=HnSeq,kappa=rep(kappa,each=n),ECDF=0)
resDFMisesExt<-data.frame(Hn=HnSeq,kappa=rep(kappa,each=n),ECDF=0)

for(i in 1:length(kappa)){
  
  row <- (i-1)*n+1
  
  if(kappa[i]<Inf){
    Qs<-ruars(n,rcayley,kappa=kappa[i])
    Hn<-discord(Qs,type='e')
    ecdfHN<-ecdf(Hn)
    resDFCayExt$ECDF[c(row:(row+(n-1)))]<-ecdfHN(HnSeq)

    Qs<-ruars(n,rvmises,kappa=kappa[i])
    Hn<-discord(Qs,type='e')
    ecdfHN<-ecdf(Hn)
    resDFMisesExt$ECDF[c(row:(row+(n-1)))]<-ecdfHN(HnSeq)

    
  }else{
    resDFCayExt$ECDF[c(row:(row+(n-1)))] <- pf(HnSeq,3,3*(n-2))
    
    resDFMisesExt$ECDF[c(row:(row+(n-1)))] <- pf(HnSeq,1,(n-2))
  }
}

resDFCayExt$kappa<-factor(resDFCayExt$kappa,levels=rev(kappa))
resDFMisesExt$kappa<-factor(resDFMisesExt$kappa,levels=rev(kappa))

qplot(Hn,ECDF,data=resDFCayExt,group=kappa,colour=kappa,geom='line',lwd=I(1.2))+
  ylab(expression(F[n](x)))+xlab(expression(x))+xlim(c(0,4))+theme_bw()+coord_fixed(4)+
  scale_color_grey(labels=Caylabs,name="")
#ggsave("C:/Users/Sta36z/Dropbox/SO3_Papers/OutlierID/Figures/CayleyHnExtrinsicECDF.pdf",height=5,width=5)

qplot(Hn,ECDF,data=resDFMisesExt,group=kappa,colour=kappa,geom='line',lwd=I(1.2))+
  ylab(expression(F[n](x)))+xlab(expression(x))+xlim(c(0,4))+theme_bw()+
  coord_fixed(4)+scale_color_grey(labels=vMiseslabs)
#ggsave("C:/Users/Sta36z/Dropbox/SO3_Papers/OutlierID/Figures/vonMisesHnExtrinsicECDF.pdf",height=5,width=5)

