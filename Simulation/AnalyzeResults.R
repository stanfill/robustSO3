library(ggplot2)

extrinsic<-read.csv("Simulation/Results/ExtrinsicRobustSimulations_1000.csv")[,-1]
intrinsic<-read.csv("Simulation/Results/IntrinsicRobustSimulations_1000.csv")[,-1]

extrinsic$Method<-'Extrinsic'
intrinsic$Method<-'Intrinsic'

final<-rbind(extrinsic,intrinsic)

#Compare contamination S's within each Method
qplot(Eps,value,data=intrinsic[intrinsic$Measure=='Bias',],colour=Estimator,group=Estimator,geom='line',size=I(1.25),xlab=expression(epsilon),ylab='Bias')+
  facet_grid(Sstar~n,labeller = label_parsed)#+coord_fixed(1/3)

qplot(Eps,value,data=extrinsic[extrinsic$Measure=='Bias',],colour=Estimator,group=Estimator,geom='line',size=I(1.25),xlab=expression(epsilon),ylab='Bias')+
  facet_grid(Sstar~n,labeller = label_parsed)#+coord_fixed(1/3)

#Compare methods for a given contamination S
qplot(Eps,value,data=final[final$Measure=='Bias' & final$Sstar=='pi/2',],colour=Estimator,group=Estimator,geom='line',size=I(1.25),xlab=expression(epsilon),ylab='Bias')+
  facet_grid(Method~n,labeller = label_parsed,scales="free_y")#+coord_fixed()

qplot(Eps,value,data=final[final$Measure=='Bias' & final$Sstar=='pi',],colour=Estimator,group=Estimator,geom='line',size=I(1.25),xlab=expression(epsilon),ylab='Bias')+
  facet_grid(Method~n,labeller = label_parsed,scales="free_y")#+coord_fixed()
