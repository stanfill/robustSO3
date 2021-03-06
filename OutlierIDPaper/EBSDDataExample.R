source("OutlierIDPaper/PrepareEBSDDataSet.R")

#Grain map based on median of all 14 scans
possibles<-which(loc.stats$extp<(0.01/nrow(loc.stats)))
d <- ggplot(loc.stats, aes(xpos, ypos, color=intp))
d2 <- d + geom_point(size=4) + scale_colour_gradient(expression(P(Reject~H[0])), low="grey99", high="grey10", limits=c(0, 0.05)) + 
  theme_bw() + xlab("") + ylab("") + coord_equal() + scale_x_continuous(limits=c(0, 12.5), breaks=seq(0, 12.5, by=2.5), labels=expression(0*mu*m, 2.5*mu*m, 5*mu*m, 7.5*mu*m, 10*mu*m, 12.5*mu*m)) + 
  scale_y_continuous(limits=c(0, 10), breaks=seq(0, 10, by=2.5), labels=expression(0*mu*m, 2.5*mu*m, 5*mu*m, 7.5*mu*m, 10*mu*m)) + 
  geom_point(shape="o", colour="yellow", size=5, data=loc.stats[possibles,])  + 
  theme(plot.margin=unit(rep(0,4), "lines"))
d2

qplot(intp,extp,data=loc.stats)
which(loc.stats$intp<0.01)



possibles<-which(loc.stats$extp<(0.05/nrow(loc.stats)))
d <- ggplot(loc.stats, aes(xpos, ypos, color=dR))
d2 <- d + geom_point(size=6) + scale_colour_gradient(expression(P(Reject~H[0])), low="grey99", high="grey10", limits=c(0, pi)) + 
  theme_bw() + xlab("") + ylab("") + coord_equal() + scale_x_continuous(limits=c(0, 12.5), breaks=seq(0, 12.5, by=2.5), labels=expression(0*mu*m, 2.5*mu*m, 5*mu*m, 7.5*mu*m, 10*mu*m, 12.5*mu*m)) + 
  scale_y_continuous(limits=c(0, 10), breaks=seq(0, 10, by=2.5), labels=expression(0*mu*m, 2.5*mu*m, 5*mu*m, 7.5*mu*m, 10*mu*m)) + 
  geom_point(shape="o", colour="yellow", size=2, data=loc.stats[possibles,])  + 
  theme(plot.margin=unit(rep(0,4), "lines"))
d2

qplot(intp,extp,data=loc.stats)
which(loc.stats$intp<0.01)


##################
#Zoom in on a particular ares
xlow <- 4.5; xhigh <- 6; ylow <- 6; yhigh <- 7
areaOI <- subset(loc.stats,xpos<xhigh & xpos>xlow & ypos>ylow & ypos<yhigh)

possibles<-which(areaOI$extp<(0.05/nrow(loc.stats)))
d <- ggplot(loc.stats, aes(xpos, ypos, color=dR))
d2 <- d + geom_point(size=4) + scale_colour_gradient(expression(P(Reject~H[0])), low="grey99", high="grey10", limits=c(0, pi)) + 
  theme_bw() + xlab("") + ylab("") + coord_equal() + scale_x_continuous(limits=c(0, 12.5), breaks=seq(0, 12.5, by=2.5), labels=expression(0*mu*m, 2.5*mu*m, 5*mu*m, 7.5*mu*m, 10*mu*m, 12.5*mu*m)) + 
  scale_y_continuous(limits=c(0, 10), breaks=seq(0, 10, by=2.5), labels=expression(0*mu*m, 2.5*mu*m, 5*mu*m, 7.5*mu*m, 10*mu*m)) + 
  geom_point(shape="o", colour="yellow", size=5, data=areaOI[possibles,])  + 
  theme(plot.margin=unit(rep(0,4), "lines"))
d2