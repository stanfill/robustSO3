#This gets the dataset in a unable condition

library(rotations)
library(plyr)
library(reshape2)
library(gridExtra)
library(xtable)

#g_legend will strip the legend and add it back so they can share one legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#load("/Users/stanfill/Dropbox/Rotation matrices/Melissa data/datasetnickel.RData")
load("C:/Users/Sta36z/Dropbox/Rotation matrices/Melissa data/datasetnickel.RData")

dat.out <- adply(data, .margins= c(1,3), function(x) {
  as.vector(x)
})
names(dat.out)[c(1,2)] <- c("rep", "location")
dat.out$xpos <- xpos[dat.out$location]
dat.out$ypos <- ypos[dat.out$location]
## convert data into a data frame including location and position information
#head(dat.out)

## are the matrices actually rotations?
checks <- adply(dat.out, .margins=1, function(x) {
  is.SO3(unlist(x[3:11]))
})
dat.out$check <- checks$V1
dat.out<-dat.out[dat.out$check==T,]

## compute the intrinsic and extrinsic test-statistic at each location
dat.discord <- dlply(dat.out, .(location), function(x) {
  res <- na.omit(x)
  res <- subset(res, check==TRUE)
  
  n <- nrow(res) 
  int <- ext  <- NULL
  if (n <= 2) {
    int <- ext  <- 0
  } else if (n > 0) {
    rots <- as.SO3(as.matrix(res[,3:11]))
    #print(rots)
    int <- max(discord(rots,type='i'))
    ext <- max(discord(rots,type='e'))
  }
  
  location <- as.numeric(as.character(unique(x$location)))
  return(list(location=location, n=n, int = int, ext = ext))
})

## find distances between estimators and angles to identity for each
loc.stats <- ldply(dat.discord, function(x) {  
  location <- as.numeric(as.character(unique(x$location)))
  if (x$n > 2)
    data.frame(location=x$location, n=x$n, 
               int=x$int, ext=x$ext,
               intp = min(1,x$n*pf(x$int,3,3*(x$n-2),lower.tail=F)), 
               extp = min(1,x$n*pf(x$ext,3,3*(x$n-2),lower.tail=F)))
})
loc.stats$xpos <- xpos[loc.stats$location]
loc.stats$ypos <- ypos[loc.stats$location]
