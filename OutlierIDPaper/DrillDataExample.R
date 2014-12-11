library(rotations)
data(drill)

nSubj <- max(drill$Subject)

subj1 <- na.omit(subset(drill,Subject==2))
his <- hes <- matrix(0,6,3)

for(i in 1:6){
  Wdat <- as.Q4(subset(subj1,Position==i & Joint=="Wrist")[,-c(1:4)])
  his[i,1] <- max(discord(Wdat,type='i'))
  hes[i,1] <- max(discord(Wdat,type='e'))
  
  Edat <- as.Q4(subset(subj1,Position==i & Joint=="Elbow")[,-c(1:4)])
  his[i,2] <- max(discord(Edat,type='i'))
  hes[i,2] <- max(discord(Edat,type='e'))
  
  Sdat <- as.Q4(subset(subj1,Position==i & Joint=="Shoulder")[,-c(1:4)])
  his[i,3] <- max(discord(Sdat,type='i'))
  hes[i,3] <- max(discord(Sdat,type='e'))
}

(intp  <- 5*pf(his,3,9,lower.tail=F))
(extp  <- 5*pf(hes,3,9,lower.tail=F))

