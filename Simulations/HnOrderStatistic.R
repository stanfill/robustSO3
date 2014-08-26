##################
#To test for outliers, want to compare the largest Hi value
#to the distribution for the largest value, not to the F

dorderF <- function(x,n,k=n,df1,df2,ncp=0){
  #x - the value at which to evaluate the pdf
  #n - the sample size
  #k - the kth observation, k in {1,...,n}
  #df1 - numerator df
  #df2 - denominator df
  #ncp - non-centrality parameter
  
  if(k<1 || k>n){
    stop("k is not in the range 1,...,n")
  }
  
  p1 <- factorial(n)/(factorial(k-1)*factorial(n-k))
  p2 <- pf(x,df1,df2,ncp)^(k-1)
  p3 <- (1-pf(x,df1,df2,ncp))^(n-k)
  p4 <- df(x,df1,df2,ncp)
  return(p1*p2*p3*p4)
}

porderF <- function(x,n,k=n,df1,df2,ncp=0,lower.tail=TRUE){
  
  lt <- rep(0,length(x))
  for(i in 1:length(x)){
    lt[i] <- integrate(dorderF,0,x[i],n=n,k=k,df1=df1,df2=df2,ncp=ncp)$value
  }
  if(lower.tail){
    return(lt)
  }else{
    return(1-lt)
  }
  
}

############################
#Look at shape of dist of dorderF
rs <- seq(0,10,length=100)
plot(rs,dorderF(rs,n=20,df1=3,df2=3*(18)),type='l')
lines(rs,df(rs,3,3*18),type='l',col=2)

plot(rs,porderF(rs,n=20,df1=3,df2=3*(18)),type='l')
lines(rs,pf(rs,3,3*18),col=2)

######
#Compare dist of H_(n) to porderF


