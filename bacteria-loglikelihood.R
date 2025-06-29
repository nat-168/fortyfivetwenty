#Perform Log Likelihood
bacterialoglik_helper<-\(theta,v,n,x){
  stopifnot(length(v)==length(n),length(n)==length(x))
  k<-length(n)
  
  #∑log(ni choose xi)
  #pt1 <- sum(log(choose(n,x)))
  
  #∑-xi*theta*vi
  pt2 <- sum(-x*theta*v)
  
  #∑(ni-xi)log(1-e^-thetavi)
  pt3 <- sum((n-x)*log(1-exp(-theta*v)))
  
  return(pt2+pt3)
}

#Vectorize
bacterialoglik <- \(theta,v,n,x){
  sapply(theta, \(thet)bacterialoglik_helper(thet,v,n,x))
}

v<-c(1/100, 1/10, 1)
n<-c(5,5,2)
x<-c(4,1,0)

plot(bacterialoglik(161:180/10,v,n,x))

#unvectorized original tests:
#bacterialoglik(1,v,n,xx)
#plot(sapply(0:500/10,\(x)bacterialoglik(x,v,n,xx)))
#plot(sapply(101:200/10,\(x)bacterialoglik(x,v,n,xx)))
#plot(sapply(161:180/10,\(x)bacterialoglik(x,v,n,xx)))
