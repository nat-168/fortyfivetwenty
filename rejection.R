#Author: River Strumwasser
#rejection.R

#2a. simplereject(n, f, a, b, k, check.k=FALSE) n samples with pdf f(·)
simplereject <- \(n, f, a, b, k, check.k){
  
  #method
  U <- runif(n,a,b)
  V <- runif(n,0,k)
  
  Unew <- U[f(U)>V]
  #V <- V[f(U)>V]
  
  if(check.k) if(mean(f(Unew)>k)>0) {warning("k is less than f height") ; stop()}
  #plot(Unew,V)
  return(Unew)
}

#hist(simplereject(10000,dnorm,-5,5,1,T))

#2b. generalreject same but gen method
generalreject <- \(n, f, k, gpdf, gsample, check.k){
  
  #method
  X <- gsample(n)
  Y <- runif(n,0,k*gpdf(X))
  
  Xnew <- X[f(X)>Y]
  #Y <- Y[f(X>Y)]
  
  if(check.k) if(mean(f(Xnew)>k*gpdf(Xnew))>0) {warning("k*g is less than f height") ; stop()}
  #plot(Xnew,Y)
  return(Xnew)
}

#testing:

g <- \(x) 0.5*exp(-abs(x))

rlaplace <- \(n){
  ifelse(rbinom(n,1,0.5), rexp(n,1), -rexp(n,1))
}