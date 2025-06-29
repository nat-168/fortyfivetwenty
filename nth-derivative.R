#nth-derivative.R by River Strumwasser

P <- \(x,n){
  if (n==0) {
    return(1)
  }
  if (n==1) {
    return(-x)
  }
  if (n>1) {
    return(-(n-1)*P(x,n-2)-x*P(x,n-1))
  }
}

phi <- \(x,n){
  return(P(x,n)*dnorm(x))
}

P6 <- \(x)(x^6-15*x^4+45*x^2-15)
#curves curve(phi(x,6),-5,5) and curve(P6(x)*dnorm(x),-5,5) are same!

psi <- \(x,n,mean,sd){
  return(phi((x-mean)/sd, n)/(sd^(n+1)))
}