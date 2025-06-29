#Author: River Strumwasser
#mixmode.R


#pdf, deriv pdf, 2nd deriv pdf:
mmf <- \(x,mu0,mu1,sig0,sig1,p){
  p*dnorm(x,mean=mu1,sd=sig1)+(1-p)*dnorm(x,mean=mu0,sd=sig0)
}

mmfp <- \(x,mu0,mu1,sig0,sig1,p){
  -p*((x-mu1)/sig1^2)*dnorm(x,mean=mu1,sd=sig1)-(1-p)*((x-mu0)/sig0^2)*dnorm(x,mean=mu0,sd=sig0)
}

mmfp2 <- \(x,mu0,mu1,sig0,sig1,p){
  p*((x-mu1)^2/sig1^4-1/sig1^2)*dnorm(x,mean=mu1,sd=sig1)+(1-p)*((x-mu0)^2/sig0^4-1/sig0^2)*dnorm(x,mean=mu0,sd=sig0)
}


#optimization function
mm.newtraph <- \(muinit,mu0,mu1,sig0,sig1,p,maxit){
  mustar0 <- muinit
  mustar1 <- mustar0-mmfp(mustar0,mu0,mu1,sig0,sig1,p)/mmfp2(mustar0,mu0,mu1,sig0,sig1,p)
  
  maxit <- maxit-1
  if (maxit<=0) {warning("You have reached the maximum number of iterations"); stop()}
  
  while (abs(mustar1-mustar0)>1e-8 | abs(mmfp(mustar1,mu0,mu1,sig0,sig1,p))>1e-8) {
    mustar0 <- mustar1
    mustar1 <- mustar0-mmfp(mustar0,mu0,mu1,sig0,sig1,p)/mmfp2(mustar0,mu0,mu1,sig0,sig1,p)

    # account for newt raph breaking:
    if (mustar1 < min(mu0,mu1) | mustar1 > max(mu0,mu1) && 
        mustar0 < min(mu0,mu1) | mustar0 > max(mu0,mu1)) return(c())
    
    maxit <- maxit-1
    if (maxit<=0) {warning("You have reached the maximum number of iterations"); stop()}
  }
  return(mustar1)
}


#final answer
mixmode <- \(mu0,mu1,sig0,sig1,p,global=TRUE,maxit=100){
  
  if (p==1) return(mu1)
  if (p==0) return(mu0)
  
  results <- c()
  results <- c(results,mm.newtraph(mu0,mu0,mu1,sig0,sig1,p,maxit))
  results <- c(results,mm.newtraph(mu1,mu0,mu1,sig0,sig1,p,maxit))
  
  if (length(results)==1) return(results)
  
  if (global == TRUE){
    
    if (abs(results[1]-results[2])<1e-8){
      #If theyre the same, return either
      return(results[1])
      
    } else if (abs(mmf(results[1],mu0,mu1,sig0,sig1,p)-mmf(results[2],mu0,mu1,sig0,sig1,p))<1e-8) {
      #If theyre different but mmf(each) same, return both
      return(results)
      
    } else {
      #If theyre different but mmf(each) different, return higher mmf
      return(results[1+(mmf(results[1],mu0,mu1,sig0,sig1,p)<mmf(results[2],mu0,mu1,sig0,sig1,p))])
      #If true, coerced to a 1 for results 2. so true when mmf(results[2]) is greater
    }
    
    #global==FALSE below
  } else {
    
    if (abs(results[1]-results[2])<1e-8) {
      #If theyre the same, return either
      return(results[1])
    } else {
      #Otherwise, return results
      return(results)
    }
  }
}

#test graphs
#curve(mmf(x,-1,1,0.25,0.25,.5),-5,5)
#curve(mmfp(x,-1,1,0.25,0.25,.5),-5,5)
#curve(mmfp2(x,-1,1,0.25,0.25,.5),-5,5)