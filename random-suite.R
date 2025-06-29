#Author: River Strumwasser
#random-suite.R

#PROBLEM 1: Copied Over Inverse CDF
rinvcdf <-\(n,Finv){
  Y <- runif(n,0,1)
  X <- sapply(Y,Finv)
  return(X)
}

#PROBLEM 2: Copied Over Rejection (both)
simplereject <- \(n, f, a, b, k, check.k){
  U <- runif(n,a,b)
  V <- runif(n,0,k)
  Unew <- U[f(U)>V]
  if(check.k) if(mean(f(Unew)>k)>0) {warning("k is less than f height") ; stop()}
  return(Unew)
}

generalreject <- \(n, f, k, gpdf, gsample, check.k){
  X <- gsample(n)
  Y <- runif(n,0,k*gpdf(X))
  Xnew <- X[f(X)>Y]
  if(check.k) if(mean(f(Xnew)>k*gpdf(Xnew))>0) {warning("k*g is less than f height") ; stop()}
  return(Xnew)
}

#3a. (1 point) randexp(n, lambda): n samples from Exponential(λ)
#y=1-e^lx, 1-y=-e^lx, log(1-y)=-lx, x=-log(1-y)/l, x=-log(y)/l
fa <- \(x,lambda) {
  -log(x)/lambda
}
randexp<- \(n, lambda){
  rinvcdf(n,\(x)fa(x,lambda))
}

#3b. (1 point) randnorm(n, mu, sigma): n samples from N(µ,σ2)

#rlaplace <- \(n){
#  ifelse(binom(n,1,0.5), exp(n,1), -exp(n,1))
#}
#this basic format (from class/rejection.R) uses stats package funcs
#replace rexp with randexp from 3a, & create rbinom replacement from scratch:

randbernoulli <- \(n,p){
  X <- runif(n,0,1)
  Y <- ifelse(X<p,1,0)
  return(Y)
}

glaplace <- \(x) 0.5*exp(-abs(x))

randlaplace <- \(n){
  ifelse(randbernoulli(n,0.5), randexp(n,1), -randexp(n,1))
}

funcnorm <- \(x,mu,sigma){
  1/(sqrt(2*pi)*sigma)*exp(-(x-mu)^2/(2*sigma^2))
}

randnorm <- \(n,mu,sigma){
  gennums <- c()
  repeat{ #loop here because of rejection to get n samples instead of ~n/k
    n01 <- generalreject(n,\(x)funcnorm(x,0,1),2,glaplace,randlaplace,TRUE)
    nms <- mu+sigma*n01
    gennums <- c(gennums,nms)
    
    if(length(gennums) >= n) return(gennums[1:n])
  }
}

#3c. (1 point) randchisq(n, df): n samples from χ2 df
randchisq <- \(n,df){
  X <- randnorm(n*df,0,1)^2
  Xmat <- matrix(X,nrow=n,ncol=df)
  Xmatsummed <- apply(Xmat,1,sum)
  return(Xmatsummed)
}

#hist(randchisq(1000,50),breaks=100,freq=FALSE)
#curve(dchisq(x,50),add=TRUE)

#3d. (1 point) randgamma(n, alpha, beta): n samples from a Gamma(α,β)

# chisq df k makes gamma alpha df/2, so chisq df = 2*alpha
# chisq gen beta 2, so div 2 and * beta for scale

randgamma <- \(n, alpha, beta){
  X <- randchisq(n,2*alpha)
  Y <- X/(2*beta)
  return(Y)
}

#hist(randgamma(10000,5,40),breaks=50,freq=FALSE)
#curve(dgamma(x,5,scale=40),add=TRUE)

#3e. (1 point) randt(n, df): n samples from Student’s t-distribution, df

#T = Z/sqrt(U/df), indep.

randt <- \(n, df){
  Z <- randnorm(n,0,1)
  U <- randchisq(n,df)
  X <- Z/sqrt(U/df)
  return(X)
}

#hist(randt(10000,1),breaks=60000,freq=FALSE,xlim=c(-5,5))
#curve(dt(x,1),add=TRUE)

#3f. randgammaOptional(n, alpha, beta): nsamples from a Gamma(α,β) α>0
randgammaOptional <- \(n, alpha, beta){
  if(alpha<=1) {
    return(runif(n,0,1)*randgammaOptional(n,alpha+1,beta))
  } else {
    f <- \(x)(beta**alpha)*(x**(alpha-1))*exp((-beta)*x)/gamma(alpha)
    #gpdf <- \(x)(d_gamma(x,floor(alpha),beta)+d_gamma(x,floor(alpha)+1,beta))/2
    #mixed gamma gpdf, scale *2 do encompass any alpha between floor and ceiling
    k <- 2
    gsample <- \(x)c(randgamma(x/2,floor(alpha),beta),randgamma(x/2,floor(alpha)+1,beta))
    
    gpdf <- \(x)(
      ((beta**floor(alpha))*(x**(floor(alpha)-1))*exp((-beta)*x)/gamma(floor(alpha)))+
      ((beta**floor(alpha)+1)*(x**(floor(alpha)))*exp((-beta)*x)/gamma(floor(alpha)+1))
      )/2
    
    gennums<-c()
    repeat{ #loop here because of rejection to get n samples instead of ~n/k
      gamsamp <- generalreject(n, f, k, gpdf, gsample, FALSE)
      #print(length(gamsamp))
      gennums <- c(gennums,gamsamp)
      
      if(length(gennums) >= n) return(gennums[1:n])
    }
  }
}

#3g. randbetaOptional(n, alpha, beta): nsamples from a Beta(α,β) distribution.
#𝑋∼Gamma(𝛼,1) and 𝑌∼Gamma(𝛽,1)
# X/(X+Y)
randbetaOptional <- \(n, alpha, beta){
  X <- randgammaOptional(n,alpha,1)
  Y <- randgammaOptional(n,beta,1)
  return(X/(X+Y))
}

#Test before done w/ gammaOptional:
#randbetaOptional <- \(n, alpha, beta){
#  X <- randgamma(n,alpha,1)
#  Y <- randgamma(n,beta,1)
#  return(X/(X+Y))
#}
