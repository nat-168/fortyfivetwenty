## Test your sampler with these values for the parameters of the prior gamma density
alpha <- 0.7
beta <- 0.01
## Test your sampler with these data (these are the same ones from HW4)
v <- c(1/100, 1/10, 1)
n <- c(5, 5, 2)
x <- c(4, 1, 0)
## Number of MCMC samples (don't confuse this with the number of sterile tubes `n` above)
Nmcmc <- 1e4


### Put your code below this line

#from edstem:

bacterialoglik <- function(theta, v, n, x) {
  stopifnot(length(v) == length(n) && length(n) == length(x) )
  sum(lchoose(n,x)) +
    sapply(theta, \(thetaj)
           - thetaj*sum(v*x) + sum( (n - x)*log(1 - exp(-thetaj*v)) ))
}

alpha0 <- 1

f <- \(theta,v,n,x) exp(bacterialoglik(theta,v,n,x))*dgamma(theta,alpha,beta)
p <- \(x,y) dgamma(x,alpha,alpha/y)

alphamcmc <- \(xt,y,v,n,x) (f(y,v,n,x)*p(xt,y))/(f(xt,v,n,x)*p(y,xt))
X <- rep(0,Nmcmc+1)
X[1] <- 20

for (t in 1:Nmcmc){
  Y <- rgamma(1,alpha,alpha/X[t])
  U <- runif(1)
  #print(f(Y)*p(X[t],Y))
  #print(f(X[t])*p(Y,X[t]))
  #print(alphamcmc(X[t],Y))
  if(U<alphamcmc(X[t],Y,v,n,x)) {
    X[t+1] <- Y #accept proposal
  } else {
    X[t+1] <- X[t] #stay where were
  }
}

#hist(X,freq=FALSE,breaks=50)
#curve(dgamma(x,2,1),add=TRUE)