#Author: River Strumwasser
#mle.R

#l(theta)
#Borrowed from edstem
bacterialoglik <- function(theta, v, n, x) {
  stopifnot(length(v) == length(n) && length(n) == length(x) )
  sum(lchoose(n,x)) +
    sapply(theta, \(thetaj)
           - thetaj*sum(v*x) + sum( (n - x)*log(1 - exp(-thetaj*v)) ))
}

#l'(theta)
bacterialoglikp <- function(theta, v, n, x) {
  stopifnot(length(v) == length(n) && length(n) == length(x) )
  sapply(theta, \(thetaj)
         - sum(v*x) + sum( (n - x)*v*exp(-thetaj*v)/(1 - exp(-thetaj*v)) ))
}

#l''(theta), using quotient rule
bacterialoglikp2 <- function(theta, v, n, x) {
  stopifnot(length(v) == length(n) && length(n) == length(x) )
  sapply(theta, \(thetaj) {
    -sum((n - x)*v^2*
          exp(-thetaj*v) / (1 - exp(-thetaj*v))^2)
  })
}


#Maximization Function

#initial guess derivation:
# ∑xi = ∑ni*e^(-theta*vi)
# ∑xi/∑ni = ∑ni*e^(-theta*vi)/∑ni --> weighted avg, by jensens & since ni in ∑
# ∑xi/∑ni ~ e^(-theta*vbar), vbar is ∑ni*vi/∑ni
# so initial guess by resub: theta = -∑ni/∑ni*vi *log(∑xi/∑ni)
# in r: -sum(n)/sum(n*v)*log(sum(x)/sum(n))

mle<- \(v,n,x,maxit=100) {
  #special case 1 - x all 0s
  if (all(x==0)) return(Inf) #shown on graph
  
  #special case 2 - xi = ni for all i --> graph shows max at 0
  if (all(x==n)) return(0) #all() found from ?"==" to ?all.equal to ?all
  
  #general soln
  theta0 <- -sum(n)/sum(n*v)*log(sum(x)/sum(n))
  theta1 <- theta0-bacterialoglikp(theta0,v,n,x)/bacterialoglikp2(theta0,v,n,x)
  
  maxit <- maxit-1
  if (maxit<=0) {warning("You have reached the maximum number of iterations"); stop()}
  
  while (abs(theta1-theta0)>1e-8 |abs(bacterialoglikp(theta1,v,n,x))>1e-8) {
    theta0 <- theta1
    theta1 <- theta0-bacterialoglikp(theta0,v,n,x)/bacterialoglikp2(theta0,v,n,x)
    maxit <- maxit-1
    if (maxit<=0) {warning("You have reached the maximum number of iterations"); stop()}
  }
  return(theta1)
}

#test numbers
v.ex<-c(1/100, 1/10, 1)
n.ex<-c(5,5,2)
x.ex<-c(4,1,0)

#test graphs
#curve(bacterialoglik(x,v.ex,n.ex,x.ex),5,100)
#curve(bacterialoglikp(x,v.ex,n.ex,x.ex),5,100)
#curve(bacterialoglikp2(x,v.ex,n.ex,x.ex),5,100)