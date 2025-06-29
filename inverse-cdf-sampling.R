#Author: River Strumwasser
#inverse-cdf-sampling.R


#1a. n random samples with cdf F(x)
rinvcdf <-\(n,Finv){
  Y <- runif(n,0,1)
  X <- sapply(Y,Finv)
  return(X)
}

#1b. cauchy F-1(x)

#F(x)=1/πarctan(x)+1/2
#y=1/πarctan(x)+1/2
#(y-1/2)π=arctan(x)
#tan(π(y-1/2))=x

cauchyFinv <- \(x){
  tan(pi*(x-1/2))
}

#1c. cauchyFinv(x) to test your rinvcdf(x,Finv)

cauchy <- \(x) atan(x)/pi+1/2

T <- rinvcdf(1000,cauchyFinv)
mean(T<0.3)
cauchy(0.3)

mean(T<4)
cauchy(4)


mean(T<20)
cauchy(20)

