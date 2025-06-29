#lasso.R
#Author: River Strumwasser

ftwiddle <- \(betaj, xjrj, xj2, lambda){
  #value of 3a (2) with arguments:
  #betaj = βj
  #xjrj = ∑(i=1 to n) xijrij
  #xj2 = ∑(i=1 to n) xij2
  #lambda = λ.
  
  #f~λ,j(βj) = -2βj*∑(i=1 to n) xijrij+βj^2*∑(i=1 to n) xij^2+λ*|βj|
  -2 * betaj * xjrj + betaj^2 * xj2 + lambda*abs(betaj)
}


X.ex <- matrix(c(1,2,3,4,5,6),nrow=3,byrow=T) #found with google; example of XTX-1
y.ex <- c(1,2,3)
lambda.ex <- 3

#copied from page
golden <- \(l,r,f,eps) {
  phi <- (1 + sqrt(5))/2
  m <- l + (r - l)/(1+phi)
  while( abs(r - l) > eps ) {
    if( r - m > m - l ) {
      n <- m + (r - m)/(1+phi)
    } else {
      n <- m - (m - l)/(1+phi)
    }

    if( f(m) > f(n) ) {
      if( m < n ) {
        r <- n
      } else {
        l <- n
      }
    } else {
      if( m < n ) {
        l <- m
      } else {
        r <- m
      }
      m <- n
    }
  }
  return(m)
}

#ridge from ols-and-ridge for initial beta guesses
ridge <- \(y,X,lambda){
  #conditions copied from ols:
  #error if y NOT numeric vector, found funcs through ?numeric & ?vector
  if ((is.numeric(y) && is.vector(y)) == FALSE) {warning("y not numeric vector"); stop()}
  #error if X NOT numeric matrix, found func through ?matrix
  if ((is.numeric(X) && is.matrix(X)) == FALSE) {warning("X not numeric matrix"); stop()}
  #error if number of rows X is NOT same length as y
  if (nrow(X) != length(y)) {warning("# Rows of X ≠ Length of y"); stop()}
  
  #error if lambda is NOT a single positive number
  if ((length(lambda) == 1 && lambda > 0 && is.numeric(lambda)) == FALSE) {
    warning("Lambda not single positive number"); stop()
  }
  
  #return length-p vector with ridge estimate from 1c
  #(XTX+λI)-1XTy
  
  #diag() from ??identity then ??diagonal
  A <- t(X)%*%X + lambda * diag(ncol(X)) #X nxp
  b <- t(X)%*%y
  
  ridgeest <- solve(A,b)
  attr(ridgeest, "lambda") <- lambda
  return(ridgeest)
}

lasso <- \(y, X, lambda) {
  #conditions copied from ridge:
  #error if y NOT numeric vector, found funcs through ?numeric & ?vector
  if ((is.numeric(y) && is.vector(y)) == FALSE) {warning("y not numeric vector"); stop()}
  #error if X NOT numeric matrix, found func through ?matrix
  if ((is.numeric(X) && is.matrix(X)) == FALSE) {warning("X not numeric matrix"); stop()}
  #error if number of rows X is NOT same length as y
  if (nrow(X) != length(y)) {warning("# Rows of X ≠ Length of y"); stop()}
  #error if lambda is NOT a single positive number
  if ((length(lambda) == 1 && lambda > 0 && is.numeric(lambda)) == FALSE) {
    warning("Lambda not single positive number"); stop()
  }
  
  #return length-p vector with lasso estimate from 3a
  #coordinate descent through each bj
  
  #for self check iteration
  #iter <- 0
  
  #special case: if all 0s lasso 0
  if(all(X==0)) {
    lassoest <- rep(0,ncol(X))
    attr(lassoest, "dim") <- c(ncol(X),1)
    attr(lassoest, "lambda") <- lambda
    return(lassoest)
  }
  
  #start with initial guesses
  beta <- ridge(y,X,lambda)
  betaold <- beta
  j = 0
  repeat{
    
    #for self check iteration
    #iter <- iter+1
    
    #start with counter
    if (j == length(beta)) j <- 0
    j <- j + 1
    
    #our golden section maximizes: so min -ftwiddle 
    #betaj = βj *TO BE MAXED
    
    #xjrj = ∑(i=1 to n) xijrij
    rij <- y - X %*% beta + X[,j] * beta[j] #all cols but add back j to cancel
    xjrj <- sum(X[,j]*rij)
    
    #xj2 = ∑(i=1 to n) xij2
    xj2 <- sum(X[,j]^2) #sum tested on X first to be discovered after initial for loop
    
    #lambda = λ *ALREADY ASSIGNED
    
    betaold[j] <- beta[j]
    
    #l&r bounds:
    # 0 ≤ |βlasso| ≤ |βols|
    # 0 ≤ |βlasso| ≤ |cβridge|
    # -|cβridge| ≤ βlasso ≤ |cβridge|
    # c should depend on p
    #but arbitrary 100 is fine for most cases & adds few runtime
    
    l <- -100*abs(ridge(y,X,lambda)[j])
    r <- 100*abs(ridge(y,X,lambda)[j])
    
    #curve(-ftwiddle(x,xjrj,xj2,lambda),l,r)
    beta[j] <- golden(l,r,\(bj)-ftwiddle(bj,xjrj,xj2,lambda),1e-8)
    
    #condition
    if (max(abs(betaold-beta))<1e-8) break
    #{print(iter); break} for check iter
  }
  
  #ending formula copied from ridge:
  lassoest <- beta
  attr(lassoest, "lambda") <- lambda
  return(lassoest)
}