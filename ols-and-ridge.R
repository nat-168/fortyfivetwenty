#ols-and-ridge.R
#Author: River Strumwasser

X <- matrix(c(1,2,3,4,5,6),nrow=3,byrow=T) #found with google; example of XTX-1
y <- c(1,2,3)

ols <- \(y,X){
  #error if y NOT numeric vector, found funcs through ?numeric & ?vector
  if ((is.numeric(y) && is.vector(y)) == FALSE) {warning("y not numeric vector"); stop()}
  #error if X NOT numeric matrix, found func through ?matrix
  if ((is.numeric(X) && is.matrix(X)) == FALSE) {warning("X not numeric matrix"); stop()}
  #error if number of rows X is NOT same length as y
  if (nrow(X) != length(y)) {warning("# Rows of X ≠ Length of y"); stop()}
  
  #return length-p vector with ordinary least squares estimate from 1b
  #(XTX)-1XTy
  A <- t(X)%*%X
  b <- t(X)%*%y
  solve(A,b)
}

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