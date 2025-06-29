 ################
 #### Begin Example 1 (you can do some tests with this one using a = 0, b = 5)
 #### Uncomment to use.
 ################
 # basis <- function(x) do.call(cbind,lapply(1:26,\(i0) dnorm(x,(i0-1)/5,0.2)))
 # attr(basis,"m") <- 26
 ################
 #### End Example 1
 ################

 ################
 #### Begin Example 2 (you can do some tests with this one using a = 0, b = pi)
 #### Uncomment to use.
 ################
 # basis <- function(x) do.call(cbind,lapply(1:20,\(i0) sin(i0*x)))
 # attr(basis,"m") <- 20
 ################
 #### End Example 2
 ################
 
 simpsons <- \(f,a,b,k){
   I <- 0
   h <- (b-a)/(k-1)
   
   for(j in 1:k){
     
     #setting "simpsons coefficient"
     if(j==1 || j==k) {
       sc<-1
     } else if(j%%2 == 1) {
       sc <- 2
     } else {
       sc <- 4
     }

     #adding sum
     I <- I + (h/3) * sc * f(a+(j-1)*h)
   }
   return(I)
 }
 
 midpt <- \(f,a,b,k){
   I <- 0
   h <- (b-a)/k
   for(j in 1:k){
     I <- I + h*f(a+(j-1/2)*h)
   }
   return(I)
 }

 innerprod <- \(basis, g, a, b, k){
   #quadtrature: simpsons if odd k but midpt if even k
   if(k %% 2 == 1) {quad <- simpsons} else {quad <- midpt}
   
   bvec <- c()
   for(i in 1:attr(basis,"m")){
     bvec <- c(bvec,quad(\(x)basis(x)[,i]*g(x),a,b,k))
   }
   return(bvec)
 }
 
 squareint <- \(basis, a, b, k){
   #quadtrature: simpsons if odd k but midpt if even k
   if(k %% 2 == 1) {quad <- simpsons} else {quad <- midpt}
   
   Amat <- c()
   for(i in 1:attr(basis,"m")){
     for(j in 1:attr(basis,"m")){
       Amat <- c(Amat,quad(\(x)basis(x)[,i]*basis(x)[,j],a,b,k))
     }
   }
   attr(Amat,"dim") <- c(attr(basis,"m"),attr(basis,"m"))
   return(Amat)
 }
 
 bestapprox <- \(basis, g, a, b, k){
   solve(squareint(basis,a,b,k),innerprod(basis,g,a,b,k))
 }

 #FOR TESTING
 #g <- function(x) sin(x)*exp(x/10)
 #c <- bestapprox(basis, g, 0, 5, 31)
 #f <- function(x) basis(x)%*%c
 #curve( g(x), 0, 5, n=1e3, lwd=3) ## the original function.
 #curve( f(x), n=1e3, add=TRUE, col=2, lty=2, lwd=3) ## the best approximation 