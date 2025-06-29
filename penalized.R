#penalized.R
#Author: River Strumwasser

library(splines)
knots <- c(60.084873374401099, 60.084873374401099,
          60.084873374401099, 60.084873374401099, 60.175460524357945,
          60.266047674314798, 60.356634824271644, 60.447221974228498,
          60.537809124185344, 60.628396274142197, 60.718983424099044,
          60.80957057405589, 60.900157724012743, 60.990744873969589,
          61.081332023926443, 61.171919173883289, 61.262506323840142,
          61.353093473796989, 61.443680623753842, 61.534267773710688,
          61.624854923667534, 61.715442073624388, 61.806029223581234,
          61.896616373538087, 61.987203523494934, 62.077790673451787,
          62.168377823408633, 62.25896497336548, 62.349552123322333,
          62.440139273279179, 62.530726423236032, 62.621313573192879,
          62.711900723149732, 62.802487873106578, 62.893075023063432,
          62.983662173020278, 63.074249322977124, 63.164836472933978,
          63.255423622890824, 63.346010772847677, 63.436597922804523,
          63.527185072761377, 63.617772222718223, 63.708359372675069,
          63.798946522631923, 63.889533672588769, 63.980120822545622,
          64.070707972502476, 64.161295122459322, 64.251882272416168,
          64.342469422373014, 64.433056572329861, 64.523643722286721,
          64.614230872243567, 64.704818022200413, 64.79540517215726,
          64.88599232211412, 64.976579472070966, 65.067166622027813,
          65.157753771984659, 65.248340921941505, 65.338928071898366,
          65.429515221855212, 65.520102371812058, 65.610689521768904,
          65.701276671725765, 65.791863821682611, 65.882450971639457,
          65.973038121596304, 66.06362527155315, 66.15421242151001,
          66.244799571466856, 66.335386721423703, 66.335386721423703,
          66.335386721423703, 66.335386721423703)
a <- min(knots); b <- max(knots)
basis <- \(x,derivs=0) splineDesign(knots,x,derivs=derivs)
k <- ncol(basis(a))
co2 <- read.csv("co2.csv")

#2a
ghat <- \(x,c){
  as.vector(basis(x)%*%c)
}

#2b
xi <- co2$years
yi <- co2$co2

chat_sse <- solve(t(basis(xi))%*%basis(xi),t(basis(xi))%*%yi)
plot(xi, yi)
lines(xi, ghat(xi, chat_sse))

#2c

P <- matrix(rep(0,ncol(basis(xi))^2), nrow=ncol(basis(xi)), ncol=ncol(basis(xi)))

for (j in 1:ncol(basis(xi))) {
  for (i in j:ncol(basis(xi))) { #start i means can do lower tri matrix & j≤i
    if (i-j < 4) { #cancels out under 4, no abs bc j≤i
      phiint <- \(x)basis(x, derivs = 2)[,i] * basis(x, derivs = 2)[,j]
      P[i,j] <- integrate(phiint, lower = knots[i], upper = knots[j+4])$value
      P[j,i] <- P[i,j] #lower tri symm
    }
  }
}

lambda <- 0.0005
chat_pensse <- solve(t(basis(xi))%*%basis(xi)+lambda*P,t(basis(xi))%*%yi)
plot(xi, yi)
lines(xi, ghat(xi, chat_pensse))


#3a

#c = (btb-lp)-1bty
#y = bc = b(btb-lp)-1bty
#sl = b(btb+lp)-1bt

b <- basis(xi)

S <- \(lambda){
  b%*%solve(t(b)%*%b+lambda*P,t(b))
}

#3b

ocv <- Vectorize(\(lambda){
  Sl <- S(lambda)
  yhat <- Sl%*%yi
  return(mean(((yi-yhat)/(1-diag(Sl)))^2)) #diag for ]ii
}, "lambda")

#3c

lambda3c <- optimize(ocv, interval=c(1e-6, 20))$minimum #guessed intvl & ran to check
lambda3c

#7.276763 is the ocv lambda here, which is far greater than suggested
#looks far smoother!! but doesn't align with obvious oscillation pattern