
SimPredPreyFn <-
function(Nobs, Nt=100, A=0.04, B=0.05, C=0.2, E=1){
  X = Y = matrix(NA, ncol=Nobs, nrow=Nt)
  X[1,1] = runif(1)
  Y[1,1] = runif(1)
  for(T in 1:Nobs){
    for(dT in 2:Nt){
      X[dT,T] = X[dT-1,T] + ( A*X[dT-1,T] - B*X[dT-1,T]*Y[dT-1,T] ) / Nt
      Y[dT,T] = Y[dT-1,T] + ( E*B*X[dT-1,T]*Y[dT-1,T] - C*Y[dT-1,T] ) / Nt
    }
    if(T<Nobs){
      X[1,T+1] = X[Nt,T]
      Y[1,T+1] = Y[Nt,T]
    }
  }
  List = list(X=X[1,], Y=Y[1,])
  matplot(cbind(X[1,],Y[1,]), type="l", xlab="Year", ylab="Abundance")
  return(List)
}
