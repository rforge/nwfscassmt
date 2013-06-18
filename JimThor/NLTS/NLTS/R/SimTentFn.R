
SimTentFn <-
function(Nobs, S){
  Y = numeric(Nobs)
  Y[1] = runif(1)
  for(i in 2:Nobs) Y[i] = TentFn(Y[i-1], S=S)
  return(Y)
}
