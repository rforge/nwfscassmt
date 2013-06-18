
TentFn <-
function(X, S=2){
  if(X < 0.5) Y = S*X
  if(X >= 0.5) Y = S*(1-X)
  return(Y)
}
