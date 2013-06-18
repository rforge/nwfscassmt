
NltsPred <-
function(Y, PredInterval, Nembed, PredNum, Method, Theta=NA){

  if(PredNum<=Nembed) stop("Can't predict less than or equal to Nembed")
  if( (PredInterval%%1)!=0 || PredInterval<=0 ) stop("PredInterval must be a positive integer")

  L = length(Y)
  Neighbors = Nembed+1

  # Basis expansion by lag
  Mat = matrix(NA, nrow=L, ncol=2+Nembed)
  colnames(Mat) = c("ObsNum","Y",paste("X_lag",(1:Nembed-1+PredInterval),sep=""))
  Mat[,'ObsNum'] = 1:L
  Mat[,'Y'] = Y
  for(i in 1:Nembed){
    Which = Mat[,'ObsNum']-PredInterval+1-i
    Which = ifelse(Which<1,NA,Which)
    Mat[,2+i] = Y[Which]
  }

  # Make library
  Which = which(Mat[,'ObsNum']==PredNum)
  PredObs = Mat[Which,]
  LibObs = na.omit(Mat[-Which,])

  # Identify neighbors
  Distance = sqrt( rowSums( ( LibObs[,-c(1:2),drop=FALSE] - outer(rep(1,nrow(LibObs)),PredObs[-c(1:2)]) )^2 ) )

  # Simplex
  if(Method=="Simplex"){
    WhichNeighbor = order(Distance)[1:Neighbors]
    Ypred = LibObs[WhichNeighbor,'Y']
    Pred = mean(Ypred)
  }
  
  # Smap (derived from Glaser et al. 2011 CJFAS)
  if(Method=="Smap"){
    Weight = exp(-Theta*Distance/mean(Distance))
    B = Weight * LibObs[,'Y']
    A = Weight * cbind(rep(1,nrow(LibObs)), LibObs[,-match(c('ObsNum','Y'),colnames(LibObs))])
    SVD = svd(A)   # A = SVD$u %*% diag(SVD$d) %*% t(SVD$v); solve(A) = SVD$v %*% diag(1/SVD$d) %*% t(SVD$u)
    Ainv = SVD$v %*% diag(1/SVD$d) %*% t(SVD$u)
    C = Ainv %*% B
    Pred = c(1, PredObs[-match(c('ObsNum','Y'),colnames(LibObs))]) %*% C
  }

  return(Pred)
}
