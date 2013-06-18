
NltsFn <-
function(Y, PredInterval, Nembed, Theta=NA, Method){

  # Table of leave-on-out predictions
  Table = array(NA, dim=c(length(Y),2))
  for(i in (Nembed+PredInterval):length(Y)){
    Pred = NltsPred(Y=Y, PredInterval=PredInterval, Nembed=Nembed, PredNum=i, Theta=Theta, Method=Method)
    Table[i,1:2] = c(Y[i],Pred)
  }

  # Calculate and return Correlation
  Corr = cor(Table[,1], Table[,2], use="complete.obs")
  Return = list(Table=Table, Corr=Corr)

  return(Return)
}
