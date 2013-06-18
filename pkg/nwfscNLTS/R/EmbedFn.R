
EmbedFn <-
function(Y, PredInterval, Candidates=1:5){
  Corr = numeric(length(Candidates))
  for(E in Candidates){
    Output = NltsFn(Y, PredInterval=PredInterval, Nembed=Candidates[E], Method="Simplex", Theta=NA)
    Corr[E] = Output$Corr
  }
  plot(Corr, type="b", ylim=c(0,1), xlim=range(Candidates))
  Max = which.max(Corr)
  return(Max)
}
