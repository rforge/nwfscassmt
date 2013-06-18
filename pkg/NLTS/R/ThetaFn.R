
ThetaFn <-
function(Y, PredInterval, ThetaSet=c(0,exp(seq(-3,3))), Nembed){
  Corr = numeric(length(ThetaSet))
  for(ThetaI in 1:length(ThetaSet)){
    Output = NltsFn(Y, PredInterval=PredInterval, Nembed=Nembed, Method="Smap", Theta=ThetaSet[ThetaI])
    Corr[ThetaI] = Output$Corr
  }
  plot(y=Corr, x=ThetaSet, type="b", ylim=c(0,1), xlim=range(ThetaSet))
  Max = which.max(Corr)
  return(list(Corr=rbind(ThetaSet,Corr), Max=Max))
}
