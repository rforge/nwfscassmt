StepwiseFn <-
function(SearchMat, Data, NDataSets, MinAge, MaxAge, RefAge, MaxSd, MaxExpectedAge, SaveFile, EffSampleSize=0, Intern=TRUE, InformationCriterion="AIC"){

  # Define variables
  Nages = MaxAge+1
  Nreaders = ncol(Data)-1
  ParamVecOpt = SearchMat[,1]
  Stop = FALSE
  IcRecord = NULL
  StateRecord = NULL
  OuterIndex = 0
  
  # Continue searching until Stop==TRUE
  while(Stop==FALSE){
    
    # Increment and intialize variables
    Out
