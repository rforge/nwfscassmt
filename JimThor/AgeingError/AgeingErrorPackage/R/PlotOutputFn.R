PlotOutputFn <-
function(Data, MaxAge, SaveFile, PlotType="PDF"){

  # Interpret inputs
  Nreaders = ncol(Data)-1
  Ages = Nages = MaxAge+1
  
  # Read REP file
  Rep = scan(paste(SaveFile,"agemat.rep",sep=""),comment.char="%", what="character", quiet=TRUE)
    
  # Read Misclassification rates
  Grep = grep("reader#", Rep)
    MisclassArray = array(NA, dim=c(Nreaders,Ages,Ages), dimnames=list(paste("Reader",1:Nreaders),paste("TrueAge",0:MaxAge),paste("EstAge",0:MaxAge)))   # Reader, TrueAge, EstAge
    for(i in 1:Nreaders){
      MisclassArray[i,,] = matrix(as.numeric(Rep[Grep[i]+1+1:(Ages^2)]),ncol=Ages,byrow=TRUE)
    }
  
  # Input estimated age-structure
  Grep = grep("age-structure", Rep)
    AgeStruct = matrix(as.numeric(Rep[Grep[1]+7+1:(2*Ages)]),ncol=2,byrow=TRUE)
    
  # Input reader error and bias  
  Grep = grep("age",Rep)[3]
    Temp = Rep[Grep+1:(5*Nages*Nreaders)]
    ErrorAndBiasArray = array(as.numeric(Temp), dim=c(5,Nages,Nreaders), dimnames=list(c("Reader","True_Age","CV","SD","Expected_age"),paste("Age",0:MaxAge),paste("Reader",1:Nreaders)))
  
  # Estimate unobserved age for each otolith
    # This is done
