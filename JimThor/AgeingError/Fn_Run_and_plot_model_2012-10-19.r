
# Stepwise selection for the Punt et al. (2008) ageing error model
# Version 1.2, 2011-09-13

# Maintained by:
# James Thorson
# 971-678-5683
# JimThor@uw.edu

#################### 
# Auxiliary functions
####################

# Short for "column matrix" (i.e. the default in R)
cMx=function(Input){as.matrix(Input)}    

# Short for "row matrix"
rMx=function(Input){
  if(is.vector(Input)){Output<-t(as.matrix(Input))}  
  if(!is.vector(Input)){Output<-as.matrix(Input)}
  Output
}

####################
# Run model
####################
RunFn = function(Data, SigOpt, BiasOpt, NDataSets, MinAge, MaxAge, RefAge, MinusAge, PlusAge, MaxSd, MaxExpectedAge, SaveFile, EffSampleSize=0, Intern=TRUE, AdmbFile=NULL, JustWrite=FALSE){

  # Copy ADMB file 
  if(!is.null(AdmbFile)) file.copy(from=paste(AdmbFile,"agemat.exe",sep=""), to=paste(SaveFile,"agemat.exe",sep=""), overwrite=TRUE)
  
  # Write DAT file
  write(c("# Maximum number of readers",ncol(Data)-1),file=paste(SaveFile,"agemat.dat",sep=""))
    write(c("# Number of data sets",NDataSets),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE)
    write(c("# Number of points per data set",nrow(Data)),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE)
    write(c("# Readers per data set",ncol(Data)-1),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE)
    write("# Which readers per data set",file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE)
    write.table(rMx(1:(ncol(Data)-1)),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE,row.names=FALSE,col.names=FALSE)
    write(c("# Minimum age",MinAge),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE)
    write(c("# Maximum age",MaxAge),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE)
    write(c("# Reference age",RefAge),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE)
    write(c("# Minus groups",MinusAge),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE)
    write(c("# Plus groups",PlusAge),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE)
    write("# Option for bias",file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE)
    write.table(rMx(BiasOpt),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE,row.names=FALSE,col.names=FALSE)
    write("# Option for standard deviation",file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE)
    write.table(rMx(SigOpt),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE,row.names=FALSE,col.names=FALSE)
    write(c("# Option for effective sample size",EffSampleSize),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE)
    write(c("# Use Par File (1=Yes)",0),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE) 
  # Write initial values  
    # Bias 
    write("\n# Min, Max, Init, Phase for Bias",file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE)
    for(BiasI in 1:length(BiasOpt)){
      # No bias
      if(BiasOpt[BiasI]<=0){}
      # Linear bias
      if(BiasOpt[BiasI]==1) write.table(rMx(c(0.001,3,1,2)),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE,col.names=FALSE,row.names=FALSE)      
      # Curvilinear bias = 0.5+Par1 + (Par3-Par1)/(1.0-mfexp(-Par2*(float(MaxAge)-1)))*(1.0-mfexp(-Par2*(float(Age1)-1)))
      # Starting value must be non-zero
      if(BiasOpt[BiasI]==2){
        write.table(rMx(c(0.001,10,1,2)),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE,col.names=FALSE,row.names=FALSE)      
        write.table(rMx(c(-10,1,0.1,2)),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE,col.names=FALSE,row.names=FALSE)      
        write.table(rMx(c(0.001,MaxAge*2,MaxAge,2)),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE,col.names=FALSE,row.names=FALSE)      
      }
    }
    # Sigma
    write("\n# Min, Max, Init, Phase for Sigma",file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE)
    for(SigI in 1:length(SigOpt)){
      # No error
      if(SigOpt[SigI]<=0){}
      # Linear CV
      if(SigOpt[SigI]==1) write.table(rMx(c(0.001,3,0.1,2)),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE,col.names=FALSE,row.names=FALSE)      
      # Curvilinear SD = Par1 + (Par3-Par1)/(1.0-exp(-Par2*(100-1)))*(1.0-exp(-Par2*(1:100)))
      # Starting value must be non-zero
      if(SigOpt[SigI]==2){
        write.table(rMx(c(0.001,100,1,2)),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE,col.names=FALSE,row.names=FALSE)      
        write.table(rMx(c(-10,1,0.1,2)),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE,col.names=FALSE,row.names=FALSE)      
        write.table(rMx(c(0.001,100,10,2)),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE,col.names=FALSE,row.names=FALSE)      
      }
      # Curvilinear CV
      # Curvilinear CV = Par1 + (Par3-Par1)/(1.0-exp(-Par2*(100-1)))*(1.0-exp(-Par2*(1:100)))
      # Starting value must be non-zero
      if(SigOpt[SigI]==3){
        write.table(rMx(c(0.001,3,0.1,2)),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE,col.names=FALSE,row.names=FALSE)      
        write.table(rMx(c(-10,1,0.1,2)),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE,col.names=FALSE,row.names=FALSE)      
        write.table(rMx(c(0.001,3,0.1,2)),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE,col.names=FALSE,row.names=FALSE)      
      }
    }
    # Probs (i.e. age-composition probability relative to reference age)
    write("\n# Min, Max, Phase for Probs",file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE)
    write.table(rMx(c(-20,20,2)),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE,col.names=FALSE,row.names=FALSE)
    # Slopes
    write("\n# Min, Max, Init, Phase for slopes",file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE)
    for(DataSetI in 1:NDataSets){
      if(MinAge>1) write.table(rMx(c(-10,0,0,1)),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE,col.names=FALSE,row.names=FALSE)      
      write.table(rMx(c(-10,0,0,1)),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE,col.names=FALSE,row.names=FALSE)      
    }    
  # Write dataset    
    write("\n# Data set",file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE)
    write.table(Data,file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE,row.names=FALSE, col.names=FALSE)
    write(c("# Test number",123456),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE)

  # Run ADMB file
  if(JustWrite==FALSE){
    setwd(SaveFile)
    Output = shell("agemat.exe -est",intern=Intern)
    Admb = scan(paste(SaveFile,"agemat.par",sep=""),comment.char="#",quiet=TRUE)
  }
}

######################
# Plot results
######################
PlotOutputFn = function(Data, MaxAge, SaveFile, PlotType="PDF"){

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
    # This is done by assigning each otolith to the age which has maximum posterior probability (i.e. the conditional mode, as is typically done for random effects)
  AgeProbs = array(NA, dim=c(nrow(Data),Ages), dimnames=list(paste("Otolith",1:nrow(Data)),paste("TrueAge",0:MaxAge)))
    OtI = AgeI = ReadI = 1
  for(OtI in 1:nrow(Data)){
  for(AgeI in 1:Ages){
    AgeProbs[OtI,AgeI] = AgeStruct[AgeI,2]
    for(ReadI in 1:Nreaders){
      if(Data[OtI,ReadI+1]!=-999){
        AgeRead = Data[OtI,ReadI+1]
        AgeProbs[OtI,AgeI] = AgeStruct[AgeI,2] * (MisclassArray[ReadI,AgeI,AgeRead+1])^Data[OtI,1]
      }
    }
  }}
  # Remove MaxAge before calculating "TrueAge" because the MaxAge is a plus-group, and ends up with maximum probability for most ages in the upper tail
  TrueAge = apply(AgeProbs, MARGIN=1, FUN=function(Vec){order(Vec[-length(Vec)],decreasing=TRUE)[1]})
  
  # Plot estimated age structure
  Temp = ifelse(Data[,-1]==-999,NA,Data[,-1])
    Temp = tapply(ifelse(is.na(Temp),0,1), INDEX=Temp, FUN=sum)
    Prop = rep(0,MaxAge+1)
    Prop[as.numeric(names(Temp))+1] = Temp
    cbind(0:MaxAge, Prop, round(AgeStruct[,2]*sum(Prop),1))
    cbind(0:MaxAge, round(Prop/sum(Prop),3), round(AgeStruct[,2],3))
  Plot = function(){
    par(mar=c(3,3,2,0),mgp=c(1.5,0.25,0),tck=-0.02,oma=c(0,0,0,0)+0.1)
    plot(x=AgeStruct[,1],y=AgeStruct[,2],type="s",lwd=2,xlab="Age",ylab="Prop",main="Estimated=Black, Observed=Red")
    hist(ifelse(Data[,-1]==-999,NA,Data[,-1]),add=TRUE,freq=FALSE,breaks=seq(0,MaxAge,by=1),col=rgb(red=1,green=0,blue=0,alpha=0.30))
  }
  if(PlotType=="PDF"){
    pdf(paste(SaveFile,"Estimated vs Observed Age Structure.pdf",sep=""),width=6,height=6)
      Plot()
    dev.off()
  }
  if(PlotType=="JPG"){
    jpeg(paste(SaveFile,"Estimated vs Observed Age Structure.jpg",sep=""),width=6,height=6,units="in",res=200)
      Plot()
    dev.off()
  }
  
  # Plot true age against different age reads
  Ncol=ceiling(sqrt(Nreaders))
    Nrow=ceiling(Nreaders/Ncol)
  Plot = function(){    
    par(mfrow=c(Nrow,Ncol),mar=c(3,3,2,0),mgp=c(1.5,0.25,0),tck=-0.02,oma=c(0,0,5,0)+0.1)
    for(ReadI in 1:Nreaders){
      Temp = cbind(TrueAge, Data[,ReadI+1]+0.5)   # Add 0.5 to match convention in Punt model that otoliths are read half way through year
      Temp = Temp[which(Data[,ReadI+1]!=-999),]   # Exclude rows with no read for this reader
      plot(x=Temp[,1],y=Temp[,2],ylim=c(0,MaxAge),xlim=c(0,MaxAge),col=rgb(red=0,green=0,blue=0,alpha=0.2),xlab="Mode predicted age | parameters",ylab="Read age",lwd=2,main=paste("Reader",ReadI),pch=21,cex=0.2)
      lines(x=c(0,MaxAge),y=c(0,MaxAge), lwd=1,lty="dashed")
      lines(x=ErrorAndBiasArray['True_Age',,ReadI],y=ErrorAndBiasArray['Expected_age',,ReadI],type="l",col="red",lwd=1)
      lines(x=ErrorAndBiasArray['True_Age',,ReadI],y=ErrorAndBiasArray['SD',,ReadI],type="l",col="blue",lwd=1)
      lines(x=ErrorAndBiasArray['True_Age',,ReadI],y=ErrorAndBiasArray['Expected_age',,ReadI] + 2*ErrorAndBiasArray['SD',,ReadI],type="l",col="red",lwd=1,lty="dashed")
      lines(x=ErrorAndBiasArray['True_Age',,ReadI],y=ErrorAndBiasArray['Expected_age',,ReadI] - 2*ErrorAndBiasArray['SD',,ReadI],type="l",col="red",lwd=1,lty="dashed")
    }
    mtext(side=3,outer=TRUE, text="Reads(dot), Sd(blue), expected_read(red solid line),\n and 95% CI for expected_read(red dotted line)",line=1)
  }
  if(PlotType=="PDF"){
    pdf(paste(SaveFile,"True vs Reads (by reader).pdf",sep=""),width=Ncol*3,height=Nrow*3)
      Plot()
    dev.off()
  }
  if(PlotType=="JPG"){
    jpeg(paste(SaveFile,"True vs Reads (by reader).jpg",sep=""),width=Ncol*3,height=Nrow*3,units="in",res=200)
      Plot()
    dev.off()
  }
  
  ## AIC
  Nll = as.numeric(scan(paste(SaveFile,"agemat.par",sep=""),comment.char="%", what="character", quiet=TRUE)[11])
    Df = as.numeric(scan(paste(SaveFile,"agemat.par",sep=""),comment.char="%", what="character", quiet=TRUE)[6])
    n = sum(ifelse(Data[,-1]==-999,0,1))
    Aic = 2*Nll + 2*Df
    Aicc = Aic + 2*Df*(Df+1)/(n-Df-1) 
    Bic = 2*Nll + Df*log(n)
      
  Output = list(Aic=Aic, Aicc=Aicc, Bic=Bic)
  return(Output)
}

