RunFn <-
function(Data, SigOpt, BiasOpt, NDataSets, MinAge, MaxAge, RefAge, MinusAge, PlusAge, MaxSd, MaxExpectedAge, SaveFile, EffSampleSize=0, Intern=TRUE, AdmbFile=NULL, JustWrite=FALSE){

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
    write("# Option for standard deviation",file=paste(SaveFile,"agemat.dat",sep=""),app
