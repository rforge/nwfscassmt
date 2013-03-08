###################               CHANGE LOG            #############################
######################################################################################
#v 2.4: Eric Ward added user-specified bounds for logit and log functions, 10/2/2012.
#v 2.5: JimT fixed bugs related to diagnostic plots
#v 2.6: JimT fixed bugs related to diagnostic plots
#v 2.7 JimT commented out diagnostic code for MLE comparison and correlation plots because they were not working
#v 2.8: JimT changed processData to not rely upon a hidden attach() command
#v 2.9: JimT changed strataYear for fixed-effects to have the first effect set to zero
#v 2.10: JimT changed vesselYear for fixed-effects to have the first effect set to zero
# v2.11 JimT reverted strataYear to have the first effect estimated, and instead set Year and Strata effects to zero; JimT also fixed warnings in doMCMCDiags() and fixed the RawIndex calculations to work with strata with no positive catch
# v2.12 EricW added an OS check for running in parallel, and included arguments for variance expansion (Gelman et al. 2007)
# v2.13 JimT added displays for processData() function regarding entries that are excluded; also fixed a bug for parallelizing 'prior.scale' option
# v2.14 JimT fixed the ComputeMleIndices() function and added it back to doMCMCDiags(); JimT also added a warning that using non-default catchability settings will cause indices to not have the same units as a Raw (design-based) index
# v2.15 JimT fixed another bug with ComputeMleIndices()
# 1/2/2013: EricW added a csv file for exluded tows in processData() and attached Data to fitted model objects
# 3/5/2013: EricW added an option for fixing the CV of positive tows at 1
############################################################################
# FUNCTION processData WRITTEN BY JIM THORSON & ERIC WARD, UPDATED 9/30/2012. 
# EMAIL: JAMES.THORSON@NOAA.GOV, ERIC.WARD@NOAA.GOV
############################################################################
processData = function() {
  
  print("Necessary column names for masterDat:")
    print("1. BEST_DEPTH_M -> tow depth in meters")
    print("2. BEST_LAT_DD -> tow latitude in degrees ")
    print("3. [insert species name] -> tow catch in kilograms ")
    print("4. YEAR -> calendar year, or time-strata")
    print("5. AREA_SWEPT_MSQ -> area-swept in square meters, or effort offset ")
    print("5. VESSEL -> vessel ID ")
  print("Please ensure that latitude and depth in strata.limits match the following boundaries:")
    print("Latitude: 42-49 in 0.5 increments")
    print("Depth (meters): 55, 75, 100, 125, 155, 183, 200, 250, 300, 350, 400, 450, 500, 549, 600, 700, 800, 900, 1000, 1100, 1200, 1280")

  # Give information about necessary headers
  if(!all(c('BEST_DEPTH_M','BEST_LAT_DD',species,'YEAR','AREA_SWEPT_MSQ','VESSEL')%in%colnames(masterDat))){
    print("Warning: processData() terminated unsuccessfully.")
    print("Please ensure that masterDat has appropriate column names")
    stop()
  }
  
  # set up the generic data frame for this species
  Data = data.frame('PROJECT_CYCLE'=masterDat[,'YEAR'], 'BEST_DEPTH_M'=masterDat[,'BEST_DEPTH_M'], 'BEST_LAT_DD'=masterDat[,'BEST_LAT_DD'], 'HAUL_WT_KG'=masterDat[,which(dimnames(masterDat)[[2]]==species)], 'year'=as.factor(masterDat[,'YEAR']), 'effort'=masterDat[,'AREA_SWEPT_MSQ']*0.0001, 'VESSEL'=masterDat[,'VESSEL']) 
  Data = cbind(Data, 'y'=Data[,'HAUL_WT_KG']) 
  Data = cbind(Data, 'strata'=apply(masterDat,1,strata.fn,Strata.df=strata.limits))
  Data = cbind(Data, 'isNonZeroTrawl'=ifelse(Data[,'y']>0,1,0))
  Data = cbind(Data, 'ones.vec'=rep(0,length(Data[,'y']))) # this is just for the 'ones-trick', inv Gaussian 
  Data = cbind(Data, 'logy3'=log(pi*2*(Data[,'y']^3))) # this is a constant for the invGaussian
  Data = cbind(Data, 'logeffort'=log(Data[,'effort']))
  Data = cbind(Data, 'effort2' = Data[,'effort']^2)
  Data = cbind(Data, 'logeffort2' = Data[,'logeffort']^2)  

  # Exclude tows from strata that are not included
  Exclude = which(is.na(Data$strata))
  print(paste("Excluded ",length(Exclude)," observations that were not assigned to any strata",sep=""))
  if(length(Exclude) > 0){
    print(paste("Observations that were not assigned to any strata are shown in 'Tows_outside_strata.csv'",sep=""))  
    write.table(Data[Exclude,],"Tows_outside_strata.csv",row.names=F,col.names=T,sep=",")
    Data = Data[-Exclude,]
  }
  
  # Exclude tows with some missing entry
  Exclude = which(apply(Data, MARGIN=1, FUN=function(Vec){any(is.na(Vec))}))
  print(paste("Excluded ",length(Exclude)," additional observations that had some missing data",sep=""))
  write.table(Data[Exclude,],"Tows_with_missing_data.csv",row.names=F,col.names=T,sep=",")
  if(length(Exclude) < 10 & length(Exclude) > 0) print(Data[Exclude,])
  if(length(Exclude) >= 10) print("Entries are not printed to the screen due to having 10 or more")
  if(length(Exclude) > 0) Data = Data[-Exclude,]
  
  # Count number of tows per strata and year
  TowsPerStrataYear = table(Data[,'strata'],Data[,'year'])
  print(paste("Tows per strata and year are displayed below"))
  print(TowsPerStrataYear)
  
  # Count number of positive catches per strata and year
  EncountersPerStrataYear = table(Data[,'strata'],Data[,'year'],Data[,'isNonZeroTrawl'])[,,"1"]
  print(paste("Encounters per strata and year are displayed below"))
  print(EncountersPerStrataYear)
  
  # Redefine variables that depend on year
  Data[,'year'] = factor(as.character(Data[,'year'])) # Record year factor in case years have been eliminated due to missing observations
  Data = data.frame(Data, 'vessel'=Letters[as.numeric(as.factor(as.character(Data[,'VESSEL'])))]) # Record year factor in case years have been eliminated due to missing observations
  # Define derived variables involving year
  Data = cbind(Data, 'strataYear'=factor(paste(Data[,'strata'],":",Data[,'year'],sep=""), levels=as.vector(outer(sort(unique(Data[,'strata'])),sort(unique(Data[,'year'])),FUN=paste,sep=":"))))
  Data = cbind(Data, 'vesselYear'=factor(paste(Data[,'vessel'],":",Data[,'year'],sep="")))
  # Attach and calculate other values that aren't in the data.frame
  #attach(Data)
  assign("Data", Data, envir = .GlobalEnv)
  # create null matrices for covariates
  assign("X.pos", matrix(0,1,1), envir = .GlobalEnv)
  assign("X.bin", matrix(0,1,1), envir = .GlobalEnv)  
  #nonZeros = which(Data[,'isNonZeroTrawl']==TRUE)
  assign("nonZeros", which(Data[,'isNonZeroTrawl']==TRUE), envir = .GlobalEnv)
  # Number of elements
  #nS = length(unique(strata))
  assign("nS", length(unique(Data[,'strata'])), envir = .GlobalEnv)
  #nY = length(unique(year))
  assign("nY", length(unique(Data[,'year'])), envir = .GlobalEnv)
  #nNonZeros = length(nonZeros)
  assign("nNonZeros", length(nonZeros), envir = .GlobalEnv)
  #n = length(y)
  assign("n", length(Data[,'y']), envir = .GlobalEnv)
  #nSY = length(unique(strataYear)) 
  assign("nSY", nlevels(Data[,'strataYear']), envir = .GlobalEnv)
  #nVY = length(unique(vesselYear))
  assign("nVY", length(unique(Data[,'vesselYear'])), envir = .GlobalEnv)  
  # Diagonal matrix for the wishart / correlation model
  assign("R",diag(2), envir = .GlobalEnv)
}

############################################################################
# FUNCTION fitCPUEModel WRITTEN BY ERIC WARD, UPDATED 9/30/2012. 
# EMAIL: ERIC.WARD@NOAA.GOV
############################################################################
# This function writes, and optinally runs the model
# in the BUGS/JAGS framework. The model options are 
# StrataYear.positiveTows = "random" (default), "fixed", not included = "zero", or correlated with zero tows = "correlated"
# StrataYear.zeroTows = "random" (default), "fixed", not included = "zero", or correlated with positive tows = "correlated"
# VesselYear.positiveTows = "random" (default), "fixed", not included = "zero", or correlated with zero tows = "correlated"
# VesselYear.zeroTows = "random" (default), "fixed", or not included = "zero", or correlated with positive tows = "correlated"
# Catchability.positiveTows = can be 1 = "one", or estimated as "linear" (default), or "quadratic"
# Catchability.zeroTows = can be 1 = "one", 0 = "zero" (default), or estimated as "linear", or "quadratic"
# year.deviations = "uncorrelated" (default) or "correlated"
# strata.deviations = "uncorrelated" (default) or "correlated"
# likelihood = "gamma" (default), "lognormal", "invGaussian", "lognormalECE", "gammaECE"
# model.name = "deltaGLM.txt", the name of this model file
# fit.model = can be T/F, if FALSE, the model file is written to a txt file
# mcmc.control$chains = MCMC chains
# mcmc.control$thin = MCMC thinning rate
# mcmc.control$burn = MCMC burn-in
# mcmc.control$iterToSave = MCMC samples post-burn
# Parallel = TRUE (default) or FALSE
# Species = string with species name, required
# logitBounds = bounds for logit space devs, default = c(-5,5)
# logBounds = bounds for log-link space devs, default = c(-5,5)

# For MCMC samples, the total number of iterations returned will be (chains * iterToSave) / thin rate
###########################################################################################
fitCPUEModel = function(modelStructure = list("StrataYear.positiveTows" = "random","VesselYear.positiveTows" = "random","StrataYear.zeroTows" ="random","VesselYear.zeroTows" = "random", "Catchability.positiveTows" = "one", "Catchability.zeroTows" = "zero", "year.deviations" = "uncorrelated","strata.deviations" = "uncorrelated"),covariates=list(positive=FALSE,binomial=FALSE),likelihood = "gamma", model.name = "deltaGLM.txt", fit.model=TRUE, mcmc.control = list(chains = 5, thin = 1, burn = 5000, iterToSave = 2000),Parallel=TRUE, Species = "NULL",logitBounds = c(-20,20),logBounds = c(-20,20), prior.scale = rep(25,4)) {

  if(modelStructure$Catchability.positiveTows%in%c("linear","quadratic") | modelStructure$Catchability.zeroTows%in%c("one","linear","quadratic")){
    print("Warning: index will not have comparable scale to a design-based (raw) index unless catchability.positiveTows equals 'one'  catchability.zeroTows equals 'zero'") 
  }
  
  if(.Platform$OS.type != "windows" & Parallel == TRUE) {
  	print("Warning: system OS is non-windows, and running in parallel is currently not supported. Code will run in non-parallel")
  	Parallel = FALSE
  }
  if(length(prior.scale) != 4 | length(which(is.na(prior.scale)==T)) > 0 | length(which(prior.scale <= 0)) > 0) {
  	print("Error: prior.scale needs to be specified as a 4-element vector, of positive values")
  	stop()
  }
  if(Species == "NULL") {
  	print("Error: you need to input a species name in the function, e.g. Species = \"Canary\"")
  	stop()
  }
  if(modelStructure$VesselYear.zeroTows =="correlated" & modelStructure$VesselYear.positiveTows != "correlated") {
  	print("Error: you specified only Vessel Year deviations for binomial model to be correlated. You need to either (1) specify both as correlated or (2) specify neither as correlated")
  	stop()	
  }
  if(modelStructure$VesselYear.zeroTows !="correlated" & modelStructure$VesselYear.positiveTows == "correlated") {
  	print("Error: you specified only Vessel Year deviations for positive model to be correlated. You need to either (1) specify both as correlated or (2) specify neither as correlated")
  	stop()	
  }
  if(modelStructure$StrataYear.zeroTows =="correlated" & modelStructure$StrataYear.positiveTows != "correlated") {
  	print("Error: you specified only Strata Year deviations for binomial model to be correlated. You need to either (1) specify both as correlated or (2) specify neither as correlated")
  	stop()	
  }
  if(modelStructure$StrataYear.zeroTows !="correlated" & modelStructure$StrataYear.positiveTows == "correlated") {
  	print("Error: you specified only Strata Year deviations for positive model to be correlated. You need to either (1) specify both as correlated or (2) specify neither as correlated")
  	stop()	
  }
  
  if(nVY == nY) {
  	# EW: this was a catch I put in for the darkblotched data on 5.13.12. 
  	# When there's only 1 vessel, and VY interaction needs to be set to 0
  	# because it contains redundant info as year
  	modelStructure$VesselYear.positiveTows == "zero"
  	modelStructure$VesselYear.zeroTows == "zero"
  }
  
  # check to make sure that the covariates are in matrix form - this is just a BUGS/JAGS thing
  nX.binomial = 1
  if(covariates$binomial) {
  	# check for matrix
  	if(is.matrix(X.bin)==FALSE) {
  		print("Error: Please specify the covariates for the binomial model as a matrix, or turn covariates off.")
  	  stop()
  	}
  	nX.binomial = dim(X.bin)[2]
  }
  nX.pos = 1
  if(covariates$positive) {
  	# check for matrix
  	if(is.matrix(X.pos)==FALSE) {
  		print("Error: Please specify the covariates for positive tows as a matrix, or turn covariates off.")
  	  stop()
  	}
  	nX.pos = dim(X.pos)[2]	
  }
  
  #################################################################
  # These 4 strings are all new, relevant for implementing the variance manipulated/expanded method described in Gelman et al. 2007, 
  # Gelman et al. 2006. The prior implicit on random effects is the non-central folded t distribution from Gelman et al. 2006
  SYexpanded = paste("   tau.xi[1] <- pow(",prior.scale[1],", -2);\n","   xi[1] ~ dnorm (0, tau.xi[1]);\n",
  "tau.eta[1] ~ dgamma(0.5, 0.5); # chi^2 with 1 d.f.\n","sigmaSY[1] <- abs(xi[1])/sqrt(tau.eta[1]); # derived, cauchy = normal/sqrt(chi^2)\n",
  "for (j in 1:nSY){\n",
  "   SYeta[j] ~ dnorm(0, tau.eta[1]); # hierarchical model for theta\n",
  "   SYdev[j] <- min(max(xi[1]*SYeta[j],",logBounds[1],"),",logBounds[2],");\n",
  "}\n",sep="")
  pSYexpanded = paste("   tau.xi[2] <- pow(",prior.scale[2],", -2);\n","   xi[2] ~ dnorm (0, tau.xi[2]);\n",
  "tau.eta[2] ~ dgamma(0.5, 0.5); # chi^2 with 1 d.f.\n","sigmaSY[2] <- abs(xi[2])/sqrt(tau.eta[2]); # derived, cauchy = normal/sqrt(chi^2)\n",
  "for (j in 1:nSY){\n",
  "   pSYeta[j] ~ dnorm(0, tau.eta[2]); # hierarchical model for theta\n",
  "   pSYdev[j] <- min(max(xi[2]*pSYeta[j],",logitBounds[1],"),",logitBounds[2],");\n",
  "}\n",sep="")
  VYexpanded = paste("   tau.xi[3] <- pow(",prior.scale[3],", -2);\n","   xi[3] ~ dnorm (0, tau.xi[3]);\n",
  "tau.eta[3] ~ dgamma(0.5, 0.5); # chi^2 with 1 d.f.\n","sigmaVY[1] <- abs(xi[3])/sqrt(tau.eta[3]); # derived, cauchy = normal/sqrt(chi^2)\n",
  "for (j in 1:nVY){\n",
  "   VYeta[j] ~ dnorm(0, tau.eta[3]); # hierarchical model for theta\n",
  "   VYdev[j] <- min(max(xi[3]*VYeta[j],",logBounds[1],"),",logBounds[2],");\n",
  "}\n",sep="")
  pVYexpanded = paste("   tau.xi[4] <- pow(",prior.scale[4],", -2);\n","   xi[4] ~ dnorm (0, tau.xi[4]);\n",
  "tau.eta[4] ~ dgamma(0.5, 0.5); # chi^2 with 1 d.f.\n","sigmaVY[2] <- abs(xi[4])/sqrt(tau.eta[4]); # derived, cauchy = normal/sqrt(chi^2)\n",
  "for (j in 1:nVY){\n",
  "   pVYeta[j] ~ dnorm(0, tau.eta[4]); # hierarchical model for theta\n",
  "   pVYdev[j] <- min(max(xi[4]*pVYeta[j],",logitBounds[1],"),",logitBounds[2],");\n",
  "}\n",sep="")
  
  ####################################################################
  # This section is related to strata-year and vessel-year effects
  # Strata-year interactions can be (1) fixed, (2) random, (3) randomExpanded, or (4) not estimated (set to 0)
  SYpos.string = ""
  if(modelStructure$StrataYear.positiveTows == "fixed") SYpos.string = paste("   for(i in 1:nSY) {\n      SYdev[i] ~ dunif(",logBounds[1],",",logBounds[2],");\n   }\n   strataYearTau[1,1] <- 0;\n","   strataYearTau[1,2] <- 0;\n   sigmaSY[1]<-0;\n   tauSY[1]<-0;\n",sep="")
  if(modelStructure$StrataYear.positiveTows == "random") SYpos.string = paste("   for(i in 1:nSY) {\n      SYdev[i] ~ dnorm(0,tauSY[1])T(",logBounds[1],",",logBounds[2],");\n   }\n   strataYearTau[1,1] <- 0;\n","   strataYearTau[1,2] <- 0;\n   sigmaSY[1]~dunif(0,100);\n   tauSY[1]<-pow(sigmaSY[1],-2);\n",sep="")
  if(modelStructure$StrataYear.positiveTows == "randomExpanded") SYpos.string = paste(SYexpanded, "   strataYearTau[1,1] <- 0;\n","   strataYearTau[1,2] <- 0;\n   tauSY[1]<-pow(sigmaSY[1],-2);\n",sep="")
  if(modelStructure$StrataYear.positiveTows == "zero") SYpos.string = "   for(i in 1:nSY) {\n      SYdev[i] <- 0;\n   }\n   strataYearTau[1,1] <- 0;\n   strataYearTau[1,2] <- 0;\n   sigmaSY[1]<-0;\n   tauSY[1]<-0;\n"
  
  SYzero.string = ""
  if(modelStructure$StrataYear.zeroTows == "fixed") SYzero.string = paste("   for(i in 1:nSY) {\n      pSYdev[i] ~ dunif(",logitBounds[1],",",logitBounds[2],");\n   }\n   strataYearTau[2,1] <- 0;\n","   strataYearTau[2,2] <- 0;\n","   sigmaSY[2] <- 0;\n   tauSY[2]<-0;\n",sep="")
  if(modelStructure$StrataYear.zeroTows == "random") SYzero.string = paste("   for(i in 1:nSY) {\n      pSYdev[i] ~ dnorm(0,tauSY[2])T(",logitBounds[1],",",logitBounds[2],");\n   }\n   strataYearTau[2,1] <- 0;\n","   strataYearTau[2,2] <- 0;\n   sigmaSY[2]~dunif(0,100);\n   tauSY[2]<-pow(sigmaSY[2],-2);\n",sep="")
  if(modelStructure$StrataYear.zeroTows == "randomExpanded") SYzero.string = paste(pSYexpanded, "   strataYearTau[2,1] <- 0;\n","   strataYearTau[2,2] <- 0;\n   tauSY[2]<-pow(sigmaSY[2],-2);\n",sep="")
  if(modelStructure$StrataYear.zeroTows == "zero") SYzero.string = "   for(i in 1:nSY) {\n      pSYdev[i] <- 0;\n   }\n   strataYearTau[2,1] <- 0;\n   strataYearTau[2,2] <- 0;\n   sigmaSY[2] <- 0;\n   tauSY[2]<-0;\n"
  
  # combine the strata year interactions into a string
  stratayear.string = paste(SYpos.string,SYzero.string)
  if(modelStructure$StrataYear.zeroTows == "correlated" & modelStructure$StrataYear.positiveTows == "correlated") {
  	# strata year deviations are MVN RE
  	stratayear.string = paste("sigmaSY[1] <- 0;\n   sigmaSY[2] <- 0;\n      strataYearTau[1:2,1:2] ~ dwish(R[1:2,1:2],2);\n   for(i in 1:nSY) {\n   sydevs[i,1:2] ~ dmnorm(zs[1:2],strataYearTau[1:2,1:2]);\n      SYdev[i] <- min(max(sydevs[i,1],",logBounds[1],"),",logBounds[2],");\n      pSYdev[i] <- min(max(sydevs[i,2],",logitBounds[1],"),",logitBounds[2],");\n   }\n",sep="")	
  }
  
  # Vessel-year interactions cab be (1) fixed, (2) random, (3) randomExpanded, or (4) not estimated (set to 0) 
  VYpos.string = ""
  if(modelStructure$VesselYear.positiveTows == "fixed") VYpos.string = paste("   VYdev[1] <- 0;\n   for(i in 2:nVY) {\n      VYdev[i] ~ dunif(",logBounds[1],",",logBounds[2],");\n   }\n   vesselYearTau[1,1] <- 0;\n","   vesselYearTau[1,2] <- 0;\n   sigmaVY[1]<-0;\n   tauVY[1]<-0;\n",sep="")
  if(modelStructure$VesselYear.positiveTows == "random") VYpos.string = paste("   for(i in 1:nVY) {\n      VYdev[i] ~ dnorm(0,tauVY[1])T(",logBounds[1],",",logBounds[2],");\n   }\n   vesselYearTau[1,1] <- 0;\n","   vesselYearTau[1,2] <- 0;\n   sigmaVY[1]~dunif(0,100);\n   tauVY[1]<-pow(sigmaVY[1],-2);\n",sep="")
  if(modelStructure$VesselYear.positiveTows == "randomExpanded") VYpos.string = paste(VYexpanded,"   vesselYearTau[1,1] <- 0;\n","   vesselYearTau[1,2] <- 0;\n   tauVY[1]<-pow(sigmaVY[1],-2);\n",sep="")
  if(modelStructure$VesselYear.positiveTows == "zero") VYpos.string = "   for(i in 1:nVY) {\n      VYdev[i] <- 0;\n   }\n   vesselYearTau[1,1] <- 0;\n   vesselYearTau[1,2] <- 0;\n   sigmaVY[1]<-0;\n   tauVY[1]<-0;\n"
  
  VYzero.string = ""
  if(modelStructure$VesselYear.zeroTows == "fixed") VYzero.string = paste("   pVYdev[1] <- 0;\n   for(i in 2:nVY) {\n      pVYdev[i] ~ dunif(",logitBounds[1],",",logitBounds[2],");\n   }\n   vesselYearTau[2,1] <- 0;\n","   vesselYearTau[2,2] <- 0;\n   sigmaVY[2]<-0;\n   tauVY[2]<-0;\n",sep="")
  if(modelStructure$VesselYear.zeroTows == "random") VYzero.string = paste("   for(i in 1:nVY) {\n      pVYdev[i] ~ dnorm(0,tauVY[2])T(",logitBounds[1],",",logitBounds[2],");\n   }\n   vesselYearTau[2,1] <- 0;\n","   vesselYearTau[2,2] <- 0;\n   sigmaVY[2] ~ dunif(0,100);\n   tauVY[2]<-pow(sigmaVY[2],-2);\n",sep="")
  if(modelStructure$VesselYear.zeroTows == "randomExpanded") VYzero.string = paste(pVYexpanded,"   vesselYearTau[2,1] <- 0;\n","   vesselYearTau[2,2] <- 0;\n   tauVY[2]<-pow(sigmaVY[2],-2);\n",sep="")
  if(modelStructure$VesselYear.zeroTows == "zero") VYzero.string = "   for(i in 1:nVY) {\n      pVYdev[i] <- 0;\n   }\n   vesselYearTau[2,1] <- 0;\n   vesselYearTau[2,2] <- 0;\n   sigmaVY[2]<-0;\n   tauVY[2]<-0;\n"
  
  # combine the strata year interactions into a string
  vesselyear.string = paste(VYpos.string,VYzero.string)
  #vesselyear.string = paste("   for(i in 1:nVY) {\n",VYpos.string,VYzero.string,"   }\n","   vesselYearTau[1,1] <- 0;\n","   vesselYearTau[1,2] <- 0;\n","   vesselYearTau[2,1] <- 0;\n","   vesselYearTau[2,2] <- 0;\n")

  if(modelStructure$VesselYear.zeroTows == "correlated" & modelStructure$VesselYear.positiveTows == "correlated") {
  	# vessel year deviations are MVN RE
  	vesselyear.string = paste("sigmaVY[1] <- 0;\n   sigmaVY[2] <- 0;\n   vesselYearTau[1:2,1:2] ~ dwish(R[1:2,1:2],2);\n   for(i in 1:nVY) {\n   vydevs[i,1:2] ~ dmnorm(zs[1:2],vesselYearTau[1:2,1:2]);\n      VYdev[i] <- min(max(vydevs[i,1],",logBounds[1],"),",logBounds[2],");\n      pVYdev[i] <- min(max(vydevs[i,2],",logitBounds[1],"),",logitBounds[2],");\n   }\n",sep="")	
  }
  
  ####################################################################
  # This section is related to offsets, separate for the positive and binomial models
  # catchability parameter can be set to 1 or estimated as linear or qudratic
  #if(modelStructure$Catchability.positiveTows=="zero") catch.posTows = "   B.pos[1]<-0;\n   B.pos[2] <- 0;\n"
  if(modelStructure$Catchability.positiveTows=="one") catch.posTows = "   logB.pos[1] <- 0;\n   B.pos[1]<-exp(logB.pos[1]);\n   B.pos[2] <- 0;\n"
  if(modelStructure$Catchability.positiveTows=="linear") catch.posTows = "   logB.pos[1] ~ dnorm(0,0.1);\n   B.pos[1]<-exp(logB.pos[1]);\n   B.pos[2] <- 0;\n"
  if(modelStructure$Catchability.positiveTows=="quadratic") catch.posTows = "   logB.pos[1] ~ dnorm(0,0.1);\n   B.pos[1]<-logB.pos[1];\n   logB.pos[2] ~ dnorm(0,0.1);\n   B.pos[2]<-logB.pos[2];\n"
  
  if(modelStructure$Catchability.zeroTows=="zero") catch.zeroTows = "   B.zero[1]<-0;\n   B.zero[2] <- 0;\n"
  if(modelStructure$Catchability.zeroTows=="one") catch.zeroTows = "   logB.zero[1] <- 0;\n   B.zero[1] <- exp(logB.zero[1]);\n   B.zero[2] <- 0;\n"
  if(modelStructure$Catchability.zeroTows=="linear") catch.zeroTows = "   logB.zero[1] ~ dnorm(0,0.1);\n   B.zero[1] <- exp(logB.zero[1]);\n   B.zero[2] <- 0;\n"
  if(modelStructure$Catchability.zeroTows=="quadratic") catch.zeroTows = "   logB.zero[1] ~ dnorm(0,0.1);\n   B.zero[1]<-logB.zero[1];\n   logB.zero[2] ~ dnorm(0,0.1);\n   B.zero[2]<-logB.zero[2];\n"
  
  ####################################################################
  # This section is related to likelihoods
  # likelihood: lognormal, gamma, inverse gaussian
  if(likelihood == "lognormal" | likelihood == "lognormalFixedCV") {
  	# CV is sqrt(exp(sig2)-1) which is approx ~ sigma for sigma is small (< 0.2)
  	likelihood.string = paste("      u.nz[i] <- Sdev[strata[nonZeros[i]]] + Ydev[year[nonZeros[i]]] + VYdev[vesselYear[nonZeros[i]]] + SYdev[strataYear[nonZeros[i]]] + B.pos[1]*logeffort[nonZeros[i]] + B.pos[2]*logeffort2[nonZeros[i]];\n", "      y[nonZeros[i]] ~ dlnorm(u.nz[i],tau[1]);\n",sep="")
  	
  	if(covariates$positive==TRUE) {
    	# modify the likelihood string to include covariates for the positive tows
    	likelihood.string = paste("      u.nz[i] <- inprod(C.pos[1:nX.pos],X.pos[i,]) + Sdev[strata[nonZeros[i]]] + Ydev[year[nonZeros[i]]] + VYdev[vesselYear[nonZeros[i]]] + SYdev[strataYear[nonZeros[i]]] + B.pos[1]*logeffort[nonZeros[i]] + B.pos[2]*logeffort2[nonZeros[i]];\n", "      y[nonZeros[i]] ~ dlnorm(u.nz[i],tau[1]);\n",sep="")
  	}
  	
  	prior.string = "   oneOverCV2[1] ~ dgamma(0.001,0.001);\n   CV[1] <- 1/sqrt(oneOverCV2[1]);\n   CV[2] <- 0;\n   sigma[1] <- sqrt(log(pow(CV[1],2)+1));\n   tau[1] <- pow(sigma[1],-2);\n   ratio <- 0;\n   p.ece[1] <- 0;\n   p.ece[2] <- 0;\n"
  	if(likelihood == "lognormalFixedCV") {
  		prior.string = "   oneOverCV2[1] <- 1;\n   CV[1] <- 1;\n   CV[2] <- 0;\n   sigma[1] <- 1;\n   tau[1] <- 1;\n   ratio <- 0;\n   p.ece[1] <- 0;\n   p.ece[2] <- 0;\n"	
  	}
  	
  }
  if(likelihood == "gamma" | likelihood == "gammaFixedCV") {
  	# gamma in this instance is parameterized in terms of the rate and shape, with mean = a/b, var = a/b2, and CV = 1/sqrt(a)
  	# So parameter 'a' has to be a constant, and 'b' varies by tows - b/c we calculate b = a/mean
    likelihood.string = paste("      u.nz[i] <- exp(min(Sdev[strata[nonZeros[i]]] + Ydev[year[nonZeros[i]]] + VYdev[vesselYear[nonZeros[i]]] + SYdev[strataYear[nonZeros[i]]] + B.pos[1]*logeffort[nonZeros[i]] + B.pos[2]*logeffort2[nonZeros[i]],100));\n","      b[i] <- gamma.a[1]/u.nz[i];\n","      y[nonZeros[i]] ~ dgamma(gamma.a[1],b[i]);\n")
      
  	if(covariates$positive==TRUE) {
    	# modify the likelihood string to include covariates for the positive tows
    	likelihood.string = paste("      u.nz[i] <- exp(min(inprod(C.pos[1:nX.pos],X.pos[i,]) + Sdev[strata[nonZeros[i]]] + Ydev[year[nonZeros[i]]] + VYdev[vesselYear[nonZeros[i]]] + SYdev[strataYear[nonZeros[i]]] + B.pos[1]*logeffort[nonZeros[i]] + B.pos[2]*logeffort2[nonZeros[i]],100));\n","      b[i] <- gamma.a[1]/u.nz[i];\n","      y[nonZeros[i]] ~ dgamma(gamma.a[1],b[i]);\n")
  	}    
  
    # for lognormal, gamma prior on 1/sigma2 = 1/CV2. To keep things consistent, gamma prior on a = 1/cv2, b/c CV = 1/sqrt(a)
    # then gamma.b[i] = gamma.a / u[i]
  	prior.string = "   oneOverCV2[1] ~ dgamma(0.001,0.001);\n   gamma.a[1] <- oneOverCV2[1];\n   CV[1] <- 1/sqrt(oneOverCV2[1]);\n   CV[2] <- 0;\n   ratio <- 0;\n   p.ece[1] <- 0;\n   p.ece[2] <- 0;\n"
  	if(likelihood == "gammaFixedCV") {
  	prior.string = "   oneOverCV2[1] <- 1;\n   gamma.a[1] <- oneOverCV2[1];\n   CV[1] <- 1/sqrt(oneOverCV2[1]);\n   CV[2] <- 0;\n   ratio <- 0;\n   p.ece[1] <- 0;\n   p.ece[2] <- 0;\n"  		
  	}	
  }
  if(likelihood == "invGaussian" | likelihood == "invGaussianFixedCV") {
  	# again, parameterize in terms of CV
  	likelihood.string = paste("      u.nz[i] <- exp(min(Sdev[strata[nonZeros[i]]] + Ydev[year[nonZeros[i]]] + VYdev[vesselYear[nonZeros[i]]] + SYdev[strataYear[nonZeros[i]]] + B.pos[1]*logeffort[nonZeros[i]] + B.pos[2]*logeffort2[nonZeros[i]],100));\n","      lambda[i] <- u.nz[i]*oneOverCV2[1];\n","      scaledLogLike[i] <- -(0.5*log(lambda[i]) - 0.5*logy3[nonZeros[i]] - lambda[i]*pow((y[nonZeros[i]]-u.nz[i]),2)/(2*u.nz[i]*u.nz[i]*y[nonZeros[i]])) + 10000;\n","      ones.vec[i] ~ dpois(scaledLogLike[i]);\n")
  	
  	if(covariates$positive==TRUE) {
    	# modify the likelihood string to include covariates for the positive tows
    	likelihood.string = paste("      u.nz[i] <- exp(min(inprod(C.pos[1:nX.pos],X.pos[i,]) + Sdev[strata[nonZeros[i]]] + Ydev[year[nonZeros[i]]] + VYdev[vesselYear[nonZeros[i]]] + SYdev[strataYear[nonZeros[i]]] + B.pos[1]*logeffort[nonZeros[i]] + B.pos[2]*logeffort2[nonZeros[i]],100));\n","      lambda[i] <- u.nz[i]*oneOverCV2[1];\n","      scaledLogLike[i] <- -(0.5*log(lambda[i]) - 0.5*logy3[nonZeros[i]] - lambda[i]*pow((y[nonZeros[i]]-u.nz[i]),2)/(2*u.nz[i]*u.nz[i]*y[nonZeros[i]])) + 10000;\n","      ones.vec[i] ~ dpois(scaledLogLike[i]);\n")
  	} 	
  	
  	prior.string = "   oneOverCV2[1] ~ dgamma(0.001,0.001);\n   CV[1] <- 1/sqrt(oneOverCV2[1]);\n   CV[2] <- 0;\n   ratio <- 0;\n   p.ece[1] <- 0;\n   p.ece[2] <- 0;\n"	
  	if(likelihood == "invGaussianFixedCV") prior.string = "   oneOverCV2[1] <- 1;\n   CV[1] <- 1/sqrt(oneOverCV2[1]);\n   CV[2] <- 0;\n   ratio <- 0;\n   p.ece[1] <- 0;\n   p.ece[2] <- 0;\n"
  }
  
  if(likelihood == "lognormalECE") {
  	# CV is sqrt(exp(sig2)-1) which is approx ~ sigma. Each tow is treated as discrete group (normal, ECE)
  	
  	# if 1, don't do anything 
  	# if 2, multiply by ratio 
  	likelihood.string = paste("      G[i] ~ dcat(p.ece[1:2]);\n      u.nz[i] <- (G[i]-1)*logratio + (Sdev[strata[nonZeros[i]]] + Ydev[year[nonZeros[i]]] + VYdev[vesselYear[nonZeros[i]]] + SYdev[strataYear[nonZeros[i]]] + B.pos[1]*logeffort[nonZeros[i]] + B.pos[2]*logeffort2[nonZeros[i]]);\n", "      y[nonZeros[i]] ~ dlnorm(u.nz[i],tau[G[i]]);\n",sep="")
  	
  	if(covariates$positive==TRUE) {       
    	# modify the likelihood string to include covariates for the positive tows
    	likelihood.string = paste("      G[i] ~ dcat(p.ece[1:2]);\n      u.nz[i] <- (G[i]-1)*logratio + (inprod(C.pos[1:nX.pos],X.pos[i,]) + Sdev[strata[nonZeros[i]]] + Ydev[year[nonZeros[i]]] + VYdev[vesselYear[nonZeros[i]]] + SYdev[strataYear[nonZeros[i]]] + B.pos[1]*logeffort[nonZeros[i]] + B.pos[2]*logeffort2[nonZeros[i]]);\n", "      y[nonZeros[i]] ~ dlnorm(u.nz[i],tau[G[i]]);\n",sep="")
  	} 		
  	
  	# for the ECE model, the normal and extreme distributions each get a separate variance
  	prior.string = "   oneOverCV2[1] ~ dgamma(0.001,0.001);\n   CV[1] <- 1/sqrt(oneOverCV2[1]);\n   sigma[1] <- sqrt(log(pow(CV[1],2)+1));\n   tau[1] <- pow(sigma[1],-2);\n   oneOverCV2[2] ~ dgamma(0.001,0.001);\n   CV[2] <- 1/sqrt(oneOverCV2[2]);\n   sigma[2] <- sqrt(log(pow(CV[2],2)+1));\n   tau[2] <- pow(sigma[2],-2);\n   logratio ~ dunif(0,5);\n   ratio <- exp(logratio);\n   alpha.ece[1] <- 1;\n   alpha.ece[2] <- 1;\n   p.ece[1:2] ~ ddirch(alpha.ece[1:2]);\n"
  }
  if(likelihood == "gammaECE") {
  	# gamma in this instance is parameterized in terms of the rate and shape, with mean = a/b, var = a/b2, and CV = 1/sqrt(a)
  	# So parameter 'a' has to be a constant, and 'b' varies by tows - b/c we calculate b = a/mean
    likelihood.string = paste("      G[i] ~ dcat(p.ece[1:2]);\n      u.nz[i] <- exp((G[i]-1)*logratio + min(Sdev[strata[nonZeros[i]]] + Ydev[year[nonZeros[i]]] + VYdev[vesselYear[nonZeros[i]]] + SYdev[strataYear[nonZeros[i]]] + B.pos[1]*logeffort[nonZeros[i]] + B.pos[2]*logeffort2[nonZeros[i]],100));\n","      b[i] <- gamma.a[G[i]]/u.nz[i];\n","      y[nonZeros[i]] ~ dgamma(gamma.a[G[i]],b[i]);\n")
      
  	if(covariates$positive==TRUE) {
  	# modify the likelihood string to include covariates for the positive tows
    	likelihood.string = paste("      G[i] ~ dcat(p.ece[1:2]);\n      u.nz[i] <- exp((G[i]-1)*logratio + min(inprod(C.pos[1:nX.pos],X.pos[i,]) + Sdev[strata[nonZeros[i]]] + Ydev[year[nonZeros[i]]] + VYdev[vesselYear[nonZeros[i]]] + SYdev[strataYear[nonZeros[i]]] + B.pos[1]*logeffort[nonZeros[i]] + B.pos[2]*logeffort2[nonZeros[i]],100));\n","      b[i] <- gamma.a[G[i]]/u.nz[i];\n","      y[nonZeros[i]] ~ dgamma(gamma.a[G[i]],b[i]);\n")
  	} 		    
      
      # for lognormal, gamma prior on 1/sigma2 = 1/CV2. To keep things consistent, gamma prior on a = 1/cv2, b/c CV = 1/sqrt(a)
      # then gamma.b[i] = gamma.a / u[i]
  	prior.string = "   oneOverCV2[1] ~ dgamma(0.001,0.001);\n   gamma.a[1] <- oneOverCV2[1];\n   CV[1] <- 1/sqrt(oneOverCV2[1]);\n   oneOverCV2[2] ~ dgamma(0.001,0.001);\n   gamma.a[2] <- oneOverCV2[2];\n   CV[2] <- 1/sqrt(oneOverCV2[2]);\n   logratio ~ dunif(0,5);\n   ratio <- exp(logratio);\n   alpha.ece[1] <- 1;\n   alpha.ece[2] <- 1;\n   p.ece[1:2] ~ ddirch(alpha.ece[1:2]);\n"
  }
  
  ####################################################################
  # This section is related to year deviations, default is to make them uncorrelated
  # if strata-year effects are estimated as fixed effects, parameters are redundant, so year deviations set to 0
  str1 = paste("      Ydev[i] ~ dunif(",logBounds[1],",",logBounds[2],");\n",sep="")
  if(modelStructure$StrataYear.positiveTows == "fixed") {str1 = "      Ydev[i] <- 0;\n"}
  str2 = paste("      pYdev[i] ~ dunif(",logitBounds[1],",",logitBounds[2],");\n",sep="")
  if(modelStructure$StrataYear.zeroTows == "fixed") {str2 = "      pYdev[i] <- 0;\n"}
  
  year.dev.string = paste("   for(i in 1:nY) { # year deviations, always fixed \n",str1,str2,"}\n   yearTau[1,1] <- 0;\n   yearTau[1,2] <- 0;\n   yearTau[2,1]<-0;\n   yearTau[2,2] <- 0;\n",sep="")
  if(modelStructure$year.deviations == "correlated") {
  	year.dev.string = paste("   yearTau[1:2,1:2] ~ dwish(R[1:2,1:2],2);\n   for(i in 1:nY) {\n   ydevs[i,1:2] ~ dmnorm(zs[1:2],yearTau[1:2,1:2]);\n      Ydev[i] <- min(max(ydevs[i,1],",logBounds[1],"),",logBounds[1],");\n      pYdev[i] <- min(max(ydevs[i,2],",logitBounds[1],"),",logitBounds[2],");\n   }\n",sep="")
  }
  
  ####################################################################
  # This section is related to strata deviations, default is to make them uncorrelated
  # if strata-year effects are estimated as fixed effects, parameters are redundant, so strata deviations set to 0
  str1 = paste("      Sdev[i] ~ dunif(",logBounds[1],",",logBounds[2],");\n",sep="")
  if(modelStructure$StrataYear.positiveTows == "fixed") {str1 = "      Sdev[i] <- 0;\n"}
  str2 = paste("      pSdev[i] ~ dunif(",logitBounds[1],",",logitBounds[2],");\n",sep="")
  if(modelStructure$StrataYear.zeroTows == "fixed") {str2 = "      pSdev[i] <- 0;\n"}
  strata.dev.string = paste("   for(i in 2:nS) { # strata deviations, always fixed \n", str1, str2,"}\n   strataTau[1,1] <- 0;\n   strataTau[1,2] <- 0;\n   strataTau[2,1]<-0;\n   strataTau[2,2] <- 0;\n")
  if(modelStructure$strata.deviations == "correlated") {
  	strata.dev.string = paste("   strataTau[1:2,1:2] ~ dwish(R[1:2,1:2],2);\n   for(i in 2:nS) {\n   sdevs[i,1:2] ~ dmnorm(zs[1:2],strataTau[1:2,1:2]);\n      Sdev[i] <- min(max(sdevs[i,1],",logBounds[1],"),",logBounds[2],");\n      pSdev[i] <- min(max(sdevs[i,2],",logitBounds[1],"),",logitBounds[2],");\n   }\n",sep="")
  }
  
  # This set of if() blocks is used to include optional strings for covariates...these are all the priors, vague normal
  covarString = ""
  if(covariates$binomial == FALSE) {
  	covarString = paste(covarString, "	C.bin[1] <- 0;\n",sep="")	
  } 
  if(covariates$binomial == TRUE) {
  	covarString = paste(covarString, "	for(i in 1:nX.binomial) {\n","		C.bin[i] ~ dnorm(0,0.001);\n","	}\n",sep="")		
  }
  if(covariates$positive == FALSE) {
  	covarString = paste(covarString, "	C.pos[1] <- 0;\n",sep="")	
  } 
  if(covariates$positive == TRUE) {
  	covarString = paste(covarString, "	for(i in 1:nX.pos) {\n","		C.pos[i] ~ dnorm(0,0.001);\n","	}\n",sep="")		
  }
  
  # This IF statement is just to modify the expected value of the binomial model to include covariates IF they are specified 
  logit.p.string = "logit(p.z[i]) <- pVYdev[vesselYear[i]] + pSdev[strata[i]] + pYdev[year[i]] + pSYdev[strataYear[i]] + B.zero[1]*effort[i] + B.zero[2]*effort2[i];\n"
  if(covariates$binomial == TRUE) {
  	logit.p.string = "logit(p.z[i]) <- inprod(C.bin[1:nX.binomial],X.bin[i,]) + pVYdev[vesselYear[i]] + pSdev[strata[i]] + pYdev[year[i]] + pSYdev[strataYear[i]] + B.zero[1]*effort[i] + B.zero[2]*effort2[i];\n"	
  }
  
  deltaGLM = 
  paste(
  "
  model {\n",
     prior.string,
     covarString,	
     
  "   # next deal with catchability parameters\n",
     catch.posTows,
     catch.zeroTows,
  "  # these zeros are optionally for MVN random effects
     zs[1] <- 0; 
     zs[2] <- 0;
     # first stratum is set to 0 for identifiability
     Sdev[1] <- 0;
     pSdev[1] <- 0;\n
  ",
  year.dev.string,
  strata.dev.string,   
  stratayear.string,   
  vesselyear.string,   	    	    
  "  # for each group, calculate the probability of zero tows and the mean
     # evaluate the likelihood  of non-zero trawls   
     for(i in 1:nNonZeros) {\n",
     	likelihood.string,
  "   }
     
     # evaluate the likelihood  of zero / non-zero trawls
     for(i in 1:n) {
     	 # group represents which vessel x strata combo\n",
  logit.p.string,
  "     isNonZeroTrawl[i] ~ dbern(p.z[i]);
     }
  }
  ",sep="")
  # write this to text file
  cat(deltaGLM, file = model.name)
  
  # fit the model and return the model object
  modelFit = NA
  
  if(fit.model) {
jags.params=c("Ydev","Sdev","SYdev","VYdev","pYdev","pSdev","pSYdev","pVYdev","B.zero","B.pos","sigmaSY","sigmaVY","CV","ratio","p.ece","yearTau","strataTau","strataYearTau","vesselYearTau","C.pos","C.bin")
    jags.data = list("y","logy3","effort","effort2","logeffort", "logeffort2","nonZeros", "n", "isNonZeroTrawl", "nNonZeros", "nVY", "nSY", "nS", "nY", "year", "vesselYear", "strata", "strataYear","ones.vec","R","X.pos","X.bin")
  
    if(Parallel==TRUE) {
    	modelFit= jags.parallel(jags.data, inits = NULL, parameters.to.save= jags.params, model.file=model.name, n.chains = mcmc.control$chains, n.burnin = mcmc.control$burn, n.thin = mcmc.control$thin, n.iter = as.numeric(mcmc.control$burn+mcmc.control$iterToSave), DIC = TRUE)
    }
    if(Parallel==FALSE) {
  	  modelFit = jags(jags.data, inits = NULL, parameters.to.save= jags.params, model.file=model.name, n.chains = mcmc.control$chains, n.burnin = mcmc.control$burn, n.thin = mcmc.control$thin, n.iter = as.numeric(mcmc.control$burn+mcmc.control$iterToSave), DIC = TRUE)
    }
  }
  
  ######################################################################
  # Create a list object that contains all of the properties of this run
  ######################################################################
  # modelStructure (list)
  # model.name (string)
  # fit.model (boolean)
  # mcmc.control (list)
  # Parallel (boolean)
  ######################################################################
  functionCall = list("modelStructure"=modelStructure, "model.name"=model.name, "fit.model"=fit.model, "mcmc.control"=mcmc.control, "Parallel"=Parallel,"likelihood"=likelihood,"covariates"=covariates) # species, likelihood, parameters TRUE/FALSE
  
  ######################################################################
  # Return a list of which parameters are actually estimated
  ######################################################################
  estimatedParameters = NA
  if(fit.model) {
    estimatedParameters = data.frame("Parameter" = colnames(modelFit$BUGSoutput$sims.matrix))
    estimatedParameters$Estimated = rep(TRUE, dim(modelFit$BUGSoutput$sims.matrix)[2])
    vars = apply(modelFit$BUGSoutput$sims.matrix,2,var) # figure out which cols have var = 0
    estimatedParameters$Estimated[which(vars==0)] = FALSE
  }
  
  return(c(modelFit,functionCall,estimatedParameters,Species=Species,Data=Data))

}

############################################################################
# FUNCTION logDensity WRITTEN BY ERIC WARD, UPDATED 9/30/2012. 
# EMAIL: ERIC.WARD@NOAA.GOV
############################################################################
logDensity = function(obj) {
	#################################
	# To do: ECEs
	################################# 
	attach.jags(obj) # attach object to workspace
	
    # Divide data into a zero and non-zero component 
    # For each data point (i = 1:N), calculate the appropriate
    # likelihood, and avg over all mcmc draws. Then multiply these 
    # values OR sum in log-space over data points.
     
	strata = Data$strata
	year = Data$year
	strataYear = Data$strataYear
	vesselYear = Data$vesselYear
	effort = Data$effort
	effort2 = Data$effort2
	logeffort = Data$logeffort
	logeffort2 = Data$logeffort2
	
    n.mcmc = dim(Sdev)[1]
	# calculate the predicted components for the positive model
    avg.pos = 0
	# calculate the predicted components for the positive model
	for(j in 1:length(nonZeros)) { # loop over records
	      # for the jth data point, over all mcmc draws calculate expected value
	      exp.value = Sdev[,as.numeric(strata[nonZeros[j]])] + Ydev[,as.numeric(year[nonZeros[j]])] + SYdev[,as.numeric(strataYear[nonZeros[j]])] + VYdev[,as.numeric(vesselYear[nonZeros[j]])] + logeffort[nonZeros[j]]*B.pos[,1] + logeffort2[nonZeros[j]]*B.pos[,2]
	     if(obj$likelihood=="gamma") {
	   	   gamma.a = (1/CV[,1])^2 # gamma.a is constant over data pts
	   	   gamma.b = gamma.a/exp(exp.value) # b = a/E[x]
	   	   avg.pos[j] = log(mean(dgamma(y[nonZeros[j]], shape = gamma.a, rate=gamma.b)))
	     }
	     if(obj$likelihood=="lognormal") {
	   	   sigma <- sqrt(log(CV[,1]^2+1))
	   	   avg.pos[j] = log(mean(dlnorm(y[nonZeros[j]], meanlog = exp.value, sdlog = sigma)))
	     }
	     if(obj$likelihood=="invGaussian") {
	   	   oneOverCV2 = (1/CV[,1])^2
	   	   u = exp(exp.value)
           lambda = u/oneOverCV2
           avg.pos[j] = log(mean(dinvgauss(rep(y[nonZeros[j]],n.mcmc), u, lambda)))
	     }
	 }
	 
     # do the same thing with the binomial model
     avg.bin = 0
     for(j in 1:dim(Data)[1]) {
     	# for this data point calculate posterior distribution of expectation
        exp.value = plogis(pSdev[,as.numeric(strata[j])] + pYdev[,as.numeric(year[j])] + pSYdev[,as.numeric(strataYear[j])] + pVYdev[,as.numeric(vesselYear[j])] + effort[j]*B.zero[,1] + effort2[j]*B.zero[,2])
        # calculate log of the integral over parameters, P(y|y,theta)
	    avg.bin[j] = log(mean(dbinom(rep(Data$isNonZeroTrawl[j],length=dim(Data)[1]),prob=exp.value,size=1)))    	
     }
     # return(sum(avg.bin)+sum(avg.pos))
     return(list('Presence.score'=avg.bin, 'Positive.score'=avg.pos))
}


############################################################################
# FUNCTIONS BELOW THIS POINT WRITTEN BY JIM THORSON, UPDATED 9/30/2012. 
# EMAIL: JAMES.THORSON@NOAA.GOV
#  (Some was based on code from John Wallace.)
############################################################################

######################################
#
# Useful functions
#
######################################
strata.fn <- function(x,Strata.df) {
	# function per A. Hicks, 5/5/2012
  tmpL <- as.numeric(x["BEST_LAT_DD"])>Strata.df$SLat & as.numeric(x["BEST_LAT_DD"])<=Strata.df$NLat
  tmpD <- as.numeric(x["BEST_DEPTH_M"])>Strata.df$MinDepth & as.numeric(x["BEST_DEPTH_M"])<=Strata.df$MaxDepth
  Char = as.character(Strata.df[tmpL&tmpD,"STRATA"]) 
  return(ifelse(length(Char)==0,NA,Char))
}

# -------------------------------------
LoadFn <- function (file, ...) {
    ls.ext <- function(file) {
      local({
         base::load(file)
         base::ls()
      }) }

    base::load(file, .GlobalEnv, ...)
    ls.ext(file)
  }
# -------------------------------------

cMx=function(Input){as.matrix(Input)}

# -------------------------------------

readIn <- function(ncol, nlines, ...) {
  x <- matrix(scan(,"", quiet=TRUE, nlines=nlines), ncol=ncol, byrow=TRUE)
  x <- as.data.frame(x, ..., stringsAsFactors = FALSE)
  names(x) <- x[1,]
  x <- x[-1,]
  rownames(x) <- 1:nrow(x)
  for( i in 1:ncol(x)){
    x[,i] <- if(all(is.na(suppressWarnings(as.numeric(x[,i]))))){ x[,i] }else{ as.numeric(x[,i]) }
  }
  return(x)
}

####################
# Function to calculate the correlations from the precision matrix
####################
corFunction = function(McmcArray, Parameter,this.names, Model) {
	nchains = dim(McmcArray)[2]
	cors = matrix(NA, nrow = dim(McmcArray)[1], ncol = nchains)
	# Shortcut for calculating correlation of 2x2
	# Matrix A = [a b c d]
	# Ainv = 1/(ad -bc)*[d -b -c a], cor = -b/sqrt(d*a)
    for(i in 1:nchains) {
      thisMat = McmcArray[,i,grep(this.names,Model$Parameter)]
      cors[,i] = -thisMat[,3]/sqrt(thisMat[,1]*thisMat[,4])
    }
    return(cors)
}

######################################
#
# Data visualization
#
######################################

##########
# Plot data
##########
PlotData = function(Data, FileName, Folder=NA){
  
  if(is.na(Folder)) Folder = getwd()
  
  Data = cbind(Data, 'Pres'=ifelse(Data[,'HAUL_WT_KG']>0,1,0))
  Pos = Data[which(Data[,'HAUL_WT_KG']>0),]
  Nyears = length(unique(Data[,'PROJECT_CYCLE']))
  
  # Histogram of positive catch | year
  Ncol = ceiling(sqrt(Nyears)); Nrow = ceiling(Nyears/Ncol)
  jpeg(paste(Folder,FileName,"Positive catch BY Year.jpg",sep=""),width=Ncol*3,height=Nrow*3,units="in",res=200)
    par(mfrow=c(Nrow,Ncol), mar=c(2,2,2,0), mgp=c(1.25,0.25,0), tck=-0.02, oma=c(4,4,0,0))
    for(YearI in 1:Nyears){
      Which = which(Pos[,'PROJECT_CYCLE']==unique(Pos[,'PROJECT_CYCLE'])[YearI])
      hist(Pos[Which,'HAUL_WT_KG'],main=unique(Pos[,'PROJECT_CYCLE'])[YearI], xlab="", ylab="", breaks=seq(0,max(Pos[,'HAUL_WT_KG'])+1,length=100), cex.main=1.5)
    }  
    mtext("Positive catch rates",outer=TRUE,line=2,side=1,cex=2)
    mtext("Frequency",outer=TRUE,line=2,side=2,cex=2)
  dev.off()
  
  # Scatterplot of positive catch by depth
  Ncol = Nrow = 1
  jpeg(paste(Folder,FileName,"Positive catch and depth.jpg",sep=""),width=Ncol*4,height=Nrow*4,units="in",res=200)
    par(mar=c(2,2,2,0), mgp=c(1.25,0.25,0), tck=-0.02)
    Y = Pos[,'HAUL_WT_KG']
    plot(x=Pos[,'BEST_DEPTH_M'], y=ifelse(Y==0,NA,Y), log="y", main="Positive catch by depth", xlab="Depth", ylab="Positive catch rates", pch=20, col=rgb(0,0,0,alpha=0.2))
    lines(lowess(x=Pos[,'BEST_DEPTH_M'], y=ifelse(Y==0,NA,Y)), lwd=2)
  dev.off()  
  
  # Scatterplot of positive catch by depth | Year
  Ncol = ceiling(sqrt(Nyears)); Nrow = ceiling(Nyears/Ncol)
  jpeg(paste(Folder,FileName,"Positive catch and depth BY Year.jpg",sep=""),width=Ncol*3,height=Nrow*3,units="in",res=200)
    par(mfrow=c(Nrow,Ncol), mar=c(2,2,2,0), mgp=c(1.25,0.25,0), tck=-0.02, oma=c(4,4,0,0))
    for(YearI in 1:Nyears){
      Which = which(Pos[,'PROJECT_CYCLE']==unique(Pos[,'PROJECT_CYCLE'])[YearI])
      Y = Pos[Which,'HAUL_WT_KG']
      plot(x=Pos[Which,'BEST_DEPTH_M'], y=ifelse(Y==0,NA,Y), log="y", main=unique(Pos[,'PROJECT_CYCLE'])[YearI], xlab="", ylab="", pch=20, col=rgb(0,0,0,alpha=0.2))
      lines(lowess(x=Pos[Which,'BEST_DEPTH_M'], y=ifelse(Y==0,NA,Y)), lwd=2)
    }  
    mtext("Depth",outer=TRUE,line=2,side=1,cex=2)
    mtext("Positive catch rates",outer=TRUE,line=2,side=2,cex=2)
  dev.off()
  # Line 405 from John Wallace's "Survey.Biomass.GlmmBUGS.ver.3.00.R"
  #xyplot(ifelse(HAUL_WT_KG==0,min(HAUL_WT_KG[HAUL_WT_KG!=0])/2,HAUL_WT_KG) ~ BEST_DEPTH_M | factor(PROJECT_CYCLE), data=Data, ylab = "Log of Weight (kg)", xlab="Depth (m)")
  
  # Scatterplot of positive catch by lattitude
  Ncol = Nrow = 1
  jpeg(paste(Folder,FileName,"Positive catch and latitude.jpg",sep=""),width=Ncol*4,height=Nrow*4,units="in",res=200)
    par(mar=c(2,2,2,0), mgp=c(1.25,0.25,0), tck=-0.02)
    Y = Pos[,'HAUL_WT_KG']
    plot(x=Pos[,'BEST_LAT_DD'], y=ifelse(Y==0,NA,Y), log="y", main="Positive catch by Latitude", xlab="Latitude", ylab="Positive catch rates", pch=20, col=rgb(0,0,0,alpha=0.2))
    lines(lowess(x=Pos[,'BEST_LAT_DD'], y=ifelse(Y==0,NA,Y)), lwd=2)
  dev.off()  

  # Scatterplot of positive catch by lattitude | Year
  Ncol = ceiling(sqrt(Nyears)); Nrow = ceiling(Nyears/Ncol)
  jpeg(paste(Folder,FileName,"Positive catch and latitude BY Year.jpg",sep=""),width=Ncol*3,height=Nrow*3,units="in",res=200)
    par(mfrow=c(Nrow,Ncol), mar=c(2,2,2,0), mgp=c(1.25,0.25,0), tck=-0.02, oma=c(4,4,0,0))
    for(YearI in 1:Nyears){
      Which = which(Pos[,'PROJECT_CYCLE']==unique(Pos[,'PROJECT_CYCLE'])[YearI])
      Y = Pos[Which,'HAUL_WT_KG']
      plot(x=Pos[Which,'BEST_LAT_DD'], y=ifelse(Y==0,NA,Y), log="y", main=unique(Pos[,'PROJECT_CYCLE'])[YearI], xlab="", ylab="", pch=20, col=rgb(0,0,0,alpha=0.2))
      lines(lowess(x=Pos[Which,'BEST_LAT_DD'], y=ifelse(Y==0,NA,Y)), lwd=2)
    }  
    mtext("Latitude",outer=TRUE,line=2,side=1,cex=2)
    mtext("Positive catch rates",outer=TRUE,line=2,side=2,cex=2)
  dev.off()
  
  # Bar graph of proportion positive by year
  Ncol = Nrow = 1
  jpeg(paste(Folder,FileName,"Presence and year.jpg",sep=""),width=Ncol*4,height=Nrow*4,units="in",res=200)
    par(mar=c(2,2,2,0), mgp=c(1.25,0.25,0), tck=-0.02)
    SumPres = tapply(Data[,'Pres'],INDEX=Data[,'PROJECT_CYCLE'],FUN=sum,na.rm=TRUE)
    SumAbs = tapply(1-Data[,'Pres'],INDEX=Data[,'PROJECT_CYCLE'],FUN=sum,na.rm=TRUE)
    Prop = tapply(Data[,'Pres'],INDEX=Data[,'PROJECT_CYCLE'],FUN=mean,na.rm=TRUE)
    #plot(x=unique(Data[,'PROJECT_CYCLE']), y=Prop, main="Presence/absence by year", xlab="Year", ylab="Proportion positive", pch=20, type="l", ylim=c(0,1))
    barplot(Prop, , main="Presence/absence by year", xlab="Year", ylab="Proportion positive", ylim=c(0,1))
    #barplot(-SumAbs, ylim=c(0,max(SumPres+SumAbs)),add=TRUE)
  dev.off()  
  
  # Bar graph of proportion positive by 25 meter depth bins | Year
  Ncol = ceiling(sqrt(Nyears)); Nrow = ceiling(Nyears/Ncol)
  BinWidth = 50
  jpeg(paste(Folder,FileName,"Presence and depth BY year.jpg",sep=""),width=Ncol*3,height=Nrow*3,units="in",res=200)
    par(mfrow=c(Nrow,Ncol), mar=c(2,2,2,0), mgp=c(1.25,0.25,0), tck=-0.02, oma=c(4,4,0,0))
    for(YearI in 1:Nyears){
      Which = which(Data[,'PROJECT_CYCLE']==unique(Data[,'PROJECT_CYCLE'])[YearI])
      X = unique(BinWidth*floor(Data[Which,'BEST_DEPTH_M']/BinWidth))
      Order = order(X)
      Prop = tapply(Data[Which,'Pres'],INDEX=BinWidth*floor(Data[Which,'BEST_DEPTH_M']/BinWidth),FUN=mean,na.rm=TRUE)[Order]
      #plot(x=X[Order], y=Prop, main=unique(Data[,'PROJECT_CYCLE'])[YearI], xlab="", ylab="", pch=20, type="l", ylim=c(0,1))
      barplot(Prop, main=unique(Data[,'PROJECT_CYCLE'])[YearI], xlab="", ylab="", ylim=c(0,1))
    }
    mtext("Depth bin",outer=TRUE,line=2,side=1,cex=2)
    mtext("Proportion positive",outer=TRUE,line=2,side=2,cex=2)
  dev.off()  

  # Bar graph of proportion positive by 1 degree latitude bins | Year
  Ncol = ceiling(sqrt(Nyears)); Nrow = ceiling(Nyears/Ncol)
  BinWidth = 1
  jpeg(paste(Folder,FileName,"Presence and latitude BY year.jpg",sep=""),width=Ncol*3,height=Nrow*3,units="in",res=200)
    par(mfrow=c(Nrow,Ncol), mar=c(2,2,2,0), mgp=c(1.25,0.25,0), tck=-0.02, oma=c(4,4,0,0))
    for(YearI in 1:Nyears){
      Which = which(Data[,'PROJECT_CYCLE']==unique(Data[,'PROJECT_CYCLE'])[YearI])
      X = unique(BinWidth*floor(Data[Which,'BEST_LAT_DD']/BinWidth))
      Order = order(X)
      Prop = tapply(Data[Which,'Pres'],INDEX=BinWidth*floor(Data[Which,'BEST_LAT_DD']/BinWidth),FUN=mean,na.rm=TRUE)[Order]
      #plot(x=X[Order], y=Prop, main=unique(Data[,'PROJECT_CYCLE'])[YearI], xlab="", ylab="", pch=20, type="l", ylim=c(0,1))
      barplot(Prop, main=unique(Data[,'PROJECT_CYCLE'])[YearI], xlab="", ylab="", ylim=c(0,1))
    }
    mtext("Latitude bin",outer=TRUE,line=2,side=1,cex=2)
    mtext("Proportion positive",outer=TRUE,line=2,side=2,cex=2)
  dev.off()  
  
}

##########
# Plot tows on a map
##########
MapData = function(Data, SA3, strata.limits, FileName, Folder){

  # Distilled from line 426-539 of "Survey.Biomass.GlmmBUGS.ver.3.00.R" from John Wallace's code
  jpeg(file=paste(Folder,FileName,"TowMap.jpg",sep=""),width=8,height=8,res=200,units="in")

    # Draw box
    plot(c(-55, -1280), c(32, 50.5), xlab = "Depth (m)", ylab = "Latitude", xlim=c(-1280, -55), ylim=c(32, 49), type = "n")
      abline(v= -unique(c(SA3$MIN_DEPTH_M, SA3$MAX_DEPTH_M)), h=unique(c(SA3$MIN_LAT_DD, SA3$MAX_LAT_DD)), col="grey78")
    abline(h=34.5, v=-c(30, 100, 300, 700)*1.8288, col='red')
    
    # Draw centroid of bins
    avelat <- apply(cbind(SA3$MIN_LAT_DD, SA3$MAX_LAT_DD), 1, mean)
    avedep <- apply(cbind(SA3$MIN_DEPTH_M, SA3$MAX_DEPTH_M), 1, mean)
    points(-avedep, avelat, cex=0.5)
    
    # Label areas
    text(-1235, 34.7, "Starting in 2004 NWFSC survey sampling density changes at Pt. Conception (34.5)", adj=0, col='red')
    text(-1235, 49.25, "INPFC Areas", adj=0, col='blue')
    abline(h=c(32, 36, 40.5, 43, 47.5), col='blue')
    text(-1235, 33.25, "Conception", adj=0, col='blue')
    text(-1235, 38.25, "Monterey", adj=0, col='blue')
    text(-1235, 41.75, "Eureka", adj=0, col='blue')
    text(-1235, 44.5, "Columbia", adj=0, col='blue')
    text(-1235, 48.25, "Vancouver", adj=0, col='blue')
    
    # Plot strata
    S = strata.limits
    for (i in 1:nrow(S)) {
      polygon(-c(S$MinDepth[i], S$MaxDepth[i] , S$MaxDepth[i], S$MinDepth[i]), c(S$SLat[i], S$SLat[i], S$NLat[i], S$NLat[i]), col = rainbow(nrow(S), alpha=0.3)[i])
      text(-mean(c(S$MinDepth[i], S$MaxDepth[i])), mean(c(S$SLat[i], S$NLat[i])), S$STRATA[i], cex=1.2)
    }
    
    # Plot Absence tows
    points(-Data$BEST_DEPTH_M[Data$HAUL_WT_KG==0], Data$BEST_LAT_DD[Data$HAUL_WT_KG==0], pch=16, cex=0.5, col=rgb(red=1,0,0,alpha=0.2))
    
    # Plot presence tows by year
    DataPos <- Data[Data$HAUL_WT_KG > 0 & !is.na(Data$HAUL_WT_KG),] # Temp Sp.pos - redefined below
    CU <- c("black", "green", "blue", "cyan", "purple", "grey", "orange", "hotpink", "brown", "darkolivegreen2" ,"darkslategrey",   "deepskyblue1"  , "violet", "cyan", "magenta", "lightsalmon", "gold")
    for(i in 1:length(unique(DataPos$PROJECT_CYCLE))) {
      Which = which(DataPos$PROJECT_CYCLE==unique(DataPos$PROJECT_CYCLE)[i])
      points(-DataPos$BEST_DEPTH_M[Which], DataPos$BEST_LAT_DD[Which], pch=16, cex=0.5, col= CU[i])
    }
  dev.off()
}


######################################
#
# Convergence diagnostics
#
######################################

##########
# Trace plots for each parameters
##########
ConvergencePlot = function(McmcArray, maxDims=8, parToMonitor, parnames, Nkeep=5000, FileName, Type="Trace", Folder=NA, Model){
  
  if(is.na(Folder)) Folder = paste(getwd(),"/",sep="")  
  
  Nparam = length(parToMonitor)
  Nplots = ceiling( Nparam / maxDims^2 )
  KeepSet = seq(1,dim(McmcArray)[1],length=min(dim(McmcArray)[1],Nkeep))

  if(Type%in%c("Trace","ACF")){
    for(PlotI in 1:Nplots){
      ParSet = (PlotI-1)*maxDims^2 + 1:min(Nparam-(PlotI-1)*maxDims^2,maxDims^2)
        Ncol=ceiling(sqrt(length(ParSet)))
        Nrow = ceiling(length(ParSet)/Ncol)
      jpeg(paste(Folder,FileName,"",Type,"_",PlotI,".jpg",sep=""),width=Ncol*1.5,height=Nrow*1.5,units="in",res=150)
        par(mfrow=c(Nrow,Ncol), mar=c(0,2,2,0), mgp=c(1.25,0.25,0), tck=-0.02)
        for(ParI in 1:length(ParSet)){
          if(Type=="Trace"){
            matplot(McmcArray[KeepSet,,parToMonitor[ParSet[ParI]]], type="l", lty="solid", col=rainbow(dim(McmcArray)[2],alpha=0.4), main=parnames[ParSet[ParI]], xaxt="n", xlab="",ylab="", lwd=2)
          }
          if(Type=="ACF"){
            Acf = apply(McmcArray[KeepSet,,parToMonitor[ParSet[ParI]],drop=FALSE],MARGIN=2,FUN=function(Vec){acf(Vec,plot=FALSE)$acf})
            if(!any(is.na(Acf))){ 
              matplot(Acf, type="h", lty="solid", col=rainbow(dim(McmcArray)[2],alpha=0.4), main=parnames[ParSet[ParI]], xaxt="n", xlab="",ylab="", ylim=c(0,1), lwd=2)
            }else{
              plot.new()
            }
          } 
        } # ParI loop
      dev.off()
    } # PlotI loop
  } 
  if(Type=="VarianceDensity"){
    ParSet = grep("sigma",parnames)
    Ncol=ceiling(sqrt(length(ParSet)))
      Nrow = ceiling(length(ParSet)/Ncol)
    jpeg(paste(Folder,FileName,"Variance_density.jpg",sep=""),width=Ncol*3,height=Nrow*3,units="in",res=200)
      par(mfrow=c(Nrow,Ncol), mar=c(2,2,2,0), mgp=c(1.25,0.25,0), tck=-0.02)
      for(ParI in 1:length(ParSet)){
        plot(x=-999,y=-999,,main=parnames[ParSet[ParI]],xlab="",ylab="",xlim=range(McmcArray[,,parToMonitor[ParSet[ParI]]]),ylim=c(0,10/diff(range(McmcArray[,,parToMonitor[ParSet[ParI]]])))) 
        for(ChainI in 1:dim(McmcArray)[2]){
          lines(density(McmcArray[,ChainI,parToMonitor[ParSet[ParI]]]),col=rainbow(dim(McmcArray)[2],alpha=0.4)[ChainI])
        }
      }
    dev.off()
  }
  
  if(Type=="CorrelationDensity"){
    ParSet = grep("Tau",parnames) # divide by 4 because we're only plotting correlation
    Ncol=ceiling(sqrt(length(ParSet)/4))
    Nrow = ceiling((length(ParSet)/4)/Ncol)
    corMcmc = array(NA, dim = c(dim(McmcArray)[1],dim(McmcArray)[2], length(ParSet)/4))
    corNames = c("strataYearTau","vesselYearTau","strataTau","yearTau")
    kept = 0
    keptNames = ""
    for(i in 1:4) {
    	if(length(grep(corNames[i],parnames)) > 0) {
    		kept = kept + 1
    		corMcmc[,,kept] = corFunction(McmcArray, Model$Parameter, corNames[i], Model)
    		keptNames[kept] = corNames[i]
    	}    	
    }
    plotNames = c("Strata-Year","Vessel-Year","strata","year")
    jpeg(paste(Folder,FileName,"Correlation_density.jpg",sep=""),width=Ncol*3,height=Nrow*3,units="in",res=200)
      par(mfrow=c(Nrow,Ncol), mar=c(2,2,2,0), mgp=c(1.25,0.25,0), tck=-0.02)
    # create a new subdataset from just the correlations
    for(par in 1:kept) {
    	dens = density(corMcmc[,1,par],from=-1,to=1) # this is for scaling ylim
    	plot(corMcmc[1,1,1], xlim=c(-1,1), ylim=c(0, 1.1*max(dens$y/sum(dens$y))), col = "white", main = plotNames[par], ylab = "Density",xlab = "correlation")
    	for(chainI in 1:dim(McmcArray)[2]) {
    	  # for each chain, plot a separate density plot
    	  dens = density(corMcmc[,chainI,par],from=-1,to=1)
          lines(dens$x, dens$y/sum(dens$y),col=rainbow(dim(McmcArray)[2],alpha=0.4)[chainI])    		
    	}
    }
    dev.off()
  }
  
}

##########
# Gelman-Rubin diagnostics
##########
McmcDiagnosticsPlot = function(McmcList, parToMonitor, FileName, Folder=NA, Geweke=FALSE){

  GelmanDiag = gelman.diag(McmcList[,parToMonitor])

  if(is.na(Folder)) Folder = paste(getwd(),"/",sep="")  

  # Plot
  jpeg(paste(Folder,FileName,"Gelman_points.jpg",sep=""),width=6,height=6,units="in",res=200)
    par(mfrow=c(1,1))
    # all points below 1.05 shown as points, all above given names
    mycols = rep("black",length(GelmanDiag$psrf[,1]))
    mycols[which(GelmanDiag$psrf[,1] > 1.05)] = "white"
    plot(GelmanDiag$psrf[,1],xlab="psrf",ylab="Gelman: point estimate",main="",col=mycols,lwd=2)
    lines(c(-1000,1000),c(1.05,1.05),lwd=2,col="red")
    if(length(which(GelmanDiag$psrf[,1] > 1.05)) > 0) {
    	text(seq(1,length(GelmanDiag$psrf[,1]))[which(GelmanDiag$psrf[,1] > 1.05)],y = GelmanDiag$psrf[which(GelmanDiag$psrf[,1] > 1.05),1],labels=names(GelmanDiag$psrf[,1])[which(GelmanDiag$psrf[,1] > 1.05)])
    }
  dev.off()
  
  # Plot histograms showing Gelman-Rubin scores
  jpeg(paste(Folder,FileName,"Gelman_histograms.jpg",sep=""),width=5,height=10,units="in",res=200)
    par(mfrow=c(2,1),mai=c(0.5,0.5,0.2,0.2))
    hist(GelmanDiag$psrf[,1],xlab="psrf",main="Gelman: point estimate")
    hist(GelmanDiag$psrf[,2],xlab="psrf",main="Gelman: upper C.I.")
  dev.off()

  # plot the geweke diagnostic   -- OFTEN DOESN"T CONVERGE
  if(Geweke==TRUE){
    GewekeDiag = geweke.diag(McmcList[,parToMonitor])
    jpeg(paste(Folder,FileName,"Geweke.jpg",sep=""),width=8,height=8,units="in",res=200)
      par(mfrow=c(3,1),mai=c(0.5,0.5,0.2,0.2))
      for(i in 1:3) {
      	plot(unlist(GewekeDiag[[i]]), xlab=paste("Parameter, chain: ",i),ylab="Geweke Z-score",ylim=c(-3,3))
      	lines(c(-1000,1000),c(1.96,1.96),col="red",lty=3)
      	lines(c(-1000,1000),c(-1.96,-1.96),col="red",lty=3)	
      }
    dev.off()
  }
}

#######################################
#
# Model fit diagnostics
#
#######################################

##########
# Display offset
##########
PlotOffset = function(Data, BugsList, maxDims=8, FileName, Folder=NA){
  if(is.na(Folder)) Folder = paste(getwd(),"/",sep="")  
  
  # Attach stuff
  attach(BugsList)
  #attach(Data)
  nonZeros = which(isNonZeroTrawl==TRUE)
    
  # Positive offset
  LogEffortRange = seq(min(logeffort),max(logeffort),length=1000)
  Nstrat = length(unique(strataYear[nonZeros]))
  Nplots = ceiling( Nstrat / maxDims^2 )
  for(PlotI in 1:Nplots){
    ParSet = (PlotI-1)*maxDims^2 + 1:min(Nstrat-(PlotI-1)*maxDims^2,maxDims^2)
    Ncol=ceiling(sqrt(length(ParSet)))
    Nrow = ceiling(length(ParSet)/Ncol)
    jpeg(paste(Folder,"/",FileName,"Offset_Positive_by_strata_",PlotI,".jpg",sep=""),width=Ncol*1.5,height=Nrow*1.5,res=150,units="in")
      par(mfrow=c(Nrow,Ncol), mar=c(2,2,2,0), mgp=c(1.25,0.25,0), tck=-0.02)
      for(StrataYearI in ParSet){
        StrataI = strata[which(strataYear==unique(strataYear)[StrataYearI])[1]]
        YearI = year[which(strataYear==unique(strataYear)[StrataYearI])[1]]
        Y = exp( median(cMx(Sdev)[,StrataI]) + median(Ydev[,YearI]) + median(SYdev[,StrataYearI]) + median(B.pos[,1])*LogEffortRange + median(B.pos[,2])*LogEffortRange )
        plot(x=exp(LogEffortRange),y=Y,ylab="",xlab="",main=unique(strataYear[nonZeros])[StrataI],ylim=c(0,max(Y)),xlim=c(0,max(exp(LogEffortRange))), type="l")
      }
    dev.off()
  }

  # Presence/Absence offset
  LogEffortRange = seq(min(logeffort),max(logeffort),length=1000)
  Nstrat = length(unique(strataYear[nonZeros]))
  Nplots = ceiling( Nstrat / maxDims^2 )
  for(PlotI in 1:Nplots){
    ParSet = (PlotI-1)*maxDims^2 + 1:min(Nstrat-(PlotI-1)*maxDims^2,maxDims^2)
    Ncol=ceiling(sqrt(length(ParSet)))
    Nrow = ceiling(length(ParSet)/Ncol)
    jpeg(paste(Folder,"/",FileName,"Offset_Presence-Absence_by_strata_",PlotI,".jpg",sep=""),width=Ncol*1.5,height=Nrow*1.5,res=150,units="in")
      par(mfrow=c(Nrow,Ncol), mar=c(2,2,2,0), mgp=c(1.25,0.25,0), tck=-0.02)
      for(StrataYearI in ParSet){
        StrataI = strata[which(strataYear==unique(strataYear)[StrataYearI])[1]]
        YearI = year[which(strataYear==unique(strataYear)[StrataYearI])[1]]
        Y = plogis( median(cMx(pSdev)[,StrataI]) + median(pYdev[,YearI]) + median(pSYdev[,StrataYearI]) + median(B.zero[,1])*LogEffortRange + median(B.zero[,2])*LogEffortRange )
        plot(x=exp(LogEffortRange),y=Y,ylab="",xlab="",main=unique(strataYear[nonZeros])[StrataI],ylim=c(0,1),xlim=c(0,max(exp(LogEffortRange))), type="l")
      }
    dev.off()
  }

  # Detach stuff
  #detach(Data)
  detach(BugsList)
}

##########
# Posterior predictive distribution
##########

PosteriorPredictive = function(Data, Model, maxDims=6, FileName, Folder=NA){
  if(is.na(Folder)) Folder = paste(getwd(),"/",sep="")  
  
  # Attach stuff
  attach(Model$BUGSoutput$sims.list)
  #attach(Data)
  Dist = Model$likelihood
  nonZeros = which(Data[,'isNonZeroTrawl']==TRUE)
  u.nz = matrix(NA, nrow=nrow(B.pos), ncol=length(nonZeros))
  
  # Warnings about mismatch between data and model
  if(nlevels(strata)!=ncol(cMx(Sdev)) | nlevels(year)!=ncol(Ydev) | nlevels(strataYear)!=ncol(SYdev)| nlevels(vesselYear)!=ncol(VYdev)){
    stop("Model and data do not match")
  }
  
  if(Dist=="lognormal"){
    sigma = sqrt(log(1+(CV[,1]^2)))
    for(i in 1:length(nonZeros)){
      u.nz[,i] <- exp( cMx(Sdev)[,strata[nonZeros[i]]] + Ydev[,year[nonZeros[i]]] + VYdev[,vesselYear[nonZeros[i]]] + SYdev[,strataYear[nonZeros[i]]] + B.pos[,1]*logeffort[nonZeros[i]] + B.pos[,2]*logeffort2[nonZeros[i]] )
    }
  }
  if(Dist=="gamma"){
    gamma.a = oneOverCV2 = 1/(CV[,1]^2)
    for(i in 1:length(nonZeros)){
      Temp = cMx(Sdev)[,strata[nonZeros[i]]] + Ydev[,year[nonZeros[i]]] + VYdev[,vesselYear[nonZeros[i]]] + SYdev[,strataYear[nonZeros[i]]] + B.pos[,1]*logeffort[nonZeros[i]] + B.pos[,2]*logeffort2[nonZeros[i]]
      u.nz[,i] <- exp( ifelse(Temp>100,100,Temp) ) # Eric included this in the JAGS code
    }
  }
  if(Dist=="invGaussian"){
    oneOverCV2 = 1/(CV[,1]^2)
    for(i in 1:length(nonZeros)){
      Temp = cMx(Sdev)[,strata[nonZeros[i]]] + Ydev[,year[nonZeros[i]]] + VYdev[,vesselYear[nonZeros[i]]] + SYdev[,strataYear[nonZeros[i]]] + B.pos[,1]*logeffort[nonZeros[i]] + B.pos[,2]*logeffort2[nonZeros[i]]
      u.nz[,i] <- exp( ifelse(Temp>100,100,Temp) ) # Eric included this in the JAGS code
    }
  }
  if(Dist=="lognormalECE"){
    u.nz2 = array(NA, dim=dim(u.nz))
    sigma = sqrt(log(1+(CV[,1]^2)))
    sigma2 = sqrt(log(1+(CV[,2]^2)))
    for(i in 1:length(nonZeros)){
      Temp = cMx(Sdev)[,strata[nonZeros[i]]] + Ydev[,year[nonZeros[i]]] + VYdev[,vesselYear[nonZeros[i]]] + SYdev[,strataYear[nonZeros[i]]] + B.pos[,1]*logeffort[nonZeros[i]] + B.pos[,2]*logeffort2[nonZeros[i]]
      u.nz[,i] <- exp( ifelse(Temp>100,100,Temp) ) # Eric included this in the JAGS code
    }
    for(i in 1:length(nonZeros)){
      u.nz2[,i] <- u.nz[,i] * ratio
    }
  }
  if(Dist=="gammaECE"){
    u.nz2 = array(NA, dim=dim(u.nz))
    gamma.a = 1/(CV[,1]^2)
    gamma.a2 = 1/(CV[,2]^2)    
    for(i in 1:length(nonZeros)){
      Temp = cMx(Sdev)[,strata[nonZeros[i]]] + Ydev[,year[nonZeros[i]]] + VYdev[,vesselYear[nonZeros[i]]] + SYdev[,strataYear[nonZeros[i]]] + B.pos[,1]*logeffort[nonZeros[i]] + B.pos[,2]*logeffort2[nonZeros[i]]
      u.nz[,i] <- exp( ifelse(Temp>100,100,Temp) ) # Eric included this in the JAGS code
    }
    for(i in 1:length(nonZeros)){
      u.nz2[,i] <- u.nz[,i] * ratio
    }
  }
  
  Q = rep(NA, nrow(Data)) # vector to track quantiles for each observation
  Nstrat = length(unique(strataYear[nonZeros]))
  Ncol=ceiling(sqrt(Nstrat)); Nrow = ceiling(Nstrat/Ncol)
  Nstrat = length(unique(strataYear[nonZeros]))
  Nplots = ceiling( Nstrat / maxDims^2 )
  for(PlotI in 1:Nplots){
    ParSet = (PlotI-1)*maxDims^2 + 1:min(Nstrat-(PlotI-1)*maxDims^2,maxDims^2)
    Ncol=ceiling(sqrt(length(ParSet)))
    Nrow = ceiling(length(ParSet)/Ncol)
    jpeg(paste(Folder,"/",FileName,"Posterior_Predictive_",PlotI,".jpg",sep=""),width=Ncol*1.5,height=Nrow*1.5,res=200,units="in")
      par(mfrow=c(Nrow,Ncol), mar=c(2,2,2,0), mgp=c(1.25,0.25,0), tck=-0.02)
      for(StrataYearI in ParSet){
        Which = which(strataYear[nonZeros]==unique(strataYear[nonZeros])[StrataYearI])
        plot(Data[nonZeros[Which],'HAUL_WT_KG'],ylab="",xlab="",log="y",main=unique(strataYear[nonZeros])[StrataYearI],col="blue", ylim=range(Data[nonZeros,'HAUL_WT_KG']))
        # mean(u.nz[,2])
        for(ObsI in 1:length(Which)){
          if(Dist=="lognormal"){     
            y = rlnorm(n=1000,meanlog=log(u.nz[,Which[ObsI]]),sdlog=sigma)   # Plotting in log-space
          }
          if(Dist=="gamma"){     
            b = gamma.a / u.nz[,Which[ObsI]];    
            y = rgamma(n=1000,shape=gamma.a,rate=b)
          }
          if(Dist=="invGaussian"){     
            lambda = u.nz[,Which[ObsI]]*oneOverCV2
            y = rinvgauss(n=1000,mu=u.nz[,Which[ObsI]],lambda=lambda)
          }
          if(Dist=="lognormalECE"){     
            ECE = rbinom(n=nrow(u.nz)*10, size=1, prob=p.ece[,2])
            y = rlnorm(n=1000, meanlog=log(u.nz[,Which[ObsI]])*(1-ECE)+log(u.nz2[,Which[ObsI]])*ECE, sdlog=sigma*(1-ECE)+sigma2*ECE)
          }
          if(Dist=="gammaECE"){     
            b = gamma.a / u.nz[,Which[ObsI]];    
            b2 = gamma.a2 / u.nz2[,Which[ObsI]];    
            ECE = rbinom(n=nrow(u.nz)*10, size=1, prob=p.ece[,2])
            y = rgamma(n=1000, shape=gamma.a*(1-ECE)+gamma.a2*ECE, rate=b*(1-ECE)+b2*ECE)
          }
          Q[nonZeros[Which[ObsI]]] = mean(y>Data[nonZeros[Which[ObsI]],'HAUL_WT_KG'])
          Quantiles = quantile(y,prob=c(0.025,0.25,0.75,0.975))
          lines(x=c(ObsI,ObsI), y=Quantiles[2:3], lwd=2)
          lines(x=c(ObsI,ObsI), y=Quantiles[c(1,4)], lwd=1,lty="dotted")
          if(Data[nonZeros[Which[ObsI]],'HAUL_WT_KG']>max(Quantiles) | Data[nonZeros[Which[ObsI]],'HAUL_WT_KG']<min(Quantiles)){
            points(x=ObsI,y=Data[nonZeros[Which[ObsI]],'HAUL_WT_KG'],pch=4,col="red",cex=2)
          }
        }
      }
    dev.off()
  }
  
  # Q-Q plot
  jpeg(paste(Folder,"/",FileName,"Q-Q_plot.jpg",sep=""),width=4,height=4,res=200,units="in")
    par(mfrow=c(1,1), mar=c(2,2,2,0), mgp=c(1.25,0.25,0), tck=-0.02)
    Qtemp = na.omit(Q)
    Order = order(Qtemp)
    plot(x=seq(0,1,length=length(nonZeros)), y=Qtemp[Order], main="Q-Q plot", xlab="Uniform", ylab="Empirical")
    abline(a=0,b=1)
  dev.off()
  
  # Detach stuff
  #detach(Data)
  detach(Model$BUGSoutput$sims.list)
}

#########
# Compute index of abundance
#########
ComputeIndices = function(Data, Model, FileName, maxDims=6, Folder=NA, Weights="StrataAreas", StrataTable){
  if(is.na(Folder)) Folder = paste(getwd(),"/",sep="")  
  
  # Attach stuff -- listed by search()
  attach(Model$BUGSoutput$sims.list)
  #attach(Data)
  modelStructure = Model$modelStructure
  Dist = Model$likelihood
  
  # Compute average AreaSwept - This AreaSwept used to calculate densities, and matters if Catch is a nonlinear function of AreaSwept
  MeanLogAreaSwept = mean(Data$logeffort)
  
  # New objects
  Chains = array(NA,dim=c(nrow(Sdev),length(unique(StrataTable[,'year'])),length(unique(StrataTable[,'strata'])),2))
  Year = Strata = Area = PosMedian = PresMedian = IndexMedian = IndexMedianWeighted = PosMean = PresMean = IndexMean = IndexMeanWeighted = CvMedian = SdLog = RawPos = RawPres = Raw = RawWeighted = matrix(NA,nrow=length(unique(StrataTable[,'year'])),ncol=length(unique(StrataTable[,'strata']))) 
  CvMedianYear = SdLogYear = rep(NA,length(unique(StrataTable[,'year'])))
  
  # Calculate values
  for(YearI in 1:nrow(PosMean)){
    for(StratI in 1:ncol(PosMean)){
      # Derived indicators
      StrataYearI = which(levels(strataYear)==paste(toupper(letters[StratI]),":",levels(year)[YearI],sep=""))
      AreaI = which(StrataTable[,'strataYear']==paste(levels(strata)[StratI],":",levels(year)[YearI],sep=""))
      Which = which(strata==toupper(letters[StratI]) & year==levels(year)[YearI])
      # Year, strata, and area
      Year[YearI,StratI] = levels(year)[YearI]
      Strata[YearI,StratI] = toupper(letters[StratI])
      if(Weights=="StrataAreas") Area[YearI,StratI] = StrataTable[AreaI,'Area_Hectares']
      if(Weights=="Equal") Area[YearI,StratI] = 1
      # Save Positive catch chain in normal-space and correct for transformation biases
      Chains[,YearI,StratI,1] = exp( cMx(Sdev)[,StratI] + cMx(Ydev)[,YearI] + cMx(SYdev)[,StrataYearI] + log(1)*B.pos[,1] + log(1)^2*B.pos[,2] ) # wardJAGS uses logeffort offset
        # Lognormal -- Bias correction 
        if(Dist=="lognormal"){
          Sigma = sqrt(log(CV[,1]^2+1))        
          Chains[,YearI,StratI,1] = Chains[,YearI,StratI,1] * exp(Sigma^2/2)
        }
        # LognormalECE -- Bias correction + incorporate ECE
        if(Dist=="lognormalECE"){
          Sigma = sqrt(log(CV[,1]^2+1))        
          Sigma2 = sqrt(log(CV[,2]^2+1))        
          Chains[,YearI,StratI,1] = Chains[,YearI,StratI,1]*exp(Sigma^2/2)*p.ece[,1] + ratio*Chains[,YearI,StratI,1]*exp(Sigma2^2/2)*p.ece[,2]
        }
        # GammaECE -- incorporate ECE 
        if(Dist=="gammaECE"){
          Chains[,YearI,StratI,1] = Chains[,YearI,StratI,1]*p.ece[,1] + ratio*Chains[,YearI,StratI,1]*p.ece[,2]
        }
        # Don't make mean-unbiased for unobserved vessel
        #if(modelStructure$VesselYear.positiveTows=="random") Chains[,YearI,StratI,1] = Chains[,YearI,StratI,1] * exp(sigmaVY[,1]^2/2)
      # Save Presence/Absence chain in normal-space and correct for transformation biases
      Chains[,YearI,StratI,2] = plogis( cMx(pSdev)[,StratI] + cMx(pYdev)[,YearI] + cMx(pSYdev)[,StrataYearI] + (1)*B.zero[,1] + (1)^2*B.zero[,2] )   # wardJAGS predicts the probability of 0 catch, and uses offset as effort, not logeffort
        # Don't make mean-unbiased for unobserved vessel (this was done incorrectly anyway)
        #if(modelStructure$VesselYear.zeroTows=="random") Chains[,YearI,StratI,2] = mean(plogis(rnorm(n=dim(Chains)[1]*1e3, mean=qlogis(Chains[,YearI,StratI,2]), sd=sigmaVY[,2])))
      # Index (median)
      PosMedian[YearI,StratI] = median( Chains[,YearI,StratI,1] )  # / 2e4 # Convert kilograms to metric tons
      PresMedian[YearI,StratI] = median( Chains[,YearI,StratI,2] )   
      IndexMedian[YearI,StratI] = median( Chains[,YearI,StratI,1] * Chains[,YearI,StratI,2] ) 
      IndexMedianWeighted[YearI,StratI] = IndexMedian[YearI,StratI] * Area[YearI,StratI]
      # Index (mean)
      PosMean[YearI,StratI] = mean( Chains[,YearI,StratI,1] )  # / 2e4 # Convert kilograms to metric tons
      PresMean[YearI,StratI] = mean( Chains[,YearI,StratI,2] )   
      IndexMean[YearI,StratI] = mean( Chains[,YearI,StratI,1] * Chains[,YearI,StratI,2] ) 
      IndexMeanWeighted[YearI,StratI] = IndexMean[YearI,StratI] * Area[YearI,StratI]
      # CV of median (from J. Wallace "Survey.Biomass.GlmmBUGS.ver.3.00.r)
      Temp = Area[YearI,StratI] * Chains[,YearI,StratI,1] * Chains[,YearI,StratI,2]
      CvMedian[YearI,StratI] = sqrt(var(Temp)) / median(Temp)
      # CV if median (from J. Wallace "Survey.Biomass.GlmmBUGS.ver.3.00.r)
      Temp = Area[YearI,StratI] * Chains[,YearI,StratI,1] * Chains[,YearI,StratI,2]
      SdLog[YearI,StratI] = sd(log(Temp))
      # Raw
      RawPos[YearI,StratI] = mean(ifelse(Data[Which,'HAUL_WT_KG']>0,Data[Which,'HAUL_WT_KG']/Data[Which,'effort'],NA),na.rm=TRUE) 
      RawPres[YearI,StratI] = mean(Data[Which,'HAUL_WT_KG']>0)    
      Raw[YearI,StratI] = RawPos[YearI,StratI] * RawPres[YearI,StratI]
      RawWeighted[YearI,StratI] = Raw[YearI,StratI] * Area[YearI,StratI]
    }
    # CV of median (from J. Wallace "Survey.Biomass.GlmmBUGS.ver.3.00.r")
    Temp = Area[YearI,StratI] * rowSums( cMx(Chains[,YearI,,1]) * cMx(Chains[,YearI,,2]) )
    CvMedianYear[YearI] = sqrt(var(Temp)) / median(Temp)
    # SD of log of index (from J. Wallace "Survey.Biomass.GlmmBUGS.ver.3.00.r")
    Temp = Area[YearI,StratI] * rowSums( cMx(Chains[,YearI,,1]) * cMx(Chains[,YearI,,2]) )
    SdLogYear[YearI] = sd(log(Temp))
  } # 1085-115

  # Plot MCMC chains for each Strata:Year
  for(Type in c("Trace","ACF")){
    for(ChainI in 1:2){
      Nstrat = length(unique(strataYear[nonZeros]))
      Nplots = ceiling( Nstrat / maxDims^2 )
      for(PlotI in 1:Nplots){
        ParSet = (PlotI-1)*maxDims^2 + 1:min(Nstrat-(PlotI-1)*maxDims^2,maxDims^2)
        Ncol=ceiling(sqrt(length(ParSet)))
        Nrow = ceiling(length(ParSet)/Ncol)
        jpeg(paste(Folder,"/",FileName,"Chain_",c("Positive","Presence")[ChainI],"_",Type,"_by_StrataYear_",PlotI,".jpg",sep=""),width=Ncol*1.5,height=Nrow*1.5,res=150,units="in")
          par(mfrow=c(Nrow,Ncol), mar=c(2,2,2,0), mgp=c(1.25,0.25,0), tck=-0.02)
          for(StrataYearI in ParSet){
            StrataI = strata[which(strataYear==unique(strataYear)[StrataYearI])[1]]
            YearI = year[which(strataYear==unique(strataYear)[StrataYearI])[1]]
            Mat = matrix(Chains[,YearI,StrataI,ChainI],ncol=Model$mcmc.control$chains,byrow=FALSE)
            if(Type=="Trace"){
              matplot(Mat, type="l", lty="solid", col=rainbow(Model$mcmc.control$chains,alpha=0.4), main=paste(YearI,StrataI), xaxt="n", xlab="",ylab="")
            }
            if(Type=="ACF"){
              Acf = apply(Mat,MARGIN=2,FUN=function(Vec){acf(Vec,plot=FALSE)$acf})
              matplot(Acf, type="h", lty="solid", col=rainbow(Model$mcmc.control$chains,alpha=0.4), main=paste(YearI,StrataI), xaxt="n", xlab="",ylab="", ylim=c(0,1), lwd=2)
            }
          }
        dev.off()
      }
    }
  }
  
  # Compile into matrices
  Results1 = data.frame(Year=as.vector(Year), Strata=as.vector(Strata), Raw=as.vector(RawWeighted), IndexMedian=as.vector(IndexMedianWeighted), IndexMean=as.vector(IndexMeanWeighted), CvMedian=as.vector(CvMedian), SdLog=as.vector(SdLog), Area=as.vector(Area), PosMedian=as.vector(PosMedian), PresMedian=as.vector(PresMedian), PosMean=as.vector(PosMean), PresMean=as.vector(PresMean), RawPos=as.vector(RawPos), RawPres=as.vector(RawPres))
  Results2 = data.frame(Year=Year[,1], Raw=rowSums(RawWeighted,na.rm=TRUE), IndexMedian=rowSums(IndexMedianWeighted,na.rm=TRUE), IndexMean=rowSums(IndexMeanWeighted,na.rm=TRUE), CvMedian=CvMedianYear, SdLog=SdLogYear)
  
  # Detach stuff -- listed by search()
  #detach(Data)
  detach(Model$BUGSoutput$sims.list)
  
  # Write and print output
  write.csv(Results1,file=paste(Folder,"/",FileName,"ResultsByYearAndStrata.csv",sep=""))
  write.csv(Results2,file=paste(Folder,"/",FileName,"ResultsByYear.csv",sep=""))

  # Return output
  Return = list(Results1=Results1, Results2=Results2)
  return(Return)
}


#########
# Compute MLE index of abundance
#########
ComputeMleIndices = function(Data, Model, FileName, Folder=NA, Weights="StrataAreas", StrataTable, Run=TRUE){
  
  # Make folder
  if(is.na(Folder)) Folder = paste(getwd(),"/",sep="")  
  
  # Attach stuff -- listed by search()
  attach(Model$BUGSoutput$sims.list)
  #attach(Data)
  modelStructure = Model$modelStructure
  Dist = Model$likelihood

  # Estimate marginal means
  Mat = matrix(NA,nrow=length(levels(year)),ncol=length(levels(strata))) 
  DataNew = data.frame(year=levels(year)[as.vector(row(Mat))], strata=levels(strata)[as.vector(col(Mat))])
  DataNew = data.frame(DataNew, logeffort=rep(log(1),nrow(DataNew)), vesselYear=rep(999,nrow(DataNew)), ones.vec=rep(log(1),nrow(DataNew)))   # I need to include ones.vec as zero for some reason

  # Only run if there's no random effects and the distribution is either Gamma and Lognormal
  if( (modelStructure$VesselYear.zeroTows%in%c("zero","fixed")) & (modelStructure$VesselYear.positiveTows%in%c("zero","fixed")) & (modelStructure$StrataYear.zeroTows%in%c("zero","fixed")) & (modelStructure$StrataYear.positiveTows%in%c("zero","fixed"))  & (Dist=="gamma" | Dist=="lognormal") & Run==TRUE){
    # Default formulae
    FormulaPres = " ~ 0 + factor(year)"
    if(nlevels(strata)>1) FormulaPres = paste(FormulaPres," + factor(strata)",sep="")
    FormulaPos = " ~ 0 + factor(year)"
    if(nlevels(strata)>1) FormulaPos = paste(FormulaPos," + factor(strata)",sep="")
  
    # Modified formulae
    if(modelStructure$StrataYear.zeroTows=="fixed" & nlevels(strata)>1){FormulaPres = paste(FormulaPres," + factor(strata):factor(year)",sep="")}  
    #if(modelStructure$VesselYear.zeroTows=="fixed"){FormulaPres = paste(FormulaPres," + factor(vesselYear)",sep="")}  
    #if(modelStructure$Catchability.zeroTows=="linear"){FormulaPres = paste(FormulaPres," + logeffort",sep="")}  
    #  if(modelStructure$Catchability.zeroTows=="quadratic"){FormulaPres = paste(FormulaPres," + logeffort + logeffort2",sep="")}  
    if(modelStructure$StrataYear.positiveTows=="fixed" & nlevels(strata)>1){FormulaPos = paste(FormulaPos," + factor(strata):factor(year)",sep="")}  
    #if(modelStructure$VesselYear.positiveTows=="fixed"){FormulaPos = paste(FormulaPos," + factor(vesselYear)",sep="")}  
    #if(modelStructure$Catchability.positiveTows=="linear"){FormulaPos = paste(FormulaPos," + logeffort",sep="")}  
    #  if(modelStructure$Catchability.positiveTows=="quadratic"){FormulaPos = paste(FormulaPos," + logeffort + logeffort2",sep="")}  
  
    #### Presence/absence
    GlmPres <- glm(as.formula(paste("ifelse(Data[,'HAUL_WT_KG']>0,1,0)",FormulaPres)), family=binomial, control=glm.control(epsilon=1e-8, maxit=1000, trace=FALSE))
    #### Positive catches
    #OffsetPos = list(ones.vec,logeffort)[[ifelse(modelStructure$Catchability.positiveTows=="one",2,1)]]
    if(Dist=="gamma") GlmPos <- glm(as.formula(paste("HAUL_WT_KG",FormulaPos)), family=Gamma(link="log"), offset=logeffort, subset=which(Data[,'HAUL_WT_KG']>0), control=glm.control(epsilon=1e-8, maxit=1000, trace=FALSE))
    if(Dist=="lognormal") GlmPos <- lm(as.formula(paste("log(HAUL_WT_KG)",FormulaPos)), offset=logeffort, subset=which(Data[,'HAUL_WT_KG']>0))
    
    # Calculate strata areas
    Area = matrix(NA,nrow=length(unique(StrataTable[,'year'])),ncol=length(unique(StrataTable[,'strata']))) 
    for(YearI in 1:nrow(Area)){
    for(StratI in 1:ncol(Area)){
      AreaI = which(StrataTable[,'strataYear']==paste(levels(strata)[StratI],":",levels(year)[YearI],sep=""))
      if(Weights=="StrataAreas") Area[YearI,StratI] = StrataTable[AreaI,'Area_Hectares']
      if(Weights=="Equal") Area[YearI,StratI] = 1
    }}
    
    # Predict indices
    IndexPres = array(predict(GlmPres, newdata=DataNew, type="response"), dim=dim(Mat))
    IndexPos = array(predict(GlmPos, newdata=DataNew, type="response"), dim=dim(Mat))
    if(Dist=="lognormal") IndexPos = exp(IndexPos + summary(GlmPos)$sigma^2/2)  # Bias correction
    Index = Area * IndexPres * IndexPos
  
  }else{
    Area = IndexPres = IndexPos = Index = matrix(NA,nrow=length(unique(StrataTable[,'year'])),ncol=length(unique(StrataTable[,'strata']))) 
  }

  # EJ's code
  #DglmData = data.frame(Catch=Data[,'HAUL_WT_KG']/Data[,'effort'], Year=Data[,'year'])
  #Dglm = dglm(DglmData, dist="gamma", J=TRUE)      # J=TRUE calculates jacknife estimates of CV
  
  # Return results
  Results1 = data.frame(year=levels(year)[as.vector(row(Mat))], strata=levels(strata)[as.vector(col(Mat))], Index=as.vector(Index), Pres=as.vector(IndexPres), Pos=as.vector(IndexPos))
  Results2 = data.frame(year=levels(year), Index=rowSums(Index), Pres=rowSums(IndexPres), Pos=rowSums(IndexPos))
 
  # Write and print output
  write.csv(Results1,file=paste(Folder,"/",FileName,"ResultsByYearAndStrata_MLE.csv",sep=""))
  write.csv(Results2,file=paste(Folder,"/",FileName,"ResultsByYear_MLE.csv",sep=""))

  # Detach stuff -- listed by search()
  #detach(Data)
  detach(Model$BUGSoutput$sims.list)

  return(list(Results1=Results1, Results2=Results2))
}    
    

########################################################
####### This block of code is related to processing output
doMCMCDiags = function(directory, mods, McmcDiagnostics=FALSE) {  

  # Identify strata and year for StratYear values
  StrataTable = data.frame( 'strataYear'=levels(strataYear), 'strata'=sapply(levels(strataYear),FUN=function(Char){strsplit(Char,":")[[1]][1]}), 'year'=sapply(levels(strataYear),FUN=function(Char){strsplit(Char,":")[[1]][2]}), 'Area_Hectares'=rep(NA,nlevels(strataYear)))
  SA3 = read.csv(paste(directory,"SA3.csv",sep=""))
  for(i in 1:nrow(StrataTable)){
    Row = which(strata.limits[,'STRATA']==StrataTable[i,'strata'])
    StrataTable[i,'Area_Hectares'] = sum(SA3[SA3[,'MAX_LAT_DD']<=strata.limits[Row,'NLat'] & SA3[,'MIN_LAT_DD']>=strata.limits[Row,'SLat'] & SA3[,'MIN_DEPTH_M']>=strata.limits[Row,'MinDepth'] & SA3[,'MAX_DEPTH_M']<=strata.limits[Row,'MaxDepth'],'AREA_HECTARES'])
  }
  
  # Make folder for plots
  Species = species = mods[[1]]$Species
  SpeciesFolder = paste(getwd(),"/",Species,"_FinalDiagnostics/",sep="")
    dir.create(SpeciesFolder, showWarnings=FALSE)
  
  # Make objects for saving output
  Indices = array(NA, dim=c(nlevels(year),length(mods),2))
  IndicesByStrata = array(NA, dim=c(nlevels(year),nlevels(strata),length(mods),2))
  
  ######################
  # Display data
  ######################
    
  # Plot data by year, depth, and latitude
  PlotData(Data=Data, FileName="", Folder=SpeciesFolder)
  # Plot location of data
  MapData(Data=Data, strata.limits=strata.limits, SA3=SA3, FileName="", Folder=SpeciesFolder)
  # Save mods for later usage
  Save = list(mods=mods, Data=Data)
  save(Save, file=paste(SpeciesFolder,"Save.RData",sep=""))
  # Load old mods (if you want to)
  #load(file=paste(SpeciesFolder,"Save.RData",sep=""))
  #attach(Save)
    
  # Format data
  ModelNumber = 1
  for(ModelNumber in 1:length(mods)){
    
    # Make folder
    Folder = paste(SpeciesFolder,"/Model=",ModelNumber,"/",sep="")
    dir.create(Folder, showWarnings=FALSE)
    
    # Unpack results
    Model = mods[[ModelNumber]]
    McmcArray = as.array(Model$BUGSoutput$sims.array)
    McmcList = vector("list",length=dim(McmcArray)[2])
    for(i in 1:length(McmcList)) McmcList[[i]] = as.mcmc(McmcArray[,i,])
    McmcList = mcmc.list(McmcList)
    BugsList = Model$BUGSoutput$sims.list
                                                                       
    # Record details
    capture.output(Model$mcmc.control, file=paste(Folder,"mcmc_control.txt",sep=""))
    capture.output(list(Model$likelihood,Model$modelStructure), file=paste(Folder,"Model_Structure.txt",sep=""))
    capture.output(Model$BUGSoutput$summary, file=paste(Folder,"BUGSoutput_summary.txt",sep=""))
    capture.output(Model$model, file=paste(Folder,"deltaGLM.txt",sep=""))
    save(Model, file=paste(Folder,"mods_for_this_run.RData",sep=""))
    
    # Load parToMinitor from Eric's code
    parToMonitor = which(Model$Estimated)
    parnames = Model$Parameter[parToMonitor] 
    
    ######################
    # Convergence diagnostics
    ######################
    
    # TracePlot
      # maxDims=10; Nkeep=5000; parToMonitor=parToMonitor; parnames=parnames[parToMonitor]; FileName=""; Type="Trace"; Folder=Folder
    ConvergencePlot(McmcArray, maxDims=8, parToMonitor=parToMonitor, parnames=parnames, FileName="", Type="Trace", Folder=Folder, Model=Model)
    # ACF plot
    ConvergencePlot(McmcArray, maxDims=8, Nkeep=500, parToMonitor=parToMonitor, parnames=parnames, FileName="", Type="ACF", Folder=Folder, Model=Model)
    # Variance density plots
    if(length(grep("sigma",parnames)>0)){
      ConvergencePlot(McmcArray, maxDims=8, Nkeep=500, parToMonitor=parToMonitor, parnames=parnames, FileName="", Type="VarianceDensity", Folder=Folder, Model=Model)
    }
    # Correlation density plots
    #if(Model$modelStructure$StrataYear.positiveTows == "correlated" || Model$modelStructure$VesselYear.positiveTows == "correlated") ConvergencePlot(McmcArray, maxDims=8, Nkeep=500, parToMonitor=parToMonitor, parnames=parnames, FileName="", Type="CorrelationDensity", Folder=Folder)
   
    # Convergence statistics -- This doesn't always converge for either Geweke or Gelman-Rubin statistics
    if(McmcDiagnostics==TRUE){
      McmcDiagnosticsPlot(McmcList, parToMonitor, FileName="", Folder=Folder, Geweke=FALSE)
    }
    
    #######################
    # Model fit diagnostics
    #######################
    
    # Visualize the realized offset
    PlotOffset(Data=Data, BugsList=BugsList, FileName="", Folder=Folder)
    # Posterior predictive distribution for positive catches
    PosteriorPredictive(Data=Data, Model=Model, FileName="", Folder=Folder)
    # JAGS indices of abundance
    McmcIndices = ComputeIndices(Data=Data, Model=Model, FileName="", Folder=Folder, Weights=c("StrataAreas","Equal")[1], StrataTable=StrataTable)
    # MLE indices of abundance
    MleIndices = try(ComputeMleIndices(Data=Data, Model=Model, FileName="", Folder=Folder, Weights=c("StrataAreas","Equal")[1], StrataTable=StrataTable, Run=TRUE), silent=TRUE)
    if(inherits(MleIndices, "try-error")==TRUE){
      MleIndices = ComputeMleIndices(Data=Data, Model=Model, FileName="", Folder=Folder, Weights=c("StrataAreas","Equal")[1], StrataTable=StrataTable, Run=FALSE)
    }
    
    # Compare JAGS and MLE
    jpeg(paste(Folder,"/","","Index_Comparison.jpg",sep=""),width=2*3,height=2*3,res=200,units="in")
      par(mfrow=c(2,2), mar=c(2.5,2,2,0), mgp=c(1.25,0.25,0), tck=-0.02)
      matplot(cbind(McmcIndices$Results1[,c('PresMedian','RawPres')], MleIndices$Results1$Pres), col=c("black","red","blue"), lty="solid", type="l", xlab="Strata and/or Year", ylab="Index or component",main="StrataYear Presence", ylim=c(0,max(cbind(McmcIndices$Results1[,c('PresMedian','RawPres')], MleIndices$Results1$Pres),na.rm=TRUE)))
      matplot(cbind(McmcIndices$Results1[,c('PosMedian','RawPos')], MleIndices$Results1$Pos), col=c("black","red","blue"), lty="solid", type="l", xlab="Strata and/or Year", ylab="Index or component",main="StrataYear Positive catch", ylim=c(0,max(cbind(McmcIndices$Results1[,c('PosMedian','RawPos')], MleIndices$Results1$Pos),na.rm=TRUE)))
      matplot(cbind(McmcIndices$Results1[,c('IndexMedian','Raw')], MleIndices$Results1$Index), col=c("black","red","blue"), lty="solid", type="l", xlab="Strata and/or Year", ylab="Index or component",main="StrataYear Index", ylim=c(0,max(cbind(McmcIndices$Results1[,c('IndexMedian','Raw')], MleIndices$Results1$Index),na.rm=TRUE)))
      matplot(cbind(McmcIndices$Results2[,c('IndexMedian','Raw')], MleIndices$Results2$Index), col=c("black","red","blue"), lty="solid", type="l", xlab="Strata and/or Year", ylab="Index or component",main="Year Index", ylim=c(0,max(cbind(McmcIndices$Results2[,c('IndexMedian','Raw')], MleIndices$Results2$Index),na.rm=TRUE)))
    dev.off()
    
    # Save McmcIndices and CV
    Indices[,ModelNumber,1] = McmcIndices$Results2$IndexMean
    Indices[,ModelNumber,2] = McmcIndices$Results2$SdLog  
    for(StratI in 1:nlevels(strata)){
      Which = which(McmcIndices$Results1$Strata==levels(McmcIndices$Results1$Strata)[StratI])
      IndicesByStrata[,StratI,ModelNumber,1] = McmcIndices$Results1$IndexMean[Which]
      IndicesByStrata[,StratI,ModelNumber,2] = McmcIndices$Results1$SdLog[Which]
    }
  }
  
  # Plot Index and CV for all model configurations
  jpeg(paste(SpeciesFolder,"Index_Comparison.jpg",sep=""),width=1*4,height=2*4,res=200,units="in")
    par(mfrow=c(2,1), mgp=c(1.25,0.25,0), mar=c(3,3,1,0), tck=-0.02)
    matplot(Indices[,,1], col="black", lty="solid", type="b", xlab="Year", ylab="Biomass index", ylim=c(0,max(Indices[,,1],na.rm=T)))
    matplot(Indices[,,2], col="black", lty="solid", type="b", xlab="Year", ylab="Index CV", ylim=c(0,max(Indices[,,2],na.rm=T)))
  dev.off()
  
  # Plot Index and CV | Strata for all model configurations
  Log = ""
  jpeg(paste(SpeciesFolder,"Index_Comparison_by_strata.jpg",sep=""),width=2*3,height=nlevels(strata)*3,res=200,units="in")
    par(mfrow=c(nlevels(strata),2), mgp=c(1.25,0.25,0), mar=c(3,3,1,0), tck=-0.02, oma=c(0,2,2,0))
    for(StratI in 1:nlevels(strata)){  
      matplot(IndicesByStrata[,StratI,,1], col="black", log=Log, lty="solid", type="b", xlab="Year", ylab="Biomass index", ylim=list( c(0,max(IndicesByStrata[,,,1])), range(IndicesByStrata[,,,1]) )[[ifelse(Log=="",1,2)]])
      if(StratI==1) mtext(side=3, outer=FALSE, line=1, text="Biomass index", cex=1.5) 
      mtext(side=2, outer=FALSE, line=2, text=levels(strata)[StratI], cex=1.5) 
      matplot(IndicesByStrata[,StratI,,2], col="black", lty="solid", type="b", xlab="Year", ylab="Index CV", ylim=c(0,max(IndicesByStrata[,,,2],na.rm=T)))
      if(StratI==1) mtext(side=3, outer=FALSE, line=1, text="Index CV", cex=1.5)
    }
  dev.off()

}
