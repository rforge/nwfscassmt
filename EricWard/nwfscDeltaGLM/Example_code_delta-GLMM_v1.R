library(stats)
library(runjags)
library(coda)
library(superdiag)
library(R2jags)
library(pscl)
library(statmod)
load.module("glm")
runif(1)

options(stringsAsFactors=TRUE)

my.wd<- SourceFile <- "C:/Users/James.Thorson/Desktop/NWFSC_SVN/EricWard/nwfscDeltaGLM/"
source(paste(my.wd,"bayesGLM v2.15.R",sep=""))
Letters = apply(MARGIN=1,FUN=paste,collapse="",expand.grid(letters,letters))

# read in the master data file
DataFile = SourceFile
setwd(SourceFile)

# Load data and strata
masterDat = read.csv(paste(DataFile,"Example_Species.csv",sep=""))
#LoadFn(paste(DataFile,"RandomSpecies.dmp",sep=""))
strata.limits <- readIn(ncol=5,nlines=6)
  STRATA  NLat SLat MinDepth MaxDepth
  A      49.0 36.0  55        183
  B      49.0 36.0  183       549
  C      49.0 34.5  549       700
  D      36.0 34.5  55        549
  E      34.5 32.0  55        549

# Modify data slightly
species = "Example_Species"
names(masterDat)[9] = species

# Preliminary data processing
processData()

# Define settings
mcmc.control = list(chains=2, thin=2, burnin=2e3, iterToSave=2e3)
#Parallel = TRUE   # If having trouble, try turning off parallel
Parallel = FALSE   # If having trouble, try turning off parallel
#modelStructure1 = list("StrataYear.positiveTows"="fixed", "VesselYear.positiveTows"="random", "StrataYear.zeroTows"="fixed", "VesselYear.zeroTows"="random", "Catchability.positiveTows"="one", "Catchability.zeroTows"="zero", "year.deviations"="fixed", "strata.deviations"="fixed")
modelStructure1 = list("StrataYear.positiveTows"="fixed", "VesselYear.positiveTows"="zero", "Vessel.positiveTows"="zero", "StrataYear.zeroTows"="fixed", "VesselYear.zeroTows"="random", "Vessel.zeroTows"="zero", "Catchability.positiveTows"="one", "Catchability.zeroTows"="zero", "year.deviations"="fixed", "strata.deviations"="fixed")
#modelStructure2 = list("StrataYear.positiveTows"="correlated", "VesselYear.positiveTows"="correlated", "StrataYear.zeroTows"="correlated", "VesselYear.zeroTows"="correlated", "Catchability.positiveTows"="one", "Catchability.zeroTows"="zero", "year.deviations"="fixed", "strata.deviations"="fixed")

# Define models
attach(Data)
mods = list()
mods[[1]] = fitCPUEModel(modelStructure=modelStructure1, mcmc.control=mcmc.control, Parallel=Parallel, Species=species)
#mods[[2]] = fitCPUEModel(modelStructure=modelStructure2, mcmc.control=mcmc.control, Parallel=Parallel, Species=species)

# Process MCMC output
# Make sure that Data is attached prior to running
doMCMCDiags(my.wd,mods)

# Extract DIC:
# SEE NWFSC Assessment Handbook for reasons why J. Thorson and E. Ward think DIC is inappropriate: 
# https://docs.google.com/a/noaa.gov/document/d/1KhQs8Q6q8iPDKdTjEvM7k4t3v0FIWBZmJsJKSMeXSXo/edit#heading=h.4iwd0uxza7v9
mods[[1]]$BUGSoutput$DIC
