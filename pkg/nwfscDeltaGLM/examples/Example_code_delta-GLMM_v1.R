

my.wd<- SourceFile <- "C:/Users/James.Thorson/Desktop/"

install.packages("nwfscDeltaGLM", repos="http://R-Forge.R-project.org")
library(nwfscDeltaGLM)

# read in the master data file
DataFile = SourceFile
setwd(SourceFile)

# Load data and strata
data(Example_Species)
masterDat = Example_Species
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
mcmc.control = list(chains=2, thin=1, burnin=1e3, iterToSave=1e3)
Parallel = FALSE   # If having trouble, try turning off parallel
modelStructure1 = list("StrataYear.positiveTows"="fixed", "VesselYear.positiveTows"="random", "Vessel.positiveTows"="zero", "StrataYear.zeroTows"="fixed", "VesselYear.zeroTows"="random", "Vessel.zeroTows"="zero", "Catchability.positiveTows"="one", "Catchability.zeroTows"="zero", "year.deviations"="fixed", "strata.deviations"="fixed")

# Define models
attach(Data)
mods = list()
mods[[1]] = fitDeltaGLM(modelStructure=modelStructure1, mcmc.control=mcmc.control,Parallel=Parallel, Species=species)

# Process MCMC output
# Make sure that Data is attached prior to running
data(SA3)
doMCMCDiags(my.wd,mods)
