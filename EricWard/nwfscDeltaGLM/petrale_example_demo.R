rm(list=ls())                                    
library(stats)
library(runjags)
library(R2jags)
library(coda)
library(superdiag)
library(R2jags)
library(pscl)
load.module("glm") # this is loading a specific library in JAGS / BUGS that implements a conditional sampler
runif(1) # needed for PCs

setwd("/Users/xwarder/Dropbox/delta GLMM project")
my.wd = "bayesGLM_final/" # project wd

# read in the master data file
masterDat = read.csv(paste(my.wd,"Data.csv",sep=""))
masterDat = masterDat[which(masterDat$SURVEY==3),]
strata.limits = read.csv(paste(my.wd,"defaultLimits.csv",sep=""))
  
names(masterDat) # check species available to be run
species = "petrale"
load("bayesGLM_final/petrale.Rdata")

setwd("/Users/xwarder/Dropbox/delta GLMM project")
source(paste(my.wd,"bayesGLM v2.12.r",sep="")) # source the scripts

# call the function to process the data frame:
processData()  
# Set MCMC parameters	 
# Note: the total iterations saved will be chains*iterToSave/thin
mcmc.control = list(chains = 1, thin = 1, burn = 200, iterToSave = 500)
# Set Parallel argument – for PCs only	
Parallel = TRUE
  
# Model 1 contains correlated positive components for strata-year and vessel-year interactions, 
mods = list()
mods[[1]] = fitCPUEModel(modelStructure=list("StrataYear.positiveTows" = "randomExpanded","VesselYear.positiveTows" = "randomExpanded","StrataYear.zeroTows" = "randomExpanded","VesselYear.zeroTows" = "randomExpanded", "Catchability.positiveTows" = "linear", "Catchability.zeroTows" = "linear", "year.deviations" = " uncorrelated ","strata.deviations" = "uncorrelated"),   
mcmc.control=mcmc.control, Parallel=Parallel, Species=species)

mods[[2]] = fitCPUEModel(modelStructure=list("StrataYear.positiveTows" = "random","VesselYear.positiveTows" = "random","StrataYear.zeroTows" = "random","VesselYear.zeroTows" = "random", "Catchability.positiveTows" = "linear", "Catchability.zeroTows" = "linear", "year.deviations" = " uncorrelated ","strata.deviations" = "uncorrelated"),   
mcmc.control=mcmc.control, Parallel=Parallel, Species=species)

# Model 3 is a simple model, with only strata and year effects estimated 
mods[[3]] = fitCPUEModel(modelStructure=list("StrataYear.positiveTows" = "zero","VesselYear.positiveTows" = "zero","StrataYear.zeroTows" = "zero","VesselYear.zeroTows" = "zero", "Catchability.positiveTows" = "one", "Catchability.zeroTows" = "zero", "year.deviations" = " uncorrelated ","strata.deviations" = "uncorrelated"),   
mcmc.control=mcmc.control, Parallel=Parallel, Species=species)

