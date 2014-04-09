##############################################################################
#
#  Sample program to set up and run PacFIN length and age-at-length comps.
#  Petrale is the sample species
#
#  Andi Stephens, 2012
#
##############################################################################

# Split output to go to console and logfile

sink("Petrale.log.text", split=T)

source("PacFIN_utilities.R")

load("Petrale/PacFIN.PTRL.bds.09.Jun.2011.dmp")

Pdata = PacFIN.PTRL.bds.09.Jun.2011

Catch = read.csv("PacFIN.PTRL.Catch.csv")

#
# If Petrale, add CALCOM data
#
# load or read the data
# then call
#
# Pdata = CombineCALCOM( Pdata, Caldata )
#

# Look at diagnostic summaries and plots

doDiags( Pdata )

# Filter the data and remove bad records

Pdata = cleanPacFIN( Pdata, CLEAN=T)

# Assign state according to method desired.

Pdata = getState( Pdata, source="SOURCE_AGID", CLEAN=T )

# Assign season.  For now, only Petrale seasons defined

Pdata = getSeason( Pdata, seas_type=1 )

# Assign fleet id based on GEAR.  This creates four fleets, the fourth composed
# of two gears.  Any records with other GEAR will be assigned to fleet 0.

Pdata$fleet = Stratify( Pdata$GEAR, splits=list("GFS", c("MDT","TB"), "TR", c("LGL","FTS")) )

# Stratify by DEPTH_AVG range

Pdata$use_depth = Stratify( Pdata$DEPTH_AVG, splits=c(0, 50, 100, 250, 500), range=T)

# Stage 1 expansion

Pdata = getExpansion_1 (Pdata, Sweight="use_wt", Tweight="use_total_wgt", minExp=1, maxExp=0.95)

# Stage 2 expansion

Pdata = getExpansion_2 ( Pdata, Catch, target = c("fleet") )

# Designate the stratification(s)

Pdata$depthrange = Stratify(Pdata$DEPTH_AVG, range=T, splits=c(0,150,350,500))

strat = c("state", "depthrange")

# Explicitly setup expansion variable.
# Choose one of:  exp1, exp2, or exp1*exp2

Pdata$Final_Sample_Sizes = Pdata$Expansion_Factor_1 *  Pdata$Expansion_Factor_2

# Get the length comps before filtering for age.

lenComps = getComps(Pdata, strat=strat, Comps="LEN")

# Simplest sex ratio just assigns them 50/50.

lenComps = doSexRatio(lenComps)
writeComps(lenComps, lbins=xxx, ...)

# Extra filtering for ages

Pdata = cleanAges(Pdata)

# Get the Age Comps according to desired stratification

ageComps = getComps(Pdata, strat=strat, Comps="AGE")
ageComps = doSexRatio(ageComps, ...)
writeComps(ageComps, ageBins=xxx, ...)

# Get the Age-at-Length Comps

aalComps = getComps(Pdata, strat, Comps="AAL")
aalComps = doSexRatio(aalComps, ...)
writeComps(aalComps)


# Close the log file

sink()

#
#  That's All, Folks!
#
##############################################################################
