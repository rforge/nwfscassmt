##############################################################################
#
#  Utilities for working up PacFIN Data.
#
#  Andi Stephens, 2012
#
##############################################################################

##############################################################################
#
# Function paste.col
# 
# Converts a row to a string, "pasting" the columns together.
#
#############################################################################

paste.col <- function(x) {

  # If it's just a vector, return each value as a string.

  if (is.null(dim(x))) {

    return(paste(as.character(x))) 

  } # End if

  # Otherwise, it's a matrix.

  # Get the value of the first column in character form

  out <- paste(as.character(x[, 1])) 

  # Add on each of the rest of the columns
  
  for (i in 2:ncol(x)) { 

   out <- paste(out, as.character(x[, i])) 

  } # End for 

  return(out) 

} # End paste.col
 

#############################################################################
#
# Function find.matching.rows
#
#
#   AUTHOR:  John R. Wallace (John.Wallace@noaa.gov) 
#   REVISED: Andi Stephens, 2010.
# 
#   Takes two tables with a shared primary key, and 
#   returns the rows of the second table for which the 
#   keys match those of the first.
#
#   NOTE:  The way this is used assumes that the second table is a
#          superset of the first (i.e., that each value is matched).
#
#   Changes:  
#
#        Changed name from original "match.f" to "find.matching.rows".
#
#        Removed sub-function 'paste.col' and made it standalone.
#
#        The matching function no longer modifies it's inputs, just
#        returns the values to be 'cbound' in the calling function.
#
#
# Using the primary keys in columns named 'findex' and 'tindex', finds the
# matching values for 'file' in 'table' and returns 'table' column(s) 'tcol'.
#
# Note that no test is made to determine if there are unmatched rows.
#
#############################################################################

find.matching.rows <- function(file, table, findex = 1, tindex = 1, tcol = 2, round. = T) {

  # Coerce a vector argument into a matrix

  if (is.null(dim(file))) {  dim(file) <- c(length(file), 1) }

  # If the primary keys are numeric, round them.

  if (round.) { 

    if (is.numeric(file[, findex])) { file[, findex] <- round(file[, findex]) }

    if (is.numeric(table[, tindex])) { table[, tindex] <- round(table[, tindex]) }

  } # End if round.

  # Convert the indices to character strings for comparison, and get the 
  # positions of the 'file' values in the 'table' values.

  matched.rows = match(paste.col(file[, findex]), paste.col(table[, tindex]))

  # Return the 'tcol' values in the rows of the 'table' that matched.

  return(table[matched.rows, tcol, drop = FALSE]) 

} # End function find.matching.rows


##############################################################################
#
#  printIf outputs comments unless VERBOSE is FALSE. If VERBOSE does not exist,
#  sets it to TRUE in the global environment.
#
##############################################################################

printIf = function( string ) {

  if (! exists("VERBOSE")) assign("VERBOSE", TRUE, pos=.GlobalEnv)
  if (VERBOSE) cat ("\n", string, "\n" )

} # End printIf


##############################################################################
#
#  cleanPacFIN filters out unusable data, and fixes lengths, ages
#  weights, and agemethods.
#
#  Fields to use post-cleaning are:
#
#  length,  age,  agemethod
#
#  The original fields (in the retained data) are left untouched for
#  diagnostics.
#
##############################################################################

cleanPacFIN = function( Pdata, 
                        USINPFC = FALSE, 
                        badRecords=NULL,
                        method=c("B","S",""),
                        sample_type = FALSE,
                        sample_method = FALSE,
                        CLEAN=TRUE) {

  printIf( "Cleaning data" )

  if (!CLEAN) {

    cat("\nGenerating data report only.  No data will be removed.\n")

    Original_data = Pdata
  } 
  
  # Define legal areas

  # This is a legacy comment from Owen's POP code.  Need to investigate.
  # NOTE only a few with "" from 2005 and 2010. NEARLY 100,000 are from CHR 
  # and VCN (through 1978) - a large proportion of early data.

  USinpfc = c("VUS","CL","VN","COL","NC","SC","EU","","CP","EK","MT","PS ")

  # Define fishyr, fleet, fishery and season  -- some assessments manipulate these.

  Pdata$fishyr = Pdata$SAMPLE_YEAR
  Pdata$fleet = 1
  Pdata$fishery = 1
  Pdata$season = 1

  # Fix Lengths

  Pdata$FORK_LENGTH[is.na(Pdata$FORK_LENGTH)] = -1
  Pdata$length = ifelse(Pdata$FISH_LENGTH > -1, Pdata$FISH_LENGTH, Pdata$FORK_LENGTH)

  # Convert mm to cm

  Pdata$lengthcm = floor(Pdata$length / 10)

  # Fix EXP_WT:  Used in expansions

  Pdata$EXP_WT[is.na(Pdata$EXP_WT)] = 0

  # Convert to metric tonnes

  #Pdata$cluster_wgt = Pdata$CLUSTER_WGT / 2205
  #Pdata$total_wgt = Pdata$TOTAL_WGT / 2205
  #Pdata$exp_wt = Pdata$EXP_WT / 2205

  Pdata$cluster_wgt = Pdata$CLUSTER_WGT 
  Pdata$total_wgt = Pdata$TOTAL_WGT 
  Pdata$exp_wt = Pdata$EXP_WT 


  # Fix Ages

  Pdata$age = Pdata$FISH_AGE_YEARS_FINAL
  Pdata$age = ifelse(Pdata$age > 0, Pdata$age, Pdata$age1)
  Pdata$age[is.na(Pdata$age)] = -1
  Pdata$age = ifelse(Pdata$age > 0, Pdata$age, Pdata$age2)
  Pdata$age[is.na(Pdata$age)] = -1
  Pdata$age = ifelse(Pdata$age > 0, Pdata$age, Pdata$age3)
  Pdata$age[is.na(Pdata$age)] = -1

  # Fix Age Methods

  Pdata$agemethod = Pdata$AGE_METHOD
  Pdata$agemethod[Pdata$agemethod == "1"] = "B"
  Pdata$agemethod[Pdata$agemethod == "2"] = "S"
  Pdata$agemethod[is.na(Pdata$agemethod)] = -1

  # We don't want no stinkin' NAs!

  Pdata$SEX[is.na(Pdata$SEX)] = "U"
  Pdata$SEX[Pdata$SEX == 0 ] = "U"

  # Flag records without a SAMPLE_NO

  Pdata$sample = Pdata$SAMPLE_NO
  Pdata$sample[is.na(Pdata$sample)] = "x"

  # Remove records

  Rec_summary = rep(0,6)

  Rec_summary[1] = nrow(Pdata)

  if (USINPFC == TRUE) { Pdata = Pdata[Pdata$INPFC_AREA %in% USinpfc,] }

  Rec_summary[2] = nrow(Pdata)

  Pdata = Pdata[!Pdata$sample %in% badRecords,]

  Rec_summary[3] = nrow(Pdata)

  if (sample_type == TRUE) { Pdata = Pdata[Pdata$SAMPLE_TYPE %in% c(NA,"","M"),] }

  Rec_summary[4] = nrow(Pdata)

  if (sample_method == TRUE) { Pdata = Pdata[Pdata$SAMPLE_METHOD == "R",] } 

  Rec_summary[5] = nrow(Pdata)

  Pdata = Pdata[!is.na(Pdata$length),]

  Rec_summary[6] = nrow(Pdata)

  # Report removals

  cat("\nRemoval Report\n\n")
  cat("Records in input:                 ", Rec_summary[1], "\n")
  cat("Records not in US INPFC_AREA:     ", Rec_summary[1] - Rec_summary[2], "\n")
  cat("Records in badRecords list:       ", Rec_summary[2] - Rec_summary[3], "\n")
  cat("Records with bad SAMPLE_TYPE      ", Rec_summary[3] - Rec_summary[4], "\n")
  cat("Records with bad SAMPLE_METHOD    ", Rec_summary[4] - Rec_summary[5], "\n")
  cat("Records with no usable length     ", Rec_summary[5] - Rec_summary[6], "\n")
  cat("Records remaining:                ", nrow(Pdata), "\n")

  if (CLEAN) {

    return(Pdata)

  } else {

    cat("\n\nReturning original data because CLEAN=FALSE\n\n")

    return(Original_data)

 } # End if-else

} # End cleanPacFIN

#############################################################################
#
# cleanAges removes the samples with bad ages or agemethods.  Depends on the
# data having first been filtered with cleanPacFIN().
#
#############################################################################

cleanAges = function( Pdata, method=c("B","S",""), minAge=-1, maxAge=NULL, CLEAN=TRUE) {

  printIf( "Cleaning data for Ages")

  if (!CLEAN) {

    cat("\nGenerating data report only.  No data will be removed.\n")

    Original_data = Pdata

  }

  # Remove bad records

  Rec_summary = rep(0,3)

  Rec_summary[1] = nrow(Pdata)

  # Records with bad ages

  Pdata = Pdata[Pdata$age > minAge,]
  Rec_summary[2] = nrow(Pdata)

  Pdata = Pdata[Pdata$agemethod %in% method,]
  Rec_summary[3] = nrow(Pdata)

  if ( ! is.null(maxAge) ) {

    cat("\nSetting maximum age to", maxAge, "\n")

    Pdata$age[Pdata$age > maxAge] = maxAge

  } # End if

  # Report removals

  cat("\nRemoval report\n\n")
  cat("Records in input:                  ", Rec_summary[1], "\n")
  cat("Records with age less than min:    ", Rec_summary[1] - Rec_summary[2], "\n")
  cat("Records with bad agemethods:       ", Rec_summary[2] - Rec_summary[3], "\n")
  cat("Records remaining:                 ", nrow(Pdata), "\n")

  if (CLEAN) {

    return(Pdata)

  } else {

    cat("\nReturning original data because CLEAN=FALSE\n\n")
    return(Original_data)

  } # End if-else

} # End cleanAges

##############################################################################
#
#  getState creates a state field that is developed from SOURCE_AGID, from
#  PSMFC_ARID, or from SAMPLE_AGENCY depending on the option selected.
#
#  If CLEAN is TRUE, records for which a state could not be assigned are
#  removed.
#
##############################################################################

getState = function ( Pdata, source="SOURCE_AGID", CLEAN=TRUE) {

  printIf( paste("Getting state information from", source ))

  sources = c("SOURCE_AGID", "PSMFC_ARID", "SAMPLE_AGENCY" )

  if ( ! source %in% sources) {

    stop( "Legitimate sources for getState are: ", sources, "\n" )

  } # End if

  if (!CLEAN) {

    cat("\nGenerating data report only.  No data will be removed.\n\n")

    Original_data = Pdata

  }

  Pdata$state = Pdata[, source]

  if ( source == "SOURCE_AGID" | source == "SAMPLE_AGENCY" ) {

    Pdata$state[Pdata$state == "C"] = "CA"
    Pdata$state[Pdata$state == "O"] = "OR"
    Pdata$state[Pdata$state == "W"] = "WA"

  } # End if

  if ( source == "PSMFC_ARID" ) {

    Pdata$state[Pdata$state %in% c("3A","3B","3S")] = "WA"
    Pdata$state[Pdata$state %in% c("2A","2B","2C")] = "OR"
    Pdata$state[Pdata$state %in% c("1A","1B","1C","CAL")] = "CA"

  } # End if PSMFC_ARID

  states = c("OR","CA","WA")

  nostate = sum(! Pdata$state %in% states )

  if (CLEAN) {

    Pdata = Pdata[Pdata$state %in% states,]

    cat("\n\n", nostate, " records were removed because no state id could be assigned\n\n")

  } else {

    cat("There are ", nostate, " records for which no state id could be assigned\n")
    cat("\nReturning original data because CLEAN=FALSE\n\n")

  } # End if CLEAN

  return( Pdata )

} # end getState

##############################################################################
#
#  getSeason adds a field for season to the data.  Several seasonal
#  schemes may be provided, including the Petrale seasons ( 1 = winter months,
#  2 else ).
#
#  Others as needed.
#
##############################################################################

getSeason = function ( Pdata, seas_type=-1) {

  Pdata$season = 1

  if ( seas_type == 0 ) {

    printIf( "Assigning season from SAMPLE_MONTH")

    Pdata$season = as.numeric(Pdata$SAMPLE_MONTH)

  } # End if

  # Petrale seasons

  if ( seas_type == 1 ) {

    printIf( "Assigning season ala Petrale")

    Pdata$season = 2
    Pdata$season[Pdata$SAMPLE_MONTH %in%  c(11,12,1,2)] = 1

  } # End if Petrale

  return( Pdata )

} # End getSeason 

##############################################################################
#
# Stratify takes an input vector and list of values used to designate
# the values returned in the strats vector.
#
# Any values not given in splits will be assigned to stratum 0.
#
# If ranges=T, this defaults to a call to findInterval, and stratum 0
# will be assigned to values smaller than the first value in splits.
#
# Note that this can be used to designate fleet or for stratification
# based on depth or INPFC area.
#
##############################################################################

Stratify = function ( inVector=NULL, splits=NULL, range=F ) {

  strats = rep(0, length(inVector))

  if (!range) {

    inVector = as.character(inVector)

    for ( i in 1:length(splits) ) {

      splitby = as.character(splits[[i]])
      strats[inVector %in% splitby] = i

    } # End for

  } else {

    strats = findInterval(inVector, splits)

  } # End if-else

  return(strats)

} # End function Stratify

##############################################################################
#
#  CombineCALCOM combines CALCOM and PacFIN data (only used for Petrale).
#
##############################################################################

CombineCALCOM = function ( Pdata, CALCOM ) {

  printIf( "Combining CALCOM and PacFIN data" )

  # Fix dates

  CALCOM$SAMPLE_DATE = as.character(CALCOM$SAMPLE_DATE)

  # Break out Year, Month, Day from vector formatted "2/23/2012"

  trueDate = as.Date(CALCOM$SAMPLE_DATE, format="%m/%d/%Y")

  CALCOM$SAMPLE_YEAR = as.numeric(format(trueDate, format="%Y"))
  CALCOM$SAMPLE_MONTH = as.numeric(format(trueDate, format="%m"))
  CALCOM$SAMPLE_DAY = as.numeric(format(trueDate, format="%d"))

  # Fix Areas

  CALCOM$PSMFC_AREA = NA
  CALCOM$PSMFC_AREA[CALCOM$PORT %in% c("ERK","CRS")] = "1C"
  CALCOM$PSMFC_AREA[CALCOM$PORT %in% c("BRG","OSF","MNT","1")] = "1B"
  CALCOM$PSMFC_AREA[CALCOM$PORT %in% c("OSB","MRO","2")] = "1A"

  # Create PacFIN format matrix and fill

  CAL.dat = as.data.frame(matrix(data=NA, nrow = nrow(CALCOM) , ncol = ncol(Pdata)))

  names(CAL.dat) = names(Pdata)
  
  CAL.dat$SPID         = CALCOM$SPECIES
  CAL.dat$SAMPLE_NO    = CALCOM$SAMPLE_NO
  CAL.dat$FISH_NO      = as.numeric(CALCOM$FISH_NO)
  CAL.dat$FISH_LENGTH  = as.numeric(CALCOM$TLENGTH)
  CAL.dat$SEX          = CALCOM$SEX
  CAL.dat$DEPTH_AVG    = as.numeric(CALCOM$DEPTH)
  CAL.dat$TOTAL_WGT    = as.numeric(CALCOM$TOTAL_WGT)
  CAL.dat$PORT         = CALCOM$PORT_COMPLEX
  CAL.dat$SAMPLE_YEAR  = as.numeric(CALCOM$SAMPLE_YEAR)
  CAL.dat$SAMPLE_MONTH = as.numeric(CALCOM$SAMPLE_MONTH)
  CAL.dat$SAMPLE_DAY   = as.numeric(CALCOM$SAMPLE_DAY)
  CAL.dat$SOURCE_AGID  = "C"
  CAL.dat$age1         = as.numeric(CALCOM$AGE)
  CAL.dat$FISH_AGE_YEARS_FINAL = as.numeric(CALCOM$AGE)

  # Fix SEX.  2=females, 1=males based on length distributions

  CAL.dat$SEX[CAL.dat$SEX=="1"] = "M"
  CAL.dat$SEX[CAL.dat$SEX=="2"] = "F"
  
  # Done.
  
  return(rbind(Pdata,CAL.dat))
  
} # End CombineCALCOM
  
##############################################################################
#
# getExpansion_1 prepares the values that are used to calculate the level-1
# expansion factor.  Indiv_Wgts controls whether or not individual weight
# calculations will be used for the expansion, which is intended to account
# for the unsampled fish in the tow. maxExp is either a number or a quantile
# determining the maximum value the expansion factor can take.
#
# The weight of the fish is calculated three ways.
#
# First, FEMALES_WGT and MALES_WGT are summed per SAMPLE_NO.
#
# The SPECIES_WGT is summed across all clusters in a sample to provide
# the second per SAMPLE_NO weight.
#
# Finally, weights of the male and female fish are calculated from their 
# lengths.  These are summed per SAMPLE_NO to provide a per-sample weight.
# Zero weights might occur for some fish. These are filled in with the 
# median weight of fish in the sample.
#
# Use_acs is all_cluster_sums with missing values filled in with state-
# specific medians.
#
# Species_Percent_Sampled is the percentage of this species in the samples,
# calculated as Trip_Sampled_Lbs divided by Use_acs.
#
# Trip_Sampled_Lbs is calculated differently for each state:
#
#     For California, Trip_Sampled_Lbs = Species_Percent_Sampled * TOTAL_WGT.
#     For Oregon, Trip_Sampled_Lbs = EXP_WT.  Where missing, use Species_Percent
#     _Sampled, as above.  For Washington, use RWT_LBS, TOTAL_WGT, median(RWT_LBS)
#     or median(TOTAL_WGT).
#
#     If all else fails, use state-specific medians.
#
# No original data columns are altered in this process.  The original
# data are supplemented with new columns containing the values to be used.
#
# New columns:
#
#    Wt_Sampled_1	per-SAMPLE_NO summed MALES_WGT + FEMALES_WGT
#    Wt_Sampled_2	per-SAMPLE_NO summed cluster_wt
#    LW_Calc_Wt		individual weights predicted from L-W
#    Wt_Sampled_3	per-SAMPLE_NO summed LW_Calc_Wts
#
#    Wt_Sampled		per-SAMPLE_NO combined weights.
#                       This is preferentially Wt_Sampled_1, with NAs
#                       replaced with values from Wt_Sampled_2, and NAs
#			remaining replaced with values from Wt_Sampled_3.
#
#    Wt_Method		Denotes by which method the Wt_Sampled was filled.
#    Numlens		per-SAMPLE_NO number of lengths
#
#    Use_acs  		all_cluster_sum with NAs replaced by state-specific medians.
#    Species_Percent_Sampled  Percentage of this species in samples.
#    Trip_Sampled_Lbs	State-specific sampled weights.
#
#
##############################################################################
  
getExpansion_1 = function( Pdata, maxExp = 0.95, Indiv_Wgts = T,
                           fa=2e-06, fb=3.5, ma=2e-06, mb=3.5, ua=2e-06, ub=3.5) {


  # Assigned 'state' if it's not already there.

  if (length(which(names(Pdata) == "state")) == 0) { Pdata = getState(Pdata, CLEAN=T) }

  # Everything is calculated in terms of unique samples.

  tows = Pdata[!duplicated(Pdata$SAMPLE_NO),]
 
  printIf("Getting the summed weights of the sexed fish for each SAMPLE_NO")
  
  tows$Wt_Sampled_1 = tows$MALES_WGT + tows$FEMALES_WGT

  Pdata$Wt_Sampled_1 = tows$Wt_Sampled_1[match(Pdata$SAMPLE_NO, tows$SAMPLE_NO)]

  printIf("Summing the SPECIES_WEIGHTS for all clusters in a SAMPLE_NO")

  Pdata$SAMP_CLUST = paste(Pdata$SAMPLE_NO, Pdata$CLUSTER_NO, sep="_")
  uniqueClusters = Pdata[!duplicated(Pdata$SAMP_CLUST),]

  tmp = aggregate(uniqueClusters$SPECIES_WGT, list(uniqueClusters$SAMPLE_NO), sum, na.rm=T)
  names(tmp) = c("SAMPLE_NO", "wgt")

  # Might have 0 values that should actually be NA.  Aggregate doesn't generate NAs.

  tmp$wgt[tmp$wgt == 0] = NA

  tows$Wt_Sampled_2 = tmp$wgt[match(tows$SAMPLE_NO,tmp$SAMPLE_NO)]

  Pdata$Wt_Sampled_2 = tmp$wgt[match(Pdata$SAMPLE_NO,tmp$SAMPLE_NO)]

  # Only if there are individual weight factor and coefficients available

  if (Indiv_Wgts) {

    printIf( paste( "Creating expansion factor based on M/F weights"))

    ############################################################################
    #
    # Create a predicted fish weight based on sex and length
    # (use mm for fitted regression coefficients!)
    # these will be summed to give the sample weight
    #
    ############################################################################
  
    Pdata$LW_Calc_Wt = NA
    Pdata$LW_Calc_Wt[Pdata$SEX=="F"] = fa*Pdata$length[Pdata$SEX=="F"]^fb
    Pdata$LW_Calc_Wt[Pdata$SEX=="M"] = ma*Pdata$length[Pdata$SEX=="M"]^mb
    Pdata$LW_Calc_Wt[Pdata$SEX=="U"] = ua*Pdata$length[Pdata$SEX=="U"]^ub
  
    # Convert to pounds
  
    Pdata$LW_Calc_Wt = Pdata$LW_Calc_Wt * 0.00220462
  
    # Get the number of observed lengths and weights to use for each sample
  
    tmp = as.data.frame(table(Pdata$SAMPLE_NO))
    names(tmp) = c("SAMPLE_NO","numlens")
  
    tows$numlens = tmp$numlens[match(tows$SAMPLE_NO, tmp$SAMPLE_NO)]
  
    # Note:  if all else fails, fill 0 weights with the median in that sample.

    tmp_wt = aggregate(Pdata$LW_Calc_Wt, list(Pdata$SAMPLE_NO), sum, na.rm=T)
    med_wt = aggregate(Pdata$LW_Calc_Wt, list(Pdata$SAMPLE_NO), median, na.rm=T)
    tmp_wt = cbind(tmp_wt, med_wt[,2])
    names(tmp_wt) = c("SAMPLE_NO","Wt_Sampled_3","Median")

    tmp_wt$Wt_Sampled_3[tmp_wt$Wt_Sampled_3 == 0] = tmp_wt$Median[tmp_wt$Wt_Sampled_3 == 0]
  
    tows$Wt_Sampled_3 = tmp_wt$Wt_Sampled_3[match(tows$SAMPLE_NO, tmp_wt$SAMPLE_NO)]

    Pdata$Wt_Sampled_3 = tows$Wt_Sampled_3[match(Pdata$SAMPLE_NO, tows$SAMPLE_NO)] 
    Pdata$numlens = tows$numlens[match(Pdata$SAMPLE_NO, tows$SAMPLE_NO)] 
  
  } # End if Indiv_Wgts


  ############################################################################
  #
  # Use calculated weights for Wt_Sampled.
  #
  ############################################################################

  tows$Wt_Sampled = tows$Wt_Sampled_1
  tows$Wt_Method = 1

  tows$Wt_Method[is.na(tows$Wt_Sampled)] = 2
  tows$Wt_Sampled[is.na(tows$Wt_Sampled)] = tows$Wt_Sampled_2[is.na(tows$Wt_Sampled)]

  if (Indiv_Wgts) {

    tows$Wt_Method[is.na(tows$Wt_Sampled)] = 3
    tows$Wt_Sampled[is.na(tows$Wt_Sampled)] = tows$Wt_Sampled_3[is.na(tows$Wt_Sampled)]

  } # End if

  tows$Wt_Method[is.na(tows$Wt_Sampled)] = NA

  Pdata$Wt_Sampled = tows$Wt_Sampled[match(Pdata$SAMPLE_NO,tows$SAMPLE_NO)]

  printIf("Done calculating sample weights")

  print(summary(Pdata$Wt_Sampled))

  printIf("Removing records with NA sample weights")

  Pdata = Pdata[!is.na(Pdata$Wt_Sampled),]

  ############################################################################ 
  #
  # Now do the all_cluster_sums, replacing NA with medians.
  #
  ############################################################################# 

  printIf("Getting cluster sums")

  tows$Use_acs = tows$all_cluster_sum

  tows$Use_acs[tows$Use_acs == 0] = NA

  # Use median from the right agency where total landed weight missing

  medallcls.or = median(tows$Use_acs[tows$state == "OR"], na.rm=T)
  medallcls.ca = median(tows$Use_acs[tows$state == "CA"], na.rm=T)
  medallcls.wa = median(tows$Use_acs[tows$state == "WA"], na.rm=T)

  tows$Use_acs[is.na(tows$Use_acs) & tows$state=="OR"] = medallcls.or
  tows$Use_acs[is.na(tows$Use_acs) & tows$state=="CA"] = medallcls.ca
  tows$Use_acs[is.na(tows$Use_acs) & tows$state=="WA"] = medallcls.wa

  # Might be all NA for one state, get the overall median

  medallcls = median(c(medallcls.or, medallcls.ca, medallcls.wa), na.rm=T)

  tows$Use_acs[is.na(tows$Use_acs)] = medallcls

  Pdata$Use_acs = tows$Use_acs[match(Pdata$SAMPLE_NO, tows$SAMPLE_NO)]
 
  ############################################################################ 
  #
  # Calculate Species_Percent_Sampled.
  #
  ############################################################################ 

  tows$Species_Percent_Sampled = tows$Wt_Sampled/tows$Use_acs

  tows$Use_Percent = tows$Species_Percent_Sampled * tows$TOTAL_WGT

  ############################################################################ 
  #
  # Get total weight per SAMPLE_NO, replacing NAs with medians.  For Oregon,
  # this is exp_wt.  For CA and OR, fill NAs with Species_Percent_Sampled *
  # TOTAL_WGT.
  #
  ############################################################################ 

  printIf("Getting total weights per sample")

  tows$Trip_Sampled_Lbs = tows$TOTAL_WGT

  # California

  tows$Trip_Sampled_Lbs[tows$state=="CA"] = 
                        tows$Use_Percent[tows$state=="CA"]

  medtotwgt.ca = median(tows$Trip_Sampled_Lbs[tows$state == "CA"], na.rm=T)
  tows$Trip_Sampled_Lbs[is.na(tows$Trip_Sampled_Lbs) & tows$state=="CA"] = medtotwgt.ca

  # Washinton

  tows$Trip_Sampled_Lbs[tows$state=="WA"] = tows$RWT_LBS[tows$state=="WA"]
  tows$Trip_Sampled_Lbs[tows$state=="WA" & is.na(tows$RWT_LBS)] = 
                       tows$TOTAL_WGT[tows$state=="WA" & is.na(tows$RWT_LBS)]

  medtotwgt1.wa = median(tows$RWT_LBS[tows$state == "WA"], na.rm=T)
  medtotwgt2.wa = median(tows$TOTAL_WGT[tows$state == "WA"], na.rm=T)
  
  tows$Trip_Sampled_Lbs[is.na(tows$Trip_Sampled_Lbs) & tows$state=="WA"] = medtotwgt1.wa
  tows$Trip_Sampled_Lbs[is.na(tows$Trip_Sampled_Lbs) & tows$state=="WA"] = medtotwgt2.wa

  # Oregon 
  # OR data is always missing total weights between 1965-1970
  # Use Species_Percent_Sampled, as for CA.

  tows$Trip_Sampled_Lbs[tows$state == "OR" & tows$exp_wt > 0] = 
                 tows$exp_wt[tows$state == "OR" & tows$exp_wt > 0]

  tows$Trip_Sampled_Lbs[tows$Trip_Sampled_Lbs == 0] = NA

  tows$Trip_Sampled_Lbs[tows$state=="OR" & is.na(tows$Trip_Sampled_Lbs)] = 
			tows$Use_Percent[tows$state=="OR" & is.na(tows$Trip_Sampled_Lbs)]

  medtotwgt.or = median(tows$Trip_Sampled_Lbs[tows$state == "OR"], na.rm=T)
  tows$Trip_Sampled_Lbs[is.na(tows$Trip_Sampled_Lbs) & tows$state=="OR"] = medtotwgt.or

  # Match Trip_Sampled_Lbs to the larger dataset.

  Pdata$Trip_Sampled_Lbs = tows$Trip_Sampled_Lbs[match(Pdata$SAMPLE_NO, tows$SAMPLE_NO)]

  print(summary(Pdata$Trip_Sampled_Lbs)); flush.console()

  Pdata$Expansion_Factor_1 = Pdata$Trip_Sampled_Lbs / Pdata$Wt_Sampled

  Pdata$Expansion_Factor_1[Pdata$Expansion_Factor_1 < 1] = 1

  if ( maxExp > 1 ) {

    max.val = maxExp

  } else {

    max.val = quantile(Pdata$Expansion_Factor_1, maxExp, na.rm=T)

  } # End if-else


  Pdata$Expansion_Factor_1[Pdata$Expansion_Factor_1 > max.val] = max.val

  printIf("Expansion Factor 1:")

  print(summary(Pdata$Expansion_Factor_1))

  # Diagnostic plots, summaries Needs to be rewritten.

  #Exp_1_Diags(tows) 

  return(Pdata)

} # End function getExpansion_1

##############################################################################
#
#  A few boxplots for expansion factors
#  Not rewritten since the expansion was rewritten -- needs work.
#
##############################################################################

Exp_1_Diags = function(tows) {

  # Diagnostics for PrepFinal_Sample_Sizes

  boxplot(tows[,c("Wt_Sampled","Wt_Sampled_3","Wt_Sampled_1","SPECIES_WGT","LW_Calc_Wt")])
  print(summary(tows[,c("Wt_Sampled","Wt_Sampled_3","Wt_Sampled_1","SPECIES_WGT","LW_Calc_Wt")]))
  boxplot(tows[,c("total_wgt","Trip_Sampled_Lbs","all_cluster_sum","Use_acs")])
  print(summary(tows[,c("total_wgt","Trip_Sampled_Lbs","all_cluster_sum","Use_acs")]))

} # End function Exp_1_Diags

#############################################################################
#
# getExpansion_2 is the second-stage expansion
#
# Calculate the stratified sampled biomass, All_Trips_Sampled_Lbs by summing
# Trip_Sampled_Lbs. Calculate the stratified catch by summing MT * 2204.
#
# The per-trip, per-stratum Expansion_Factor_2 is the catch / sampled catch.
#
#############################################################################

getExpansion_2 = function( Pdata, Catch, target, maxExp=0.95 ) {

  Pdata$EF2 = 1

  # Per individual sample

  tows = Pdata[!duplicated(Pdata$SAMPLE_NO),]

  # Get the stratified totals for each level of stratification

  for ( i in sort(unique(tows[,target])) ) {

    for ( yr in sort(unique(tows$fishyr)) ) {

      #printIf(paste("Working on Year, level", yr, i, sep="  "))

      thiscatch = Catch[Catch[,target] == i & Catch$Year == yr,]

      if ( nrow(thiscatch) == 0 ) {

          printIf(paste("No records for Year, level", yr, i, sep="  "))
          next()

      } # End if

      for ( s in sort(unique(tows$state)) ) {

        sumSampled = sum(tows$Trip_Sampled_Lbs[tows[,target] == i &
                                               tows$fishyr == yr &
                                               tows$state == s], na.rm=T)

        printIf(paste("Thiscatch, sumSampled, state", thiscatch[,s], sumSampled, s, sep="  "))

        # Convert MT to Lbs by multiplying by 2204

        tows$EF2[tows[,target] == i & 
                 tows$fishyr == yr & tows$state == s] = 2204 * thiscatch[,s] / sumSampled 

      } # End for s

    } # End for yr

  } # End for i

  if ( maxExp > 1 ) {

    max.val = maxExp

  } else {

    max.val = quantile(tows$EF2, maxExp)

  } # End if-else

  tows$EF2[tows$EF2 > max.val] = max.val
  tows$EF2[tows$EF2 < 1] = 1

  print(summary(tows$EF2))

  # Match EF2 to the larger dataset.

  Pdata$Expansion_Factor_2 = tows$EF2[match(Pdata$SAMPLE_NO, tows$SAMPLE_NO)]

  return(Pdata)
  
} # End function getExpansion_2


#############################################################################
#
# doSexRatio determines gender for unsexed fish depending on whether the 
# Rvector is a single ratio, or the vector of ratios per length.
#
# The input data expected are data that have already been aggregated.
#
# Lbins are the bins corresponding to the vector of ratios (unused for
# a single ratio).
#
# If NOT using a vector of lengths and corresponding ratios, then ratioU
# is the ratio applied to unsexed fish less than or equal to maxsizeU.
# This is usually 0.5, since these are fish so small that sexing them is
# probably not done correctly.
#
# GTsizeU is the size above which the ratio is assumed to be 1.0 (big mamas).
#
#############################################################################

doSexRatio = function( CompData, Rvector=0.5, Lbins=NULL, ratioU=.5,
                       maxsizeU=0, GTsizeU=Inf ) {

  # Fix arithmetic

  CompList = c("male","msamps","mtows","female","fsamps","ftows","unsexed","usamps","utows")
  tmp = CompData[,CompList]
  tmp[is.na(tmp)] = 0
  CompData[,CompList] = tmp

  if ( maxsizeU != 0 ) {

    # Create vectors from the appropriate values.

    Rvector = c(ratioU, Rvector, 1)
    Lbins = c(0, maxsizeU, GTsizeU)

  } # End if

  if ( length(Rvector) > 1 ) {

    tmp = paste(Rvector, collapse = ",")
    printIf(paste("Applying sex ratio:", tmp, "to numbers, samples and tows"))

    # Recode lengths to correspond to bins given

    lens = findInterval(CompData$length, Lbins, rightmost.closed=T)

    # Should be unnecessary

    if ( min(lens) == 0 ) {

      lens[lens == 0] = 1
      cat("Adjusting lengths below lbins up into first bin")

    } # End if

    for ( i in sort(unique(lens)) ) {

      CompData$female[lens == i] = CompData$female[lens == i] +
                                   CompData$unsexed[lens == i] *
                                   Rvector[i] 

      CompData$fsamps[lens == i] = CompData$fsamps[lens == i] +
                                   CompData$usamps[lens == i] *
                                   Rvector[i] 

      CompData$ftows[lens == i] = CompData$ftows[lens == i] +
                                  CompData$utows[lens == i] *
                                  Rvector[i] 

      CompData$male[lens == i] = CompData$male[lens == i] +
                                 CompData$unsexed[lens == i] *
                                 (1 - Rvector[i])

      CompData$msamps[lens == i] = CompData$msamps[lens == i] +
                                   CompData$usamps[lens == i] *
                                   (1 - Rvector[i])

      CompData$mtows[lens == i] = CompData$mtows[lens == i] +
                                  CompData$utows[lens == i] * 
                                  (1 - Rvector[i])

    } # End for

  } else { 

    tmp = paste(Rvector, collapse = ",")
    printIf(paste("Applying sex ratio:", tmp, "to numbers, samples and tows"))

    # apply a single ratio over all lengths

    CompData$female = CompData$female + Rvector * CompData$unsexed
    CompData$fsamps = CompData$fsamps + Rvector * CompData$usamps
    CompData$ftows = CompData$ftows + Rvector * CompData$utows

    CompData$male = CompData$male + (1 - Rvector) * CompData$unsexed
    CompData$msamps = CompData$msamps + (1 - Rvector) * CompData$usamps
    CompData$mtows = CompData$mtows + (1 - Rvector) * CompData$utows

  } # End if-else

  return(CompData)

} # End doSexRatio

##############################################################################
#
# Write out comps
#
# Expecting composition data with columns for male and female numbers.
#
# May be ID'd as to state, depth, etc. factors.
#
##############################################################################

writeComps = function(inComps, fname="out.csv", abins=NULL, lbins=NULL, 
                      maxAge=Inf, partition=0, ageErr=0) {

  # Unsexed fish should represent males + females

  cat(paste("Writing comps to file", fname, "\n", sep=" "))
  flush.console()

  inComps$unsexed = inComps$male + inComps$female
  inComps$usamps = inComps$msamps + inComps$fsamps
  inComps$utows = inComps$mtows + inComps$ftows

  # Which comps are we doing?

  Names = names(inComps)
  AGE = which(Names == "age")
  LEN = which(Names == "lengthcm")

  # Fix length bins

  if ( !is.null(inComps$lengthcm) ) {

    if ( is.null(lbins) ) {

      printIf("No length bins provided, using data as-is")

      lbins = sort(unique(inComps$lengthcm))

    } # End if

    # Re-code actual lengths to be lbins

    if ( min(lbins) > 0 ) { lbins = c(0,lbins) }

    lbins = c(lbins, (max(lbins) + 1))
    
    inComps$lbin = findInterval(inComps$lengthcm, lbins, all.inside=T)

    LbinLo = lbins
    LbinHi = lbins[-1] - 1

    # Last length bin is supposed to include the largest bin, right?

    LbinHi[length(LbinHi)] = max(lbins)

    # Note that the last bin is a dummy bin

    LbinHi = c(LbinHi, (max(LbinLo) + 1))

    cat("Lbins:\n\n")
    cat(LbinLo, "\n")
    cat(LbinHi, "\n\n")
    cat("Note that last bin is a dummy bin\n\n")

  } # End if

  # Fix age bins

  if ( !is.null(inComps$age) ) {

    if ( is.null(abins) ) {

      printIf("No age bins provided, using data as-is")

      abins = sort(unique(inComps$age))
    
      abins = abins[abins < maxAge]

    } # End if

    # Re-code actual ages to be abins

    if ( min(abins) > 0 ) { abins = c(0,abins) }

    # add extra, dummy bin because all.inside=T

    abins = c(abins, (max(abins) + 1))
   
    inComps$abin = findInterval(inComps$age, abins, all.inside=T)

    cat("Abins:\n\n")
    cat(abins, "\n\n")
    cat("Note that last bin is a dummy bin\n\n")

  } # End if

  AAL = FALSE

  if ( length(AGE) > 0 ) {

    target = "abin"

    STRAT = AGE-1

    KeyNames = c(Names[1:STRAT])
    inComps$key = paste.col(inComps[,KeyNames])

    # matrix will be Ages, Ntows, Nsamps.
    # it gets re-ordered later.

    NCOLS = 2 + length(abins)
    OutNames = c(paste("A",abins, sep=""), "Ntows","Nsamps")

    if ( length(LEN) > 0 ) {

      AAL = TRUE

      STRAT = AGE-2

      KeyNames = c(Names[1:STRAT], "lbin")
      inComps$key = paste.col(inComps[,KeyNames])

      cat("\n\nAge-at-length takes awhile to assemble.  Be patient!\n")
      flush.console()

      # matrix will be Ages, LbinLo, LbinHi, Ntows, Nsamps.
      # it gets re-ordered later.

      NCOLS =  4 + length(abins)
      OutNames = c(paste("A",abins,sep=""), "lbin","Ntows","Nsamps")

    } # End if

    # Note that ages run from 0, but output columns numbers start at 1.
    # Adjust the "abins" to match the columns.  The column names are
    # already correct.

    inComps$abin = inComps$abin + 1 

  } else {

    target = "lbin"

    STRAT = LEN-1

    KeyNames = c(Names[1:STRAT])
    inComps$key = paste.col(inComps[,KeyNames])

    # matrix will have Lbins, Ntows, Nsamps
    # it gets re-ordered later.

    NCOLS = 2 + length(lbins)
    OutNames = c(paste("L",lbins, sep=""), "Ntows","Nsamps")

  } # End if-else

  # Rename columns to be used below

  names(inComps)[which(names(inComps) == "female")] = "f"
  names(inComps)[which(names(inComps) == "male")] = "m"
  names(inComps)[which(names(inComps) == "unsexed")] = "u"

  # We'll work key by key

  uKeys = inComps$key[!duplicated(inComps$key)]

  # Save the matching stratification columns

  uStrat = inComps[!duplicated(inComps$key), 1:STRAT]

  cat(length(uKeys), "unique keys for", nrow(inComps), "records\n\n")
  flush.console()

  head(inComps)

  cat("\n\n")
  flush.console()

  # For each gender in turn

  for ( g in c("m","f","u")) {

    cat(paste("Assembling, sex is:", g, "\n", sep=" "))
    flush.console()

    tows = which(names(inComps) == paste(g,"tows",sep=""))
    samps = which(names(inComps) == paste(g,"samps",sep=""))

    # Create output matrix

    output = data.frame(matrix(nrow=length(uKeys), ncol=NCOLS, 0))

    names(output) = OutNames

    for ( k in 1:length(uKeys) ) {

      # Get the matching records

      slice = inComps[inComps$key == uKeys[k],]

      if ( AAL ) {

        output$lbin[k] = slice$lbin[1]

      } # End if

      output$Nsamps[k] = sum(slice[,samps], na.rm=T)
      output$Ntows[k] = sum(slice[,tows], na.rm=T)

      for ( s in 1:length(slice[,target]) ) {

          index = slice[s,target]
          output[k,index] = slice[s,g]

      } # End for s 

    } # End for k

    # Save and identify

    output[is.na(output)] = 0

    if ( AAL ) {
    
      output$LbinLo = LbinLo[output$lbin]
      output$LbinHi = LbinHi[output$lbin]

    } # End if

    assign(paste(g, "Comps", sep=""), output)

  } # End for g

  # Now assemble everything and write to a file
  # Note that we are stripping the last, dummy bin.
 
  NCOLS = ifelse(AAL, NCOLS-5, NCOLS-3)

  blanks = mComps[1:NCOLS]
  blanks[,] = 0

  # Stratification puts fleet before fishyr, which is the correct order
  # for the strat, but wrong for the SS input file.  Move it.
  #
  # correctOrder = c("fishyr", "season", "fleet", "gender", "partition",
  #                  "ageErr", "LbinLo", "LbinHi", "Nsamp")

  fleetWas = which(names(uStrat) == "fleet")

  tmp = uStrat$fleet
  uStrat$fleet = uStrat$fishyr
  uStrat$fishyr = tmp

  names(uStrat)[fleetWas] = "fishyr"
  names(uStrat)[(fleetWas + 1)] = "fleet"

  # Fill the rest of the values

  uStrat$gender = 0
  uStrat$partition = partition
  uStrat$ageErr = ageErr
  
  if ( AAL ) {

    # Note that until empty rows are removed, the LbinLo and LbinHi columns
    # are the same in each dataset

    uStrat$LbinLo = fComps$LbinLo
    uStrat$LbinHi = fComps$LbinHi

  }

  Alltows = rowSums(cbind(fComps$Ntows, mComps$Ntows), na.rm=T)
  Allsamps = rowSums(cbind(fComps$Nsamps, mComps$Nsamps), na.rm=T)

  FthenM = cbind(uStrat, Alltows, Allsamps, fComps[,1:NCOLS], mComps[,1:NCOLS])

  Mout = cbind(uStrat, mComps$Ntows, mComps$Nsamps, blanks, mComps[1:NCOLS])
  Fout = cbind(uStrat, fComps$Ntows, fComps$Nsamps, fComps[1:NCOLS], blanks)
  Uout = cbind(uStrat, uComps$Ntows, uComps$Nsamps, uComps[1:NCOLS], blanks)


  # Make it pretty

  index = which(names(Fout) == "fComps$Ntows")

  names(Mout)[index] = "Ntows"
  names(Fout)[index] = "Ntows"
  names(Uout)[index] = "Ntows"

  names(Mout)[index + 1] = "Nsamps"
  names(Fout)[index + 1] = "Nsamps"
  names(Uout)[index + 1] = "Nsamps"

  # Remove empty rows

  Fout = Fout[Fout$Nsamps > 0,]
  Mout = Mout[Mout$Nsamps > 0,]

  Fout$gender=1
  Mout$gender=2
  Uout$gender=3

  # Print the whole shebang out to a file.

  # Turn off meaningless warnings.
  oldwarn = options("warn")
  options("warn" = -1)

  cat("Writing FthenM, dimensions:", dim(FthenM), "\n")
  IDstring = paste("\n\n", "Females then males")
  cat(file=fname, IDstring, "\n", append=T)
  write.table(file=fname, FthenM, sep=",", col.names=T, row.names=F, append=T)

  cat("Writing F only, dimensions:", dim(Fout), "\n")
  IDstring = paste("\n\n", "Females only")
  cat(file=fname, IDstring, "\n", append=T)
  write.table(file=fname, Fout, sep=",", col.names=T, row.names=F, append=T)

  cat("Writing M only, dimensions:", dim(Mout), "\n")
  IDstring = paste("\n\n",  "Males only")
  cat(file=fname, IDstring, "\n", append=T)
  write.table(file=fname, Mout, sep=",", col.names=T, row.names=F, append=T)

  cat("Writing combined sexes as females, dimensions:", dim(Uout), " \n")
  IDstring = paste("\n\n", "Sexes combined")
  cat(file=fname, IDstring, "\n", append=T)
  write.table(file=fname, Uout, sep=",", col.names=T, row.names=F, append=T)

  # Reset warnings

  options("warn" = oldwarn[[1]])

} # End function writeComps

##############################################################################
#
# getComps aggregates by length, by age, or by age-at-length according to the
# given stratification.
# 
# Note that the aggregation is of the Pdata$Final_Sample_Size value, which should be
# set to the desired expansion, e.g. Pdata$Final_Sample_Size = Pdata$Expansion_Factor_1
# or Pdata$Final_Sample_Size = Pdata$Expansion_Factor_1 * Pdata$Expansion_Factor_2
#
# The default stratification is by fleet, fishyr, and season.
# The lengthcm, age or both are appended depending on the "Comps" argument.
# The "strat" argument is prepended to this list.
#
##############################################################################

getComps = function( Pdata, strat=NULL, Comps="AAL" ) {

  # Set up stratification

  usualSuspects = c("fleet", "fishyr", "season")

  if (Comps == "LENGTH" | Comps == "LEN") {

    usualSuspects = c(usualSuspects, "lengthcm")

  } else if (Comps == "AGE") {

    usualSuspects = c(usualSuspects, "agemethod", "age")

  } else {

    usualSuspects = c(usualSuspects, "agemethod", "lengthcm", "age")

  } # End if-else-else

  strat = strat[!strat %in% usualSuspects]

  Cstrat = c(strat, usualSuspects)

  printIf( paste("Aggregating, stratification is by", paste(Cstrat, collapse=", ")) )

  # Used to get the number of SAMPLE_NOs per aggregation

  lenique = function(x) { return(length(unique(x))) }

  tmp = Pdata[Pdata$SEX == "M",]
  maleAgeComps = aggregate(tmp$Final_Sample_Size, tmp[,Cstrat], sum, na.rm=T)
  maleSamples = aggregate(tmp$FREQ, tmp[,Cstrat], sum, na.rm=T)
  maleTows = aggregate(tmp$SAMPLE_NO, tmp[,Cstrat], lenique)


  tmp = Pdata[Pdata$SEX == "F",]
  femaleAgeComps = aggregate(tmp$Final_Sample_Size, tmp[,Cstrat], sum, na.rm=T)
  femaleSamples = aggregate(tmp$FREQ, tmp[,Cstrat], sum, na.rm=T)
  femaleTows = aggregate(tmp$SAMPLE_NO, tmp[,Cstrat], lenique)


  tmp = Pdata[Pdata$SEX == "U",]
  unSexedAgeComps = aggregate(tmp$Final_Sample_Size, tmp[,Cstrat], sum, na.rm=T)
  unSamples = aggregate(tmp$FREQ, tmp[,Cstrat], sum, na.rm=T)
  unTows = aggregate(tmp$SAMPLE_NO, tmp[,Cstrat], lenique)

  names(maleAgeComps) = c(Cstrat, "male")
  names(femaleAgeComps) = c(Cstrat, "female")
  names(unSexedAgeComps) = c(Cstrat, "unsexed")

  names(maleSamples) = c(Cstrat, "msamps")
  names(femaleSamples) = c(Cstrat, "fsamps")
  names(unSamples) = c(Cstrat, "usamps")

  names(maleTows) = c(Cstrat, "mtows")
  names(femaleTows) = c(Cstrat, "ftows")
  names(unTows) = c(Cstrat, "utows")

  # Add tows and samples to the Comps

  maleAgeComps = merge(maleAgeComps, maleSamples, by=Cstrat, all=T)
  maleAgeComps = merge(maleAgeComps, maleTows, by=Cstrat, all=T)

  femaleAgeComps = merge(femaleAgeComps, femaleSamples, by=Cstrat, all=T)
  femaleAgeComps = merge(femaleAgeComps, femaleTows, by=Cstrat, all=T)

  unSexedAgeComps = merge(unSexedAgeComps, unSamples, by=Cstrat, all=T)
  unSexedAgeComps = merge(unSexedAgeComps, unTows, by=Cstrat, all=T)

  ageComps = merge(maleAgeComps, femaleAgeComps, by=Cstrat, all=T)
  ageComps = merge(ageComps, unSexedAgeComps, by=Cstrat, all=T)

  return(ageComps)

} # End function getComps


#  cat("\nExamining the gender-specific age comps,\n")

#  if (sum(out1.ageComps[,c("F0","M")]) > 0 ) {

#    cat("\n There are observations in the F0 and/or M0 columns.\n",
#        "Sum these into the first age bin\n")

#  } else {

#    cat("\n There are no observations in the F0/M0 columns.\n",
#        "These columns are being removed.\n")

#    out1.ageComps = out1.ageComps[,!names(ageComps %in% c("F0","M0"))]

#  } # End if-else

  # Do non-gender-specific comps

#  out2.ageComps = CommAges2SS3.gn(agecomps, ageBins, gender=0,
#                             nSamps="EnterNsamps", partition=2)

#  cat("\nExamining the gender-NON-specific age comps,\n")

#  if (sum(out1.ageComps[,"U0"]) > 0 ) {

#    cat("\n There are observations in the U0 column.\n",
#        "Sum these into the first age bin\n")

#  } else {

#    cat("\n There are no observations in the U0 column.\n",
#        "This columns is being removed.\n")

#    out1.ageComps = out1.ageComps[,!names(ageComps %in% c("F0","M0"))]

#  } # End if-else

#  return(out1.ageComps, out2.ageComps)

#} # End doAAL


##############################################################################
#
# doDiags creates diagnostic plots and summaries, writing them to a file
# in addition to plotting onscreen and console.
#
##############################################################################

doDiags = function( Pdata, fname=NULL ) {

  printIf( "Running diagnostics" )

  if ( is.null(fname) ) {

    # Set up filenames for txt, pdf

    species = sort(unique(Pdata$SPID))
    outfile = paste( "Diags.", species, ".txt", sep="" )
    pdffile = paste( "Diags.", species, ".pdf", sep="")

  } else {

    outfile = paste(fname, ".txt", sep="")
    pdffile = paste(fname, ".pdf", sep="")

  } # End ifelse

  cat( "Text will be written to", outfile, "\n" )
  cat( "Plots will be written to", pdffile, "\n" )

  # Open outfile for writing


  # Develop statistics of interest

  len = Pdata[!is.na(Pdata$FISH_LENGTH),]
  len$len = floor(len$FISH_LENGTH/10)
  len$depth_mid = (len$DEPTH_MIN+len$DEPTH_MAX)/2
  ltows = len[!duplicated(len$SAMPLE_NO),]

  meanLen.yr = tapply(len$len,list(len$SAMPLE_YEAR),mean)
  meanLen = tapply(len$len,list(len$SAMPLE_NO,len$SAMPLE_YEAR),mean)

  age = Pdata[!is.na(Pdata$FISH_AGE_YEARS_FINAL),]
  age$age = age$FISH_AGE_YEARS_FINAL
  atows = age[!duplicated(age$SAMPLE_NO),]
  meanAge = tapply(age$age,list(age$SAMPLE_NO,age$SAMPLE_YEAR),mean)

  # Print tables

  cat("Lengths for which FISH_LENGTH_TYPE is T:  ")
  print(len[len$FISH_LENGTH_TYPE=="T",])
  cat("\n\n")

  cat("Records per SAMPLE_YEAR\n\n")
  print(table(Pdata$SAMPLE_YEAR,useNA="ifany"))
  cat("\n\n")

  cat("SOURCE_AGID vs. SAMPLE_AGENCY\n")
  print(table(Pdata$SOURCE_AGID,Pdata$SAMPLE_AGENCY,useNA="ifany"))
  cat("\n\n")

  cat("FISH_LENGTH_TYPE\n")
  print(table(Pdata$FISH_LENGTH_TYPE,useNA="ifany"))
  cat("\n\n")

  cat("FISH_LENGTH\n")
  print(table(Pdata$FISH_LENGTH,useNA="ifany"))
  cat("\n\n")

  cat("GEAR vs GRID\n")
  print(table(len$GEAR,len$GRID))
  cat("\n\n")

  cat("FISH_LENGTH for lengthed fish\n")
  print(table(len$FISH_LENGTH_TYPE,useNA="ifany"))
  cat("\n\n")

  cat("SAMPLE_YEAR vs SOURCE_AGID for lengthed fish\n")
  print(table(len$SAMPLE_YEAR,len$SOURCE_AGID))
  cat("\n\n")
  
  cat("Difference between FISH_LENGTH and floor(FISH_LENGTH)\n")
  print(table(len$FISH_LENGTH-floor(len$FISH_LENGTH)))
  cat("\n\n")

  cat("Difference between FISH_LENGTH/10 and floor(FISH_LENGTH/10)\n")
  print(table(round(len$FISH_LENGTH/10-floor(len$FISH_LENGTH/10),1)))
  cat("\n\n")

  cat("DEPTH_AVG for lengthed fish\n")
  print(table(is.na(len$DEPTH_AVG)))
  cat("\n\n")

  cat("SAMPLE_YEAR vs. SOURCE_AGID for SAMPLE_NOs with lengthed fish\n")
  print(table(ltows$SAMPLE_YEAR,ltows$SOURCE_AGID))
  cat("\n\n")

  cat("DEPTH_AVG for SAMPLE_NOs with lengthed fish\n")
  print(table(is.na(ltows$DEPTH_AVG),useNA="ifany"))
  cat("\n\n")

  cat("Number of aged fish\n")
  print(nrow(age))
  cat("\n\n")

  cat("SAMPLE_YEAR vs. SOURCE_AGID for SAMPLE_NOs with aged fish\n")
  print(table(atows$SAMPLE_YEAR,atows$SOURCE_AGID))
  cat("\n\n")

  cat("SAMPLE_YEAR vs. SOURCE_AGID for aged fish\n")
  print(table(age$SAMPLE_YEAR,age$SOURCE_AGID))
  cat("\n\n")

  cat("age2 vs. age3 for aged fish\n")
  print(table(age$age2,age$age3,useNA="ifany"))
  cat("\n\n")

  cat("age1 vs. age2 for aged fish\n")
  print(table(age$age1,age$age2,useNA="ifany"))
  cat("\n\n")
  
  # Close output file


  # Plots

  # Fingers crossed, works the same for Mac and PC

  # Set up device for pdf

  while (! is.null(dev.list()) ) {

    dev.off()

  } # End while

  pdf(pdffile)

  par(mfrow=c(2,2))

  hist(len$len,nclass=30, xlab="", main="FISH_LENGTH")

  barplot(table(10*round(len$FISH_LENGTH/10-floor(len$FISH_LENGTH/10),1)),
          xlab="Difference in rounded and floored lengths")

  plot(len$FISH_LENGTH,len$FORK_LENGTH,pch=16, xlab="FISH_LENGTH", ylab="FORK_LENGTH")

  plot(len$DEPTH_AVG,len$depth_mid,xlim=c(0,400),ylim=c(0,400), xlab="DEPTH_AVG", ylab="Depth_mid")
  abline(a=0,b=1)

  hist(ltows$DEPTH_AVG, xlab="", main="DEPTH_AVG")

  hist(age$age,nclass=30, xlab="", main="Age")

  par(mfrow=c(2,1))
  boxplot(as.list(as.data.frame(meanLen)),varwidth=T,main="Mean length")
  boxplot(as.list(as.data.frame(meanAge)),varwidth=T,main="Mean age")

  dev.off()

} # End doDiags
  
#
# That's All, Folks!
#
##############################################################################
