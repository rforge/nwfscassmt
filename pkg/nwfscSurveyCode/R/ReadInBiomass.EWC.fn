###############################################################################################
#Reads in the triennial survey data and filters the data into what is necessary
#It reads in catch and station data and makes sure only the species necessary are kept
#### may want to keep NA (blank in Excel) to select the zero tows
#removeCAN is a flag if you want tows in Canadian waters removed
#   need the file called foreign_hauls.csv
#
#Written by Allan Hicks
#Necessary column names
#   SPECIES_CODE
#   WEIGHT in kg
#   DISTANCE_FISHED is in km
#   NET_WIDTH is in m
###############################################################################
ReadInBiomass.EWC.fn <- function(dataFile,directory,species=c(NA),removeCAN=T,verbose=F) {
    dat <- read.csv(paste(directory,dataFile,sep="\\"))
    totRows <- nrow(dat)
    dat <- dat[dat$SPECIES_CODE %in% species,]
    if(verbose) {cat(nrow(dat),"rows kept out of",totRows,"after filtering by species\n")}
    totRows <- nrow(dat)
    
    if(removeCAN) {
        foreignHauls <- read.csv(paste(directory,"foreign_hauls.csv",sep="\\"))
        foreignInd <- !(dat$HAULJOIN %in% foreignHauls$HAULJOIN)
        dat <- dat[foreignInd,]
        if(verbose) {cat(sum(foreignInd),"rows kept (or",sum(!foreignInd),"removed) out of",totRows,"after removing foreign hauls\n")}
    }            
    
    #calculate the density (kg/km^2) using net width and distance fished
    dat$areaFished <- dat$DISTANCE_FISHED*(dat$NET_WIDTH/1000)
    tmp <- sum(is.na(dat$areaFished))
    if(tmp>0 | verbose) {cat("There are",tmp,"instances where area swept could not be calculated due to missing values.\n")}
    dat$kgPerKm2 <- dat$WEIGHT/dat$areaFished
    dat$kgPerKm2[is.na(dat$kgPerKm2)&!is.na(dat$areaFished)] <- 0                                  #the tows with no observation of the species
    return(dat)
}
