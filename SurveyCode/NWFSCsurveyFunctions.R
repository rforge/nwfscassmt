#NWFSC survey function
#Written by Allan Hicks
#allan.hicks@noaa.gov
#206-302-2435

readDataFromExcel.fn <- function(file,sheet,headerRow) {
    #reads in the data from an Excel sheet
    #the column names of the data are on row number headerRow
    #this function is meant to be used by more specific functions below
    #written in 2009 by Allan Hicks
    #note that in later, 64-bit versions of R, the RODBC package is not working
    require(RODBC)
    channel <- odbcConnectExcel(file)
    info <- sqlFetch(channel,sheet,as.is=T,max=headerRow-2)        #read in info at top of sheet. It always puts the first row as header
    xx <- sqlFetchMore(channel,as.is=T,colnames=T)
    close(channel)
    print(as.data.frame(info[,1]))
    return(xx)
}
    
#readInExcelRawBiomassData.fn <- function(file,sheet="HaulCatchWt&Effort",headerRow=9) {
#I wanted to read in the raw data shett, but it won't read in the numbers.
#it seems to be caused by the header, because when the info at the top is removed, it works fine.
#I have no idea why the Biomass abundance works
#}


readInExcelBiomass.fn <- function(file,sheet="BiomassAbundance",headerRow=6,colNames=NA) {
    #Reads in the stratum biomasses from the Excel file provided by Beth (BiomassAbundance Sheet)
    #headerRow is the row number of the column where the data start
    #it doesn't read in the column names correctly, so I put in simplified names. Make sure that these match what is in your Excel spreadsheet
    #written by Allan Hicks, 2009

    xx <- readDataFromExcel.fn(file,sheet,headerRow)
    if(is.null(colNames[1])) {
        nombres <- c("Species","ScientificName","SpeciesCode","Year","Project","StrataAreaVersion","AreaSetId","AreaName","SouthernLatitude","NorthernLatitude","DepthStrataSet","MinStratumDepth","MaxStratumDepth","StratumArea","Biomass","Abundance","CpueWeightVar","CpueCountVar","BiomassVar","CV","N","Nbio","Npos","NbioPos")
        names(xx) <- nombres
        cat("NOTE: column names have been modified from the Excel Sheet. You may want to verify that they match.\n")
    }
    if(!is.na(colNames[1])) {
        names(xx) <- colNames
        return(xx)
    }
    return(xx)
}
#PetBio <- readInExcelBiomass.fn("C:\\NOAA2009\\Petrale\\Data\\Survey\\NWFSCsurvey\\FisheryIndices2009_Petrale.xls")

readInBiomass.fn <- function(filename,headerRow=1,sep=",",
                                colNames = c("Species","ScientificName","SpeciesCode","Year","Project","StrataAreaVersion","AreaSetId","AreaName","SouthernLatitude","NorthernLatitude","DepthStrataSet","MinStratumDepth","MaxStratumDepth","StratumArea","Biomass","Abundance","CpueWeightVar","CpueCountVar","BiomassVar","CV","N","Nbio","Npos","NbioPos")                    )
{
    #Reads in the stratum biomasses from the text file saved from the Excel worksheet
    #headerRow is the line number where the data start (use to skip over header lines)
    #you can keep your column names or use the simplified names I provide. Make sure that these match what is in your Excel spreadsheet or make sense!
    #written by Allan Hicks, 2011

    xx <- read.table(filename,skip=headerRow-1,sep=sep,header=T)
    if(length(colNames) == ncol(xx)) {
        names(xx) <- colNames
        cat("NOTE: column names have been modified from the csv file. You may want to verify that they match.\n")
    }
    return(xx)
}






SS3Biomass.fn <- function(bio,fleet="EnterFleet",season=1,outputMedian=T) {
    #This function outputs a dataframe in the abundance index format needed for SS3
    #it calculates the standard error of log(B) by first finding the cv of B
    #It can output the median biomass or the mean biomass. Median biomass is bias corrected so that log(Bmedian)=mean(log(B))
    xxx <- split(bio[,c("AreaName","MinStratumDepth","MaxStratumDepth","Biomass","BiomassVar")],bio$Year)
    years <- names(xxx)
    bio <- unlist(lapply(xxx,function(x){sum(x$Biomass,na.rm=T)}))
    variance <- unlist(lapply(xxx,function(x){sum(x$BiomassVar,na.rm=T)}))
    cv <- sqrt(variance)/bio
    seLogB <- sqrt(log(cv^2+1))
    if(outputMedian==T) {
        med <- bio*exp(-0.5*seLogB^2)
        out <- data.frame(Year=years,Season=season,Fleet=fleet,Value=med,seLogB=seLogB)
    }
    else {
        out <- data.frame(Year=years,Season=season,Fleet=fleet,Value=bio,seLogB=seLogB)
    }
    row.names(out) <- NULL
    return(out)
}    

strataBiomass.fn <- function(bio,theStrata="Year",fleet="EnterFleet",season=1,outputMedian=T) {
    #This function outputs a dataframe in the abundance index format needed for SS3
    #it calculates the standard error of log(B) by first finding the cv of B
    #It can output the median biomass or the mean biomass. Median biomass is bias corrected so that log(Bmedian)=mean(log(B))
    #print(bio[,theStrata])
    xxx <- split(bio[,c(theStrata,"Biomass","BiomassVar")],bio[,theStrata])
    strat <- names(xxx)
    bio <- unlist(lapply(xxx,function(x){sum(x$Biomass,na.rm=T)}))
    variance <- unlist(lapply(xxx,function(x){sum(x$BiomassVar,na.rm=T)}))
    cv <- sqrt(variance)/bio
    seLogB <- sqrt(log(cv^2+1))
    if(outputMedian==T) {
        med <- bio*exp(-0.5*seLogB^2)
        out <- data.frame(Stratum=strat,Season=season,Fleet=fleet,Value=med,seLogB=seLogB)
    }
    else {
        out <- data.frame(Stratum=strat,Season=season,Fleet=fleet,Value=bio,seLogB=seLogB)
    }
    row.names(out) <- NULL
    return(out)
}    

GetTotalBiomassExcel.fn <- function(file,sheet="BiomassAbundance",headerRow,fleet="EnterFleet",season=1,outputMedian=T) {
    #a wrapper for to output the biomass in SS3 format, reading directly from Excel
    bio <- readInExcelBiomass.fn(file,sheet,headerRow)
    SS3Biomass.fn(bio,fleet,season)
}
GetTotalBiomass.fn <- function(file,headerRow,fleet="EnterFleet",season=1,outputMedian=T) {
    #a wrapper for to output the biomass in SS3 format, reading in from a csv file
    bio <- readInBiomass.fn(file,headerRow)
    SS3Biomass.fn(bio,fleet,season)
}
    
#GetTotalBiomass.fn("C:\\NOAA2009\\Petrale\\Data\\Survey\\NWFSCsurvey\\FisheryIndices2009_Petrale.xls",headerRow=6)


readInExcelLengthComps.fn <- function(file,sheet="LengthComps",headerRow=7) {
    #Reads in the stratum numbers from the Excel file provided by Beth (LengthComps Sheet)
    #headerRow is the row number of the column where the data start
    #it doesn't read in the column names correctly, so I put in simplified names. Make sure that these match what is in your Excel spreadsheet
    #written by Allan Hicks, 2009
    nombres <- c("SpeciesCode","ScientificName","Species","Year","Project","AreaSetIdentifier","AreaName","DepthStrataSet","MinStratumDepth","MaxStratumDepth","Length","NumF","NumM","NumUnsexed")

    xx <- readDataFromExcel.fn(file,sheet,headerRow)
    names(xx) <- nombres
    cat("\nNOTE: column names have been modified from the Excel Sheet. You may want to verify that they match.\n")
    return(xx)
}
#PetLen <- readInExcelLengthComps.fn("C:\\NOAA2009\\Petrale\\Data\\Survey\\NWFSCsurvey\\FisheryIndices2009_Petrale.xls")


readInLengthComps.fn <- function(file,headerRow=7,sep=",",
                                    colNames=c("SpeciesCode","ScientificName","Species","Year","Project","AreaSetIdentifier","AreaName","DepthStrataSet","MinStratumDepth","MaxStratumDepth","Length","NumF","NumM","NumUnsexed")
) {
    #Reads in the stratum numbers from a csv file (usually converted from Excel file provided by Beth (LengthComps Sheet))
    #headerRow is the row number of the column where the data start
    #it doesn't read in the column names correctly, so I put in simplified names. Make sure that these match what is in your Excel spreadsheet
    #written by Allan Hicks, 2011

    xx <- read.table(file,skip=headerRow-1,sep=sep,header=T)
    if(length(colNames) == ncol(xx)) {
        names(xx) <- colNames
        cat("NOTE: column names have been modified from the csv file. You may want to verify that they match.\n")
    }
    return(xx)
}


SS3LF.fn <- function(len,lgthBins=1,gender=3,nSamps="EnterNsamps",fleet="EnterFleet",season=1,partition=0,NAs2zero=T,sexRatioUnsexed=NA,maxSizeUnsexed=NA) {
    #calculates proportion at length and reformats into SS# format
    #Gender=0: sexes combined and entered in female placeholders. Male values ignored
    #Gender=1: females only. Male values ignored
    #Gender=2: males only. Female values ignored.
    #Gender=3: both sexes. Proportions over males and females sum to 1. See below for how unsexed are treated
    #lgthBins is either the interval between length bins or the actual length bins
    #note that 0 and Inf are tacked on the ends to account for lengths outside the interval
    #The largest length bin includes all lengths greater
    #NOTE: The length bin called F0 or M0 is retained to show proportion of lengths smaller than smallest bin
    #      You will want to likely add this to your first length bin and delete this before putting in SS3, or
    #       start the lgthBins argument at the 2nd length bin and F0 will be all fish smaller (hence the first length bin)
    #The NAs2zero determines if NA values will be changed to 0.0.
    #     You will want this to be T when inputting into SS3, but may want it to be false when plotting
    #     Note that an NA means there were no males and females recorded. Therefore you may get a zero if only one sex was recorded
    #Sex Ratio for unsexed
    #     sexRatioUnsexed is the proportion used to assign unsexed fish to females
    #       -- if a single number, it is only used for sizes at or below maxSizeUnsexed, and sex ratio for size bins above that are calculated using the 
    #           number of males and females observed in the same size class (or one lower if none available in the same bin)
    #       -- if a vector (maybe to be implemented), it must be the same length as the lgthBins and indicates the proportion of unsexed assigned to females for each length bin
    #       -- if it is NA, unsexed fish are omitted WHEN GENDER=3 (THIS IS THE DEFAULT)
    #     maxSizeUnsexed determines the maximum size at which the sexRatioUnsexed is applied. If sexRatioUnsexed is a vector, this is ignored
    
    if(length(lgthBins)==1) {
        Lengths <- c(0,seq(floor(min(len$Length)),ceiling(max(len$Length)),lgthBins),Inf)
    }
    else{
        Lengths <- c(0,lgthBins,Inf)        #put 0 and Inf on ends because all.inside=T in findInterval below. Treats these as minus and plus groups
    }

    len$allLs <- Lengths[findInterval(len$Length,Lengths,all.inside=T)]
    #print(table(len$allLs))
    
    if(length(sexRatioUnsexed)==1 & !is.na(sexRatioUnsexed)) {
        len$sexRatio <- len$NumF/(len$NumF+len$NumM)
        len$sexRatio[len$Length <= maxSizeUnsexed] <- sexRatioUnsexed
        #now fill in any missing ratios with ratios of that bin from other years and strata (can probably be done more efficiently)
        noRatio <- which(is.na(len$sexRatio))
        if(length(noRatio)>0) cat("\nThese are sex ratios that were filled in using observations from the same lengths from different strata and years\n")
        for(i in noRatio) {
            inds <- len$allLs==len$allLs[i]
            tmpF <- sum(len$NumF[inds])
            tmpM <- sum(len$NumM[inds])
            len$sexRatio[i] <- tmpF/(tmpF+tmpM)
            print(len[i,c("Length","allLs","NumF","NumM","sexRatio")])
        }

        noRatio <- which(is.na(len$sexRatio))
        if(length(noRatio)>0) cat("\nThese are sex ratios that were filled in using observations from nearby lengths\n")
        for(i in noRatio) {
            nearLens <- Lengths[c(which(Lengths==len$allLs[i])-1,which(Lengths==len$allLs[i])+1)]
            inds <- len$allLs %in% nearLens
            tmpF <- sum(len$NumF[inds])
            tmpM <- sum(len$NumM[inds])
            len$sexRatio[i] <- tmpF/(tmpF+tmpM)
            print(len[i,c("Length","allLs","NumF","NumM","sexRatio")])
        }
        noRatio <- which(is.na(len$sexRatio))
        if(length(noRatio)>0) cat("Some sex ratios were left unknown and omitted\n\n")
        if(length(noRatio)==0) cat("Done filling in sex ratios\n\n")

        tmpFemUnsex <- round(len$sexRatio*len$NumUnsexed)
        tmpMaleUnsex <- len$NumUnsexed - tmpFemUnsex
        len$NumF <- len$NumF + tmpFemUnsex
        len$NumM <- len$NumM + tmpMaleUnsex
        #print(unique(round(len$sexRatio,1)))
    }
    #if(length(sexRatioUnsexed) > 1)
    #    if(length(sexRatioUnsexed)!=(length(Lengths)-2)) stop("sexRatioUnsexed must be a single number or the same length as the length bins")   
    #}

    xx <- split(len[,c("Length","NumF","NumM","NumUnsexed")],len$Year)
    years <- names(xx)
    
    year.fn <- function(x,Lengths) {
        allLs <- Lengths[findInterval(x$Length,Lengths,all.inside=T)]    #finds the interval that the length falls in and floors it (so 23.2 would be in 23 if 23 was a level in Lengths, all.inside puts maximum age group into N-1 group, thus I padded with Inf.)
        totalU <- tapply(x$NumUnsexed,allLs,sum,na.rm=T)
        totalF <- tapply(x$NumF,allLs,sum,na.rm=T)
        totalM <- tapply(x$NumM,allLs,sum,na.rm=T)
        out <- data.frame(Length=Lengths,numF=rep(NA,length(Lengths)),numM=rep(NA,length(Lengths)),numU=rep(NA,length(Lengths)))
        row.names(out) <- out$Length
        out[names(totalF),"numF"] <- totalF
        out[names(totalM),"numM"] <- totalM
        out[names(totalU),"numU"] <- totalU
        out[names(totalF),"propF"] <- 100*totalF/(sum(totalF,na.rm=T)+sum(totalM,na.rm=T))
        out[names(totalM),"propM"] <- 100*totalM/(sum(totalF,na.rm=T)+sum(totalM,na.rm=T))
        out[names(totalU),"propU"] <- 100*(totalF+totalM)/sum(totalF+totalM,na.rm=T)            #unsexed have been added in above
        out <- out[-nrow(out),]   #remove last row because Inf is always NA due to inside.all=T
        return(out)
    }
    L.year <- lapply(xx,year.fn,Lengths=Lengths)

    #output SS3 format specific to the gender choice
    lgths <- as.character(L.year[[1]]$Length)
    if(gender==0) {
        Ls <- unlist(lapply(L.year,function(x){x$propU}))
        if(NAs2zero){Ls[is.na(Ls)] <- 0}
        Ls <- matrix(Ls,nrow=length(L.year),byrow=T,
            dimnames=list(NULL,paste(rep("U",length(lgths)),lgths,sep="")))
        out <- data.frame(year=as.numeric(names(L.year)),Season=season,Fleet=fleet,gender=gender,partition=partition,nSamps=nSamps,Ls,Ls)
    }
    if(gender==3) {
        Ls <- unlist(lapply(L.year,function(x){c(x$propF,x$propM)}))
        if(NAs2zero){Ls[is.na(Ls)] <- 0}
        Ls <- matrix(Ls,nrow=length(L.year),byrow=T,
            dimnames=list(NULL,paste(c(rep("F",length(lgths)),rep("M",length(lgths))),lgths,sep="")))
        out <- data.frame(year=as.numeric(names(L.year)),Season=season,Fleet=fleet,gender=gender,partition=partition,nSamps=nSamps,Ls)
    }

    cat("\nNOTE: You may need to add the column called F0 and/or M0 to your first length bin\n\tand delete that column.\n\tThese are the proportion of lengths smaller than the first length bin\n\n")
    return(out)
}    
#tmp <- SS3LF.fn(PetLen,lgthBins=2,gender=3)

GetLFs.fn <- function(file,headerRow,lgthBins=1,gender=3,nSamps="EnterNsamps",fleet="EnterFleet",season=1,partition=0,NAs2zero=T,sep=",",sexRatioUnsexed=NA,maxSizeUnsexed=NA) {
    len <- readInLengthComps.fn(file,headerRow=headerRow,sep=sep)
    SS3LF.fn(len,lgthBins=lgthBins,gender=gender,nSamps=nSamps,fleet=fleet,season=season,partition=partition,NAs2zero=NAs2zero,sexRatioUnsexed=NA,maxSizeUnsexed=NA)
}

GetLFsFromExcel.fn <- function(file,sheet="LengthComps",headerRow,lgthBins=1,gender=3,nSamps="EnterNsamps",fleet="EnterFleet",season=1,partition=0,NAs2zero=T,sexRatioUnsexed=NA,maxSizeUnsexed=NA) {
    len <- readInExcelLengthComps.fn(file,sheet=sheet,headerRow=headerRow)
    SS3LF.fn(len,lgthBins=lgthBins,gender=gender,nSamps=nSamps,fleet=fleet,season=season,partition=partition,NAs2zero,sexRatioUnsexed=NA,maxSizeUnsexed=NA)
}
#GetLFsFromExcel.fn("C:\\NOAA2009\\Petrale\\Data\\Survey\\NWFSCsurvey\\FisheryIndices2009_Petrale.xls",headerRow=7)

plotSexRatio.fn <- function(len,fn=median,circleSize=0.1) {
    ratioF <- len$NumF/(len$NumF+len$NumM)
    yF <- lapply(split(ratioF,floor(len$Length)),fn)
    x <- names(split(ratioF,floor(len$Length)))
    nobs <- unlist(lapply(split(ratioF,floor(len$Length)),length))
    plot(x,yF,type="l",col="red")
    symbols(x,yF,circles=nobs,inches=circleSize,fg="red",bg=rgb(1,0,0,alpha=0.5),add=T)
}


SS3LFstrata.fn <- function(len,Strata="Year",lgthBins=1,gender=3,nSamps="EnterNsamps",fleet="EnterFleet",season=1,partition=0,NAs2zero=T) {
    #calculates proportion at length and reformats into SS# format
    #Gender=0: sexes combined and entered in female placeholders. Male values ignored
    #Gender=1: females only. Male values ignored
    #Gender=2: males only. Female values ignored.
    #Gender=3: both sexes. Proportions over males and females sum to 1.
    #lgthBins is either the interval between length bins or the actual length bins
    #note that 0 and Inf are tacked on the ends to account for lengths outside the interval
    #The largest length bin includes all lengths greater
    #NOTE: The length bin called F0 or M0 is retained to show proportion of lengths smaller than smallest bin
    #      You will want to delete this before putting in SS3
    #The NAs2zero determines if NA values will be changed to 0.0.
    #     You will want this to be T when inputting into SS3, but may want it to be false when plotting
    #     Note that an NA means there were no males and females recorded. Therefore you may get a zero if only one sex was recorded
    # NOTICE: UNSEXED CORRECTION NOT IMPLEMENTED
    xx <- split(len[,c("Length","NumF","NumM","NumUnsexed")],len[,Strata])
    years <- names(xx)
    if(length(lgthBins)==1) {
        Lengths <- c(0,seq(floor(min(len$Length)),ceiling(max(len$Length)),lgthBins),Inf)
    }
    else{
        Lengths <- c(0,lgthBins,Inf)        #put 0 and Inf on ends because all.inside=T in findInterval below. Treats these as minus and plus groups
    }

    year.fn <- function(x,Lengths) {
        allLs <- Lengths[findInterval(x$Length,Lengths,all.inside=T)]    #finds the interval that the length falls in and floors it (so 23.2 would be in 23 if 23 was a level in Lengths, all.inside puts maximum age group into N-1 group, thus I padded with Inf.)
        totalU <- tapply(x$NumUnsexed,allLs,sum,na.rm=T)
        totalF <- tapply(x$NumF,allLs,sum,na.rm=T)
        totalM <- tapply(x$NumM,allLs,sum,na.rm=T)
        out <- data.frame(Length=Lengths,numF=rep(NA,length(Lengths)),numM=rep(NA,length(Lengths)),numU=rep(NA,length(Lengths)))
        row.names(out) <- out$Length
        out[names(totalF),"numF"] <- totalF
        out[names(totalM),"numM"] <- totalM
        out[names(totalU),"numU"] <- totalU
        out[names(totalF),"propF"] <- 100*totalF/(sum(totalF,na.rm=T)+sum(totalM,na.rm=T))
        out[names(totalM),"propM"] <- 100*totalM/(sum(totalF,na.rm=T)+sum(totalM,na.rm=T))
        out[names(totalU),"propU"] <- 100*(totalU+totalF+totalM)/sum(totalU+totalF+totalM,na.rm=T)
        out <- out[-nrow(out),]   #remove last row because Inf is always NA due to inside.all=T
        return(out)
    }
    L.year <- lapply(xx,year.fn,Lengths=Lengths)

    #output SS3 format specific to the gender choice
    lgths <- as.character(L.year[[1]]$Length)
    if(gender==0) {
        Ls <- unlist(lapply(L.year,function(x){x$propU}))
        if(NAs2zero){Ls[is.na(Ls)] <- 0}
        Ls <- matrix(Ls,nrow=length(L.year),byrow=T,
            dimnames=list(NULL,paste(rep("U",length(lgths)),lgths,sep="")))
        out <- data.frame(year=names(L.year),Season=season,Fleet=fleet,gender=gender,partition=partition,nSamps=nSamps,Ls,Ls)
    }
    if(gender==3) {
        Ls <- unlist(lapply(L.year,function(x){c(x$propF,x$propM)}))
        if(NAs2zero){Ls[is.na(Ls)] <- 0}
        Ls <- matrix(Ls,nrow=length(L.year),byrow=T,
            dimnames=list(NULL,paste(c(rep("F",length(lgths)),rep("M",length(lgths))),lgths,sep="")))
        out <- data.frame(year=names(L.year),Season=season,Fleet=fleet,gender=gender,partition=partition,nSamps=nSamps,Ls)
    }

    cat("\nNOTE: You may need to delete a column called F0 and/or M0. These are the proportion of lengths smaller than the first length bin\n\n")
    return(out)
}    
#tmp <- SS3LF.fn(PetLen,lgthBins=2,gender=3)


readInExcelAgeComps.fn <- function(file,sheet="AgeComps",headerRow=7) {
    #Reads in the stratum numbers from the Excel file provided by Beth (AgeComps Sheet)
    #headerRow is the row number of the column where the data start
    #it doesn't read in the column names correctly, so I put in simplified names. Make sure that these match what is in your Excel spreadsheet
    #written by Allan Hicks, 3/21/09
    nombres <- c("SpeciesCode","Species","Year","Project","AreaName","MinStratumDepth","MaxStratumDepth","Length","Age","NumF","NumM","NumUnsexed","LengthedAgeTally","AgeTallyF","AgeTallyM","AgeTallyU")

    xx <- readDataFromExcel.fn(file,sheet,headerRow)
    names(xx) <- nombres
    cat("\nNOTE: column names have been modified from the Excel Sheet. You may want to verify that they match.\n\n")
    return(xx)
}
#PetAge <- readInExcelAgeComps.fn("C:\\NOAA2009\\Petrale\\Data\\Survey\\NWFSCsurvey\\FisheryIndices2009_Petrale.xls")

readInAgeComps.fn <- function(file,headerRow=7,sep=",",
    colNames = c("SpeciesCode","Species","Year","Project","AreaName","MinStratumDepth","MaxStratumDepth","Length","Age","NumF","NumM","NumUnsexed","LengthedAgeTally","AgeTallyF","AgeTallyM","AgeTallyU")
) {
    #Reads in the stratum numbers from the Excel file provided by Beth (AgeComps Sheet)
    #headerRow is the row number of the column where the data start
    #it doesn't read in the column names correctly, so I put in simplified names. Make sure that these match what is in your Excel spreadsheet
    #written by Allan Hicks, 3/21/09

    xx <- read.table(file,skip=headerRow-1,sep=sep,header=T)
    if(length(colNames) == ncol(xx)) {
        names(xx) <- colNames
        cat("NOTE: column names have been modified from the csv file. You may want to verify that they match.\n")
    }
    return(xx)
}

SS3AgeAtLen.fn <- function(ages,lgthBins=1,ageBins=1,fleet="EnterFleet",season=1,partition=0,ageerr="EnterAgeErr",raw=T,sexRatioUnsexed=NA,maxSizeUnsexed=NA) {
    #calculates proportion of age at length and reformats into SS3 format
    #Uses raw numbers at length, assuming that is a random sample conditioned on length and sex.
    #To use expanded numbers, set raw=F
    #Only gender codes 1 and 2 and puts males and females on separate lines because the age@L is conditioned on sex (a sample of females of length 25cm, for example)
    #Gender=1: females only. Male values ignored
    #Gender=2: males only. Female values ignored.
    #lgthBins is either the interval between length bins or the actual length bins
    #note that 0 and Inf are tacked on the ends to account for lengths and ages outside the interval. You may want to add these in to first and last bin.
    years <- sort(unique(ages$Year))
    if(length(lgthBins)==1) {
        Lengths <- c(0,seq(floor(min(ages$Length)),ceiling(max(ages$Length)),lgthBins),Inf)
    }
    else{
        Lengths <- c(0,lgthBins,Inf)        #put 0 and Inf on ends because all.inside=T in findInterval below. Treats these as minus and plus groups
    }
    if(length(ageBins)==1) {
        Ages <- c(0,seq(floor(min(ages$Age)),ceiling(max(ages$Age)),ageBins),Inf)
    }
    else{
        Ages <- c(0,ageBins,Inf)        #put 0 and Inf on ends because all.inside=T in findInterval below. Treats these as minus and plus groups
    }
    allLs <- Lengths[findInterval(ages$Length,Lengths,all.inside=T)]   #finds the interval that the length falls in and floors it (so 23.2 would be in 23 if 23 was a level in Lengths, all.inside puts maximum age group into N-1 group, thus I padded with Inf.)

    if(length(sexRatioUnsexed)==1 & !is.na(sexRatioUnsexed)) {
        ages$allLs <- allLs
        ages$sexRatio <- ages$NumF/(ages$NumF+ages$NumM)
        ages$sexRatio[ages$Length <= maxSizeUnsexed] <- sexRatioUnsexed
        #now fill in any missing ratios with ratios of that bin from other years and strata (can probably be done more efficiently)
        noRatio <- which(is.na(ages$sexRatio))
        if(length(noRatio)>0) cat("\nThese are sex ratios that were filled in using observations from the same lengths from different strata and years\n")
        for(i in noRatio) {
            inds <- allLs==allLs[i]
            tmpF <- sum(ages$NumF[inds])
            tmpM <- sum(ages$NumM[inds])
            ages$sexRatio[i] <- tmpF/(tmpF+tmpM)
            print(ages[i,c("Length","allLs","Age","NumF","NumM","sexRatio")])
        }

        noRatio <- which(is.na(ages$sexRatio))
        if(length(noRatio)>0) cat("\nThese are sex ratios that were filled in using observations from nearby lengths\n")
        for(i in noRatio) {
            nearLens <- Lengths[c(which(Lengths==allLs[i])-1,which(Lengths==allLs[i])+1)]
            inds <- ages$allLs %in% nearLens
            tmpF <- sum(ages$NumF[inds])
            tmpM <- sum(ages$NumM[inds])
            ages$sexRatio[i] <- tmpF/(tmpF+tmpM)
            print(ages[i,c("Length","allLs","Age","NumF","NumM","sexRatio")])
        }
        noRatio <- which(is.na(ages$sexRatio))
        if(length(noRatio)>0) cat("Some sex ratios were left unknown and omitted\n\n")
        if(length(noRatio)==0) cat("Done filling in sex ratios\n\n")

        tmpFemUnsex <- round(ages$sexRatio*ages$NumUnsexed)
        tmpMaleUnsex <- ages$NumUnsexed - tmpFemUnsex
        ages$NumF <- ages$NumF + tmpFemUnsex
        ages$NumM <- ages$NumM + tmpMaleUnsex
        print(unique(round(ages$sexRatio,1)))
    }


    if(raw){xx <- split(ages[,c("Year","Length","Age","AgeTallyF","AgeTallyM")],paste(ages$Year,allLs))}
    if(!raw){
        ages[,c("AgeTallyF","AgeTallyM")] <- ages[,c("NumF","NumM")]  #use the expanded numbers
        xx <- split(ages[,c("Year","Length","Age","AgeTallyF","AgeTallyM")],paste(ages$Year,allLs))
    }

    bin.fn <- function(x,Ages) {
        allAs <- Ages[findInterval(x$Age,Ages,all.inside=T)]
        totalF <- tapply(x$AgeTallyF,allAs,sum,na.rm=T)
        totalM <- tapply(x$AgeTallyM,allAs,sum,na.rm=T)
        out <- data.frame(Age=Ages,numF=rep(NA,length(Ages)),numM=rep(NA,length(Ages)))
        row.names(out) <- out$Age
        out[names(totalF),"numF"] <- totalF
        out[names(totalM),"numM"] <- totalM
        out[names(totalF),"propF"] <- 100*totalF/sum(totalF,na.rm=T)
        out[names(totalM),"propM"] <- 100*totalM/sum(totalM,na.rm=T)
        out <- out[-nrow(out),]   #remove last row because Inf and always NA due to inside.all=T (but needed in findInterval)
        return(out)
    }
    A.bin <- lapply(xx,bin.fn,Ages=Ages)
    
    Nobs.fn <- function(x) {
        nF <- sum(x$AgeTallyF,na.rm=T)
        nM <- sum(x$AgeTallyM,na.rm=T)
        out <- c(nF,nM)
        names(out) <- c("nF","nM")
        return(out)
    }
    nobs <- lapply(xx,Nobs.fn)

    #output SS3 format with gender on separate lines
    ages <- as.character(A.bin[[1]]$Age)
    
    #gender=1 (females only, males ignored)
    nsF <- unlist(lapply(nobs,function(x){x["nF"]}))
    nsM <- unlist(lapply(nobs,function(x){x["nM"]}))
    AsF <- unlist(lapply(A.bin,function(x){x$propF}))
    AsM <- unlist(lapply(A.bin,function(x){x$propM}))
    AsF[is.na(AsF)] <- 0
    AsM[is.na(AsM)] <- 0
    AsF <- matrix(AsF,nrow=length(A.bin),byrow=T,
          dimnames=list(NULL,paste(rep("F",length(ages)),ages,sep="")))
    AsF[,"F1"] <- AsF[,"F0"]+AsF[,"F1"]     #add in all ages before the minimum age to the first age bin
    numFzero <- sum(AsF[,"F0"])
    AsF <- AsF[,-match("F0",dimnames(AsF)[[2]])]        #remove F0 column
    AsM <- matrix(AsM,nrow=length(A.bin),byrow=T,
          dimnames=list(NULL,paste(rep("M",length(ages)),ages,sep="")))
    AsM[,"M1"] <- AsM[,"M0"]+AsM[,"M1"]     #add in all ages before the minimum age to the first age bin
    numMzero <- sum(AsM[,"M0"])
    AsM <- AsM[,-match("M0",dimnames(AsM)[[2]])]

    outF <- data.frame(year=as.numeric(substring(names(A.bin),1,4)),Season=season,Fleet=fleet,gender=1,partition=partition,ageErr=ageerr,
                          LbinLo=as.numeric(substring(names(A.bin),6)),LbinHi=as.numeric(substring(names(A.bin),6)),nSamps=nsF,AsF,AsF)
    outM <- data.frame(year=as.numeric(substring(names(A.bin),1,4)),Season=season,Fleet=fleet,gender=2,partition=partition,ageErr=ageerr,
                          LbinLo=as.numeric(substring(names(A.bin),6)),LbinHi=as.numeric(substring(names(A.bin),6)),nSamps=nsM,AsM,AsM)
    indZero <- apply(outF[,-c(1:9)],1,sum)==0
    outF <- outF[!indZero,]   #remove any rows that have no female observations (they may be there because of male obs)
    indZero <- apply(outM[,-c(1:9)],1,sum)==0
    outM <- outM[!indZero,]   #remove any rows that have no male observations (they may be there because of female obs)
    rownames(outF) <- paste("F",1:nrow(outF),sep="")
    rownames(outM) <- paste("M",1:nrow(outM),sep="")

    cat("There are",numFzero,"females in the age 0 to age",ages[2],"that were added into the first age bin\n")
    cat("There are",numMzero,"males in the age 0 to age",ages[2],"that were added into the first age bin\n")
    cat("The number of fish in each sample were input into the nSamps column\nUse Beth's Excel file for the number of tows")
    return(list(female=outF,male=outM))
}
#tmp <- SS3AgeAtLen.fn(PetAge,lgthBins=2,ageBins=1:17)

GetAgesExcel.fn <- function(file,sheet="AgeComps",headerRow,lgthBins=1,ageBins=1,nSamps="EnterNsamps",fleet="EnterFleet",season=1,partition=0,ageerr="EnterAgeErr",raw=T) {
    age <- readInExcelAgeComps.fn(file,sheet=sheet,headerRow=headerRow)
    SS3AgeAtLen.fn(age,lgthBins=lgthBins,ageBins=ageBins,fleet=fleet,season=season,partition=partition,ageerr=ageerr,raw=raw)
}
#GetAges.fn("C:\\NOAA2009\\Petrale\\Data\\Survey\\NWFSCsurvey\\FisheryIndices2009_Petrale.xls",headerRow=7,lgthBins=seq(12,62,2),ageBins=1:17)

GetAges.fn <- function(file,headerRow,lgthBins=1,ageBins=1,nSamps="EnterNsamps",fleet="EnterFleet",season=1,partition=0,ageerr="EnterAgeErr",raw=T,sep=",") {
    age <- readInAgeComps.fn(file,headerRow=headerRow,sep=sep)
    SS3AgeAtLen.fn(age,lgthBins=lgthBins,ageBins=ageBins,fleet=fleet,season=season,partition=partition,ageerr=ageerr,raw=raw)
}




plotBio.fn <- function(bio,CI=0.95,scalar=1e6,gap=0.03,ylab="Biomass ('000 mt)",xlab="Year",ylim=NULL,...) {
    #Plots the biomass with confidence intervals
    #uses data in the format of SS3, so you can use the GetTotalBiomass.fn or your own data.frame
    #scalar is simply the divisor for the biomass
    #gap is a value that introduces a slight gap between the point estimate and the start of the line for the CI
    # careful because a gap too large will invert the CI, making it look huge. You should know when this happens
    y <- as.numeric(as.character(bio$Value/scalar))
    x <- as.numeric(as.character(bio$Year))
    se <- as.numeric(as.character(bio$seLogB))
    logB <- log(bio$Value)
    ci <- exp(rbind(c(logB+qnorm(1-(1-CI)/2)*se),c(logB-qnorm(1-(1-CI)/2)*se)))/scalar
    if(is.null(ylim)) {
        ylim <- c(0,1.05*max(ci))
    }
    
    gap <- gap*max(y)
    plot(x,y,,ylab=ylab,xlab=xlab,ylim=ylim,...)
    segments(x,y+gap,x,ci[1,])
    segments(x,y-gap,x,ci[2,])
}
#bio <- GetTotalBiomass.fn("C:\\NOAA2009\\Petrale\\Data\\Survey\\NWFSCsurvey\\FisheryIndices2009_Petrale.xls",headerRow=6)
#plotBio.fn(bio,pch=16)

plotBioStrata.fn <- function(bio,CI=0.95,scalar=1e6,gap=0.03,ylab="Biomass ('000 mt)",xlab="Stratum",ylim=NULL,...) {
    #Plots the biomass with confidence intervals
    #uses data in the format of SS3, so you can use the GetTotalBiomass.fn or your own data.frame
    #scalar is simply the divisor for the biomass
    #gap is a value that introduces a slight gap between the point estimate and the start of the line for the CI
    # careful because a gap too large will invert the CI, making it look huge. You should know when this happens
    y <- as.numeric(as.character(bio$Value/scalar))
    x <- 1:nrow(bio)
    se <- as.numeric(as.character(bio$seLogB))
    logB <- log(bio$Value)
    ci <- exp(rbind(c(logB+qnorm(1-(1-CI)/2)*se),c(logB-qnorm(1-(1-CI)/2)*se)))/scalar
    if(is.null(ylim)) {
        ylim <- c(0,1.05*max(ci))
    }
    
    gap <- gap*max(y)
    plot(x,y,ylab=ylab,xlab=xlab,ylim=ylim,xaxt="n",...)
    axis(1,at=x,label=as.character(bio$Stratum),padj=0,las=2,cex.axis=0.5)
    segments(x,y+gap,x,ci[1,])
    segments(x,y-gap,x,ci[2,])
}

plotFreqData.fn <- function(dat,inch=0.15,ylab="Bins",xlab="Year",zero2NAs=T,...) {
    #This function plots frequency data as bubble plots
    #You may want to change all zeros to NA's so that those observations are not plotted.
    #   If you don't then set zero2NAs=F
    x <- as.numeric(as.character(dat$year))
    gender <- dat$gender[1]
    dat <- dat[,-c(1:6)]
    if(gender==0) {
        #dat <- dat[,1:(ncol(dat)/2)]
        dat <- dat[,-match("U0.1",names(dat))]
        #print(names(dat))
        dat <- dat[,-match("U0",names(dat))]
        #print(names(dat))
    }
    if(gender==3) {
        dat <- dat[,-match("F0",names(dat))]
        dat <- dat[,-match("M0",names(dat))]
    }
    numLens <- ncol(dat)/2
    y <- as.numeric(substring(names(dat),2))
    y <- y[1:numLens]
    if(zero2NAs) {dat[dat==0] <- NA}
    
    if(gender==0) {
        z <- c(unlist(dat[,1:numLens]),max(dat))
        symbols(c(rep(x,length(y)),0),c(rep(y,each=length(x)),0),circles=z,main="Unsexed+Males+Females",inches=inch,xlab=xlab,ylab=ylab,xlim=range(x),...)
    }
    if(gender==3) {
        z <- c(unlist(dat[,1:numLens]),max(dat))
        symbols(c(rep(x,length(y)),0),c(rep(y,each=length(x)),0),circles=z,main="Female",inches=inch,xlab=xlab,ylab=ylab,xlim=range(x),...)
        z <- c(unlist(dat[,(numLens+1):ncol(dat)]),max(dat))
        symbols(c(rep(x,length(y)),0),c(rep(y,each=length(x)),0),circles=z,main="Male",inches=inch,xlab=xlab,ylab=ylab,xlim=range(x),...)
    }
}
#lfs <- GetLFs.fn("C:\\NOAA2009\\Petrale\\Data\\Survey\\NWFSCsurvey\\FisheryIndices2009_Petrale.xls",headerRow=7)
#par(mfrow=c(1,2))
#plotFreqData.fn(lfs,zero2NAs=F)
###plot only observations where there is a male and/or female
#lfs <- GetLFs.fn("C:\\NOAA2009\\Petrale\\Data\\Survey\\NWFSCsurvey\\FisheryIndices2009_Petrale.xls",headerRow=7,NAs2zero=F)
#par(mfrow=c(1,2))
#plotFreqData.fn(lfs,zero2NAs=F)
###set all zeros to NAs if you didn't want them plotted
#lfs <- GetLFs.fn("C:\\NOAA2009\\Petrale\\Data\\Survey\\NWFSCsurvey\\FisheryIndices2009_Petrale.xls",headerRow=7)
#par(mfrow=c(1,2))
#plotFreqData.fn(lfs)

plotLFs.fn <- function(dat,inch=0.15,ylab="Bins",xlab="Year",zero2NAs=T,...) {
    #a wrapper for the plotFreqData because I originally called it plotLFs.fn and want to keep it compatible with already written analyses
    plotFreqData.fn(dat,inch,ylab,xlab,zero2NAs,...)
}
#lfs <- GetLFs.fn("C:\\NOAA2009\\Petrale\\Data\\Survey\\NWFSCsurvey\\FisheryIndices2009_Petrale.xls",headerRow=7)
#par(mfrow=c(1,2))
#plotLFs.fn(lfs,zero2NAs=F)


plotFreqDataStrata.fn <- function(dat,inch=0.15,ylab="Bins",xlab="",zero2NAs=T,...) {
    #This function plots frequency data as bubble plots when strata instead of years are used
    #You may want to change all zeros to NA's so that those observations are not plotted.
    #   If you don't then set zero2NAs=F
    x <- 1:nrow(dat)
    xlabels <- as.character(dat$year)
    gender <- dat$gender[1]
    dat <- dat[,-c(1:6)]
    if(gender==0) {
        #dat <- dat[,1:(ncol(dat)/2)]
        dat <- dat[,-match("U0.1",names(dat))]
        dat <- dat[,-match("U0",names(dat))]
    }
    if(gender==3) {
        dat <- dat[,-match("F0",names(dat))]
        dat <- dat[,-match("M0",names(dat))]
    }
    numLens <- ncol(dat)/2
    y <- as.numeric(substring(names(dat),2))
    y <- y[1:numLens]
    if(zero2NAs) {dat[dat==0] <- NA}
    
    if(gender==0) {
        print(dat$year)
        z <- c(unlist(dat[,1:numLens]),max(dat))
        symbols(c(rep(x,length(y)),0),c(rep(y,each=length(x)),0),circles=z,main="Unsexed+Males+Females",inches=inch,xlab=xlab,ylab=ylab,xlim=range(x),xaxt="n",...)
        axis(1,at=x,label=xlabels,padj=0,las=2,cex.axis=0.5)
    }
    if(gender==3) {
        z <- c(unlist(dat[,1:numLens]),max(dat))
        symbols(c(rep(x,length(y)),0),c(rep(y,each=length(x)),0),circles=z,main="Female",inches=inch,xlab=xlab,ylab=ylab,xlim=range(x),xaxt="n",...)
        axis(1,at=x,label=xlabels,padj=0,las=2,cex.axis=0.5)
        z <- c(unlist(dat[,(numLens+1):ncol(dat)]),max(dat))
        symbols(c(rep(x,length(y)),0),c(rep(y,each=length(x)),0),circles=z,main="Male",inches=inch,xlab=xlab,ylab=ylab,xlim=range(x),xaxt="n",...)
        axis(1,at=x,label=xlabels,padj=0,las=2,cex.axis=0.5)
    }
}


#Treating the ages the same as the LFs just to plot the expanded age frequency
#PetA <- readInExcelAgeComps.fn("C:\\NOAA2009\\Petrale\\Data\\Survey\\NWFSCsurvey\\FisheryIndices2009_Petrale.xls")
#PetA$Length <- PetA$Age
#tmp <- SS3LF.fn(PetA,lgthBins=1:17,gender=3)
#plotLFs.fn(tmp)
