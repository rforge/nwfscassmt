\name{SS3AgeAtLen.fn}
\alias{SS3AgeAtLen.fn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SS3AgeAtLen.fn(ages, lgthBins = 1, ageBins = 1, fleet = "EnterFleet", season = 1, partition = 0, ageerr = "EnterAgeErr", raw = T, sexRatioUnsexed = NA, maxSizeUnsexed = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ages}{
%%     ~~Describe \code{ages} here~~
}
  \item{lgthBins}{
%%     ~~Describe \code{lgthBins} here~~
}
  \item{ageBins}{
%%     ~~Describe \code{ageBins} here~~
}
  \item{fleet}{
%%     ~~Describe \code{fleet} here~~
}
  \item{season}{
%%     ~~Describe \code{season} here~~
}
  \item{partition}{
%%     ~~Describe \code{partition} here~~
}
  \item{ageerr}{
%%     ~~Describe \code{ageerr} here~~
}
  \item{raw}{
%%     ~~Describe \code{raw} here~~
}
  \item{sexRatioUnsexed}{
%%     ~~Describe \code{sexRatioUnsexed} here~~
}
  \item{maxSizeUnsexed}{
%%     ~~Describe \code{maxSizeUnsexed} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (ages, lgthBins = 1, ageBins = 1, fleet = "EnterFleet", 
    season = 1, partition = 0, ageerr = "EnterAgeErr", raw = T, 
    sexRatioUnsexed = NA, maxSizeUnsexed = NA) 
{
    years <- sort(unique(ages$Year))
    if (length(lgthBins) == 1) {
        Lengths <- c(0, seq(floor(min(ages$Length)), ceiling(max(ages$Length)), 
            lgthBins), Inf)
    }
    else {
        Lengths <- c(0, lgthBins, Inf)
    }
    if (length(ageBins) == 1) {
        Ages <- c(0, seq(floor(min(ages$Age)), ceiling(max(ages$Age)), 
            ageBins), Inf)
    }
    else {
        Ages <- c(0, ageBins, Inf)
    }
    allLs <- Lengths[findInterval(ages$Length, Lengths, all.inside = T)]
    if (length(sexRatioUnsexed) == 1 & !is.na(sexRatioUnsexed)) {
        ages$allLs <- allLs
        ages$sexRatio <- ages$NumF/(ages$NumF + ages$NumM)
        ages$sexRatio[ages$Length <= maxSizeUnsexed] <- sexRatioUnsexed
        noRatio <- which(is.na(ages$sexRatio))
        if (length(noRatio) > 0) 
            cat("\nThese are sex ratios that were filled in using observations from the same lengths from different strata and years\n")
        for (i in noRatio) {
            inds <- allLs == allLs[i]
            tmpF <- sum(ages$NumF[inds])
            tmpM <- sum(ages$NumM[inds])
            ages$sexRatio[i] <- tmpF/(tmpF + tmpM)
            print(ages[i, c("Length", "allLs", "Age", "NumF", 
                "NumM", "sexRatio")])
        }
        noRatio <- which(is.na(ages$sexRatio))
        if (length(noRatio) > 0) 
            cat("\nThese are sex ratios that were filled in using observations from nearby lengths\n")
        for (i in noRatio) {
            nearLens <- Lengths[c(which(Lengths == allLs[i]) - 
                1, which(Lengths == allLs[i]) + 1)]
            inds <- ages$allLs \%in\% nearLens
            tmpF <- sum(ages$NumF[inds])
            tmpM <- sum(ages$NumM[inds])
            ages$sexRatio[i] <- tmpF/(tmpF + tmpM)
            print(ages[i, c("Length", "allLs", "Age", "NumF", 
                "NumM", "sexRatio")])
        }
        noRatio <- which(is.na(ages$sexRatio))
        if (length(noRatio) > 0) 
            cat("Some sex ratios were left unknown and omitted\n\n")
        if (length(noRatio) == 0) 
            cat("Done filling in sex ratios\n\n")
        tmpFemUnsex <- round(ages$sexRatio * ages$NumUnsexed)
        tmpMaleUnsex <- ages$NumUnsexed - tmpFemUnsex
        ages$NumF <- ages$NumF + tmpFemUnsex
        ages$NumM <- ages$NumM + tmpMaleUnsex
        print(unique(round(ages$sexRatio, 1)))
    }
    if (raw) {
        xx <- split(ages[, c("Year", "Length", "Age", "AgeTallyF", 
            "AgeTallyM")], paste(ages$Year, allLs))
    }
    if (!raw) {
        ages[, c("AgeTallyF", "AgeTallyM")] <- ages[, c("NumF", 
            "NumM")]
        xx <- split(ages[, c("Year", "Length", "Age", "AgeTallyF", 
            "AgeTallyM")], paste(ages$Year, allLs))
    }
    bin.fn <- function(x, Ages) {
        allAs <- Ages[findInterval(x$Age, Ages, all.inside = T)]
        totalF <- tapply(x$AgeTallyF, allAs, sum, na.rm = T)
        totalM <- tapply(x$AgeTallyM, allAs, sum, na.rm = T)
        out <- data.frame(Age = Ages, numF = rep(NA, length(Ages)), 
            numM = rep(NA, length(Ages)))
        row.names(out) <- out$Age
        out[names(totalF), "numF"] <- totalF
        out[names(totalM), "numM"] <- totalM
        out[names(totalF), "propF"] <- 100 * totalF/sum(totalF, 
            na.rm = T)
        out[names(totalM), "propM"] <- 100 * totalM/sum(totalM, 
            na.rm = T)
        out <- out[-nrow(out), ]
        return(out)
    }
    A.bin <- lapply(xx, bin.fn, Ages = Ages)
    Nobs.fn <- function(x) {
        nF <- sum(x$AgeTallyF, na.rm = T)
        nM <- sum(x$AgeTallyM, na.rm = T)
        out <- c(nF, nM)
        names(out) <- c("nF", "nM")
        return(out)
    }
    nobs <- lapply(xx, Nobs.fn)
    ages <- as.character(A.bin[[1]]$Age)
    nsF <- unlist(lapply(nobs, function(x) {
        x["nF"]
    }))
    nsM <- unlist(lapply(nobs, function(x) {
        x["nM"]
    }))
    AsF <- unlist(lapply(A.bin, function(x) {
        x$propF
    }))
    AsM <- unlist(lapply(A.bin, function(x) {
        x$propM
    }))
    AsF[is.na(AsF)] <- 0
    AsM[is.na(AsM)] <- 0
    AsF <- matrix(AsF, nrow = length(A.bin), byrow = T, dimnames = list(NULL, 
        paste(rep("F", length(ages)), ages, sep = "")))
    AsF[, "F1"] <- AsF[, "F0"] + AsF[, "F1"]
    numFzero <- sum(AsF[, "F0"])
    AsF <- AsF[, -match("F0", dimnames(AsF)[[2]])]
    AsM <- matrix(AsM, nrow = length(A.bin), byrow = T, dimnames = list(NULL, 
        paste(rep("M", length(ages)), ages, sep = "")))
    AsM[, "M1"] <- AsM[, "M0"] + AsM[, "M1"]
    numMzero <- sum(AsM[, "M0"])
    AsM <- AsM[, -match("M0", dimnames(AsM)[[2]])]
    outF <- data.frame(year = as.numeric(substring(names(A.bin), 
        1, 4)), Season = season, Fleet = fleet, gender = 1, partition = partition, 
        ageErr = ageerr, LbinLo = as.numeric(substring(names(A.bin), 
            6)), LbinHi = as.numeric(substring(names(A.bin), 
            6)), nSamps = nsF, AsF, AsF)
    outM <- data.frame(year = as.numeric(substring(names(A.bin), 
        1, 4)), Season = season, Fleet = fleet, gender = 2, partition = partition, 
        ageErr = ageerr, LbinLo = as.numeric(substring(names(A.bin), 
            6)), LbinHi = as.numeric(substring(names(A.bin), 
            6)), nSamps = nsM, AsM, AsM)
    indZero <- apply(outF[, -c(1:9)], 1, sum) == 0
    outF <- outF[!indZero, ]
    indZero <- apply(outM[, -c(1:9)], 1, sum) == 0
    outM <- outM[!indZero, ]
    rownames(outF) <- paste("F", 1:nrow(outF), sep = "")
    rownames(outM) <- paste("M", 1:nrow(outM), sep = "")
    cat("There are", numFzero, "females in the age 0 to age", 
        ages[2], "that were added into the first age bin\n")
    cat("There are", numMzero, "males in the age 0 to age", ages[2], 
        "that were added into the first age bin\n")
    cat("The number of fish in each sample were input into the nSamps column\nUse Beth's Excel file for the number of tows")
    return(list(female = outF, male = outM))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
