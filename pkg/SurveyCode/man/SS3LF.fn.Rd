\name{SS3LF.fn}
\alias{SS3LF.fn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SS3LF.fn(len, lgthBins = 1, gender = 3, nSamps = "EnterNsamps", fleet = "EnterFleet", season = 1, partition = 0, NAs2zero = T, sexRatioUnsexed = NA, maxSizeUnsexed = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{len}{
%%     ~~Describe \code{len} here~~
}
  \item{lgthBins}{
%%     ~~Describe \code{lgthBins} here~~
}
  \item{gender}{
%%     ~~Describe \code{gender} here~~
}
  \item{nSamps}{
%%     ~~Describe \code{nSamps} here~~
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
  \item{NAs2zero}{
%%     ~~Describe \code{NAs2zero} here~~
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
function (len, lgthBins = 1, gender = 3, nSamps = "EnterNsamps", 
    fleet = "EnterFleet", season = 1, partition = 0, NAs2zero = T, 
    sexRatioUnsexed = NA, maxSizeUnsexed = NA) 
{
    if (length(lgthBins) == 1) {
        Lengths <- c(0, seq(floor(min(len$Length)), ceiling(max(len$Length)), 
            lgthBins), Inf)
    }
    else {
        Lengths <- c(0, lgthBins, Inf)
    }
    len$allLs <- Lengths[findInterval(len$Length, Lengths, all.inside = T)]
    if (length(sexRatioUnsexed) == 1 & !is.na(sexRatioUnsexed)) {
        len$sexRatio <- len$NumF/(len$NumF + len$NumM)
        len$sexRatio[len$Length <= maxSizeUnsexed] <- sexRatioUnsexed
        noRatio <- which(is.na(len$sexRatio))
        if (length(noRatio) > 0) 
            cat("\nThese are sex ratios that were filled in using observations from the same lengths from different strata and years\n")
        for (i in noRatio) {
            inds <- len$allLs == len$allLs[i]
            tmpF <- sum(len$NumF[inds])
            tmpM <- sum(len$NumM[inds])
            len$sexRatio[i] <- tmpF/(tmpF + tmpM)
            print(len[i, c("Length", "allLs", "NumF", "NumM", 
                "sexRatio")])
        }
        noRatio <- which(is.na(len$sexRatio))
        if (length(noRatio) > 0) 
            cat("\nThese are sex ratios that were filled in using observations from nearby lengths\n")
        for (i in noRatio) {
            nearLens <- Lengths[c(which(Lengths == len$allLs[i]) - 
                1, which(Lengths == len$allLs[i]) + 1)]
            inds <- len$allLs \%in\% nearLens
            tmpF <- sum(len$NumF[inds])
            tmpM <- sum(len$NumM[inds])
            len$sexRatio[i] <- tmpF/(tmpF + tmpM)
            print(len[i, c("Length", "allLs", "NumF", "NumM", 
                "sexRatio")])
        }
        noRatio <- which(is.na(len$sexRatio))
        if (length(noRatio) > 0) 
            cat("Some sex ratios were left unknown and omitted\n\n")
        if (length(noRatio) == 0) 
            cat("Done filling in sex ratios\n\n")
        tmpFemUnsex <- round(len$sexRatio * len$NumUnsexed)
        tmpMaleUnsex <- len$NumUnsexed - tmpFemUnsex
        len$NumF <- len$NumF + tmpFemUnsex
        len$NumM <- len$NumM + tmpMaleUnsex
    }
    xx <- split(len[, c("Length", "NumF", "NumM", "NumUnsexed")], 
        len$Year)
    years <- names(xx)
    year.fn <- function(x, Lengths) {
        allLs <- Lengths[findInterval(x$Length, Lengths, all.inside = T)]
        totalU <- tapply(x$NumUnsexed, allLs, sum, na.rm = T)
        totalF <- tapply(x$NumF, allLs, sum, na.rm = T)
        totalM <- tapply(x$NumM, allLs, sum, na.rm = T)
        out <- data.frame(Length = Lengths, numF = rep(NA, length(Lengths)), 
            numM = rep(NA, length(Lengths)), numU = rep(NA, length(Lengths)))
        row.names(out) <- out$Length
        out[names(totalF), "numF"] <- totalF
        out[names(totalM), "numM"] <- totalM
        out[names(totalU), "numU"] <- totalU
        out[names(totalF), "propF"] <- 100 * totalF/(sum(totalF, 
            na.rm = T) + sum(totalM, na.rm = T))
        out[names(totalM), "propM"] <- 100 * totalM/(sum(totalF, 
            na.rm = T) + sum(totalM, na.rm = T))
        out[names(totalU), "propU"] <- 100 * (totalF + totalM)/sum(totalF + 
            totalM, na.rm = T)
        out <- out[-nrow(out), ]
        return(out)
    }
    L.year <- lapply(xx, year.fn, Lengths = Lengths)
    lgths <- as.character(L.year[[1]]$Length)
    if (gender == 0) {
        Ls <- unlist(lapply(L.year, function(x) {
            x$propU
        }))
        if (NAs2zero) {
            Ls[is.na(Ls)] <- 0
        }
        Ls <- matrix(Ls, nrow = length(L.year), byrow = T, dimnames = list(NULL, 
            paste(rep("U", length(lgths)), lgths, sep = "")))
        out <- data.frame(year = as.numeric(names(L.year)), Season = season, 
            Fleet = fleet, gender = gender, partition = partition, 
            nSamps = nSamps, Ls, Ls)
    }
    if (gender == 3) {
        Ls <- unlist(lapply(L.year, function(x) {
            c(x$propF, x$propM)
        }))
        if (NAs2zero) {
            Ls[is.na(Ls)] <- 0
        }
        Ls <- matrix(Ls, nrow = length(L.year), byrow = T, dimnames = list(NULL, 
            paste(c(rep("F", length(lgths)), rep("M", length(lgths))), 
                lgths, sep = "")))
        out <- data.frame(year = as.numeric(names(L.year)), Season = season, 
            Fleet = fleet, gender = gender, partition = partition, 
            nSamps = nSamps, Ls)
    }
    cat("\nNOTE: You may need to add the column called F0 and/or M0 to your first length bin\n\tand delete that column.\n\tThese are the proportion of lengths smaller than the first length bin\n\n")
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
