\name{multiMap}
\alias{multiMap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a multi-panel figure with maps having consistent scales
}
\description{
Create a multi-panel figure with maps having consistent scales.
}
\usage{
multiMap(nrows = 1, ncols = 1, lonLims, latLims, adjLonRange = TRUE, fnc, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nrows}{
%%     ~~Describe \code{nrows} here~~
}
  \item{ncols}{
%%     ~~Describe \code{ncols} here~~
}
  \item{lonLims}{
%%     ~~Describe \code{lonLims} here~~
}
  \item{latLims}{
%%     ~~Describe \code{latLims} here~~
}
  \item{adjLonRange}{
%%     ~~Describe \code{adjLonRange} here~~
}
  \item{fnc}{
%%     ~~Describe \code{fnc} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (nrows = 1, ncols = 1, lonLims, latLims, adjLonRange = TRUE, 
    fnc, ...) 
{
    numMaps <- nrows * ncols
    if (ncol(latLims) != 2) 
        stop("latLims has to have two columns for min latitude and max latitude")
    if (ncol(lonLims) != 2) 
        stop("lonLims has to have two columns for min longitude and max longitude")
    if (nrow(latLims) != numMaps) 
        stop("Need 1 row in latLims for each map panel")
    if (nrow(lonLims) != numMaps) 
        stop("Need 1 row in lonLims for each map panel")
    midLat <- apply(latLims, 1, mean)
    midLon <- apply(lonLims, 1, mean)
    dis <- rep(NA, numMaps)
    if (adjLonRange) {
        for (i in 1:numMaps) {
            dis[i] <- distance(lonLims[i, 1], midLat[i], lonLims[i, 
                2], midLat[i])
        }
        for (i in 1:numMaps) {
            diffDis <- max(dis) - dis[i]
            lonLims[i, 1] <- uniroot(function(x, lon2, lat) {
                distance(x, lat, lon2, lat) - diffDis/2
            }, c(-180, lonLims[i, 1]), lon2 = lonLims[i, 1], 
                lat = midLat[i])$root
            lonLims[i, 2] <- uniroot(function(x, lon2, lat) {
                distance(x, lat, lon2, lat) - diffDis/2
            }, c(0, lonLims[i, 2]), lon2 = lonLims[i, 2], lat = midLat[i])$root
        }
    }
    if (!adjLonRange) {
        stop("adjust latitude not implemented yet\n")
    }
    par(mfrow = c(nrows, ncols))
    x <- apply(latLims, 1, diff)
    x <- which(x == max(x))
    fnc(xlim = lonLims[x[1], ], ylim = latLims[x[1], ], ...)
    plt1 <- par()$plt
    par(mfrow = c(nrows, ncols))
    for (i in 1:numMaps) {
        fnc(xlim = lonLims[i, ], ylim = latLims[i, ], plt = NULL, 
            ...)
        par(plt = plt1)
    }
    invisible()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
