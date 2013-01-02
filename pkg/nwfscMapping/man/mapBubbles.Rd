\name{mapBubbles}
\alias{mapBubbles}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot a map with bubbles indicating size of attribute
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mapBubbles(dat, depths = NULL, plotMap = T, maxValue = max(dat$Z), maxSize = 0.3, fgCircle = "black", bgCircle = rgb(0, 0, 0, 0.3), type = "surface", zeroSize = 0.04, fgZero = gray(0.5), bgZero = gray(0.5), legendPos = "bottomright", legendBreaks = NULL, legendType = "vert", legendTitle = NULL, legendCex = 0.8, legendTitleCex = legendCex, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
%%     ~~Describe \code{dat} here~~
}
  \item{depths}{
%%     ~~Describe \code{depths} here~~
}
  \item{plotMap}{
%%     ~~Describe \code{plotMap} here~~
}
  \item{maxValue}{
%%     ~~Describe \code{maxValue} here~~
}
  \item{maxSize}{
%%     ~~Describe \code{maxSize} here~~
}
  \item{fgCircle}{
%%     ~~Describe \code{fgCircle} here~~
}
  \item{bgCircle}{
%%     ~~Describe \code{bgCircle} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{zeroSize}{
%%     ~~Describe \code{zeroSize} here~~
}
  \item{fgZero}{
%%     ~~Describe \code{fgZero} here~~
}
  \item{bgZero}{
%%     ~~Describe \code{bgZero} here~~
}
  \item{legendPos}{
%%     ~~Describe \code{legendPos} here~~
}
  \item{legendBreaks}{
%%     ~~Describe \code{legendBreaks} here~~
}
  \item{legendType}{
%%     ~~Describe \code{legendType} here~~
}
  \item{legendTitle}{
%%     ~~Describe \code{legendTitle} here~~
}
  \item{legendCex}{
%%     ~~Describe \code{legendCex} here~~
}
  \item{legendTitleCex}{
%%     ~~Describe \code{legendTitleCex} here~~
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
function (dat, depths = NULL, plotMap = T, maxValue = max(dat$Z), 
    maxSize = 0.3, fgCircle = "black", bgCircle = rgb(0, 0, 0, 
        0.3), type = "surface", zeroSize = 0.04, fgZero = gray(0.5), 
    bgZero = gray(0.5), legendPos = "bottomright", legendBreaks = NULL, 
    legendType = "vert", legendTitle = NULL, legendCex = 0.8, 
    legendTitleCex = legendCex, ...) 
{
    if (plotMap) {
        plotMap(westCoastLL, tck = c(-0.02), ...)
        if (exists("westCoastLL")) {
            addLines(westCoastLL)
        }
        else {
            print("Run data(westCoastLL) before calling this function to plot West Coast shoreline")
        }
        if (exists("WCstatesInlandPBS")) {
            addLines(WCstatesInlandPBS)
        }
        else {
            print("Run data(WCstatesInlandPBS) before calling this function to plot state boundaries")
        }
        if (!is.null(depths)) {
            for (i in 1:length(depths)) {
                addLines(depths[[i]], col = gray(0.5))
            }
        }
    }
    ind0 <- dat$Z <= 0
    ind <- dat$Z > 0
    symbols(dat[ind0, "Longitude"], dat[ind0, "Latitude"], squares = rep(1, 
        sum(ind0)), add = T, fg = fgZero, bg = fgZero, inches = zeroSize)
    tmp <- dat[ind, ]
    event <- as.EventData(data.frame(EID = 1:nrow(tmp), X = tmp$Longitude, 
        Y = tmp$Latitude, Z = tmp$Z))
    .addBubbles.ach(event, type = type, z.max = maxValue, max.size = maxSize, 
        symbol.fg = fgCircle, symbol.bg = bgCircle, symbol.zero = F, 
        legend.pos = legendPos, legend.breaks = legendBreaks, 
        legend.type = legendType, legend.title = legendTitle, 
        legend.cex = legendCex, legend.title.cex = legendTitleCex)
    invisible()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
