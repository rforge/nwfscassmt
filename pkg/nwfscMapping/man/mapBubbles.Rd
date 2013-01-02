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
mapBubbles(dat, depths = NULL, plotMap = T, maxValue = max(dat$Z), maxSize = 0.3, circle.fg = "black", circle.bg = rgb(0, 0, 0, 0.3), type = "surface", zeroSize = 0.04, fgZero = gray(0.5), bgZero = gray(0.5), legendPos = "bottomright", legendBreaks = NULL, legendType = "vert", legendTitle = NULL, legendCex = 0.8, legendTitleCex = legendCex, ...)
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
  \item{circle.fg}{
%%     ~~Describe \code{fgCircle} here~~
}
  \item{circle.bg}{
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
    library(PBSmapping)
    data(westCoastLL)
    data(WCstatesInlandPBS)
    #create a fake dataset
    datt <- data.frame(Latitude=c(42,43,44,45,46,47,runif(10,42,47)),Longitude=c(-126,-125.5,-128,-127.5,-126.5,-126.5,runif(10,-128,-125)),Z=c(1,2,3,4,5,6,rep(0,10)))
    #plot the dataset
    mapBubbles(datt,xlim=c(-130,-120),ylim=c(40,49),col="green4",bg="lightblue")
    #place legend elsewhere
    mapBubbles(datt,xlim=c(-130,-120),ylim=c(40,49),legendPos="bottomleft")
    #nested legend with title
    mapBubbles(datt,xlim=c(-130,-120),ylim=c(40,49),legendBreaks=c(2,6),legendType="nested",legendTitle="Abundance")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
