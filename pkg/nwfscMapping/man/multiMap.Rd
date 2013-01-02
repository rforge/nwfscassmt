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
    library(PBSmapping)
    data(westCoastLL)
    data(WCstatesInlandPBS)
    #plot the dataset
    multiMap(nrows=1,ncols=2,lonLims=matrix(c(-125,-125,-122,-123),nrow=2),latLims=matrix(c(42,45,45,46),nrow=2),fnc=plotMap,polys=westCoastLL)
    multiMap(nrows=1,ncols=3,
             lonLims=matrix(c(-123.1,-122.9,-125,-121,-125,-121),nrow=3,byrow=T),
             latLims=matrix(c(37,40,40,45,45,48),nrow=3,byrow=T),
             fnc=plotMap,polys=westCoastLL)

    #plot using mapBubbles function
    library(PBSmapping)
    data(westCoastLL)
    data(WCstatesInlandPBS)
    #create a fake dataset
    set.seed(55)
    datt <- data.frame(Latitude=c(38,43,44,45,46,47,runif(100,37,47)),
                       Longitude=c(-123.5,-125.5,-128,-127.5,-126.5,-126.5,runif(100,-128,-124)),Z=c(1,2,3,4,5,6,rep(0,100)))
    #mapBiomass.fn(datt,xlim=c(-130,-120),ylim=c(40,49))
    windows(height=7,width=10)
    multiMap(nrows=1,ncols=3,
             lonLims=matrix(c(-126,-121,-129,-122.8,-129,-122.3),nrow=3,byrow=T),
             latLims=matrix(c(37,40,40,45,45,48),nrow=3,byrow=T),
             fnc=mapBubbles,dat=datt,col="green4",bg="lightblue")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
