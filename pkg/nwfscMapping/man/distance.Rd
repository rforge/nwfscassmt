\name{distance}
\alias{distance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
calculate distance (km) between two points
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
distance(long1, lat1, long2, lat2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{long1}{
%%     ~~Describe \code{long1} here~~
}
  \item{lat1}{
%%     ~~Describe \code{lat1} here~~
}
  \item{long2}{
%%     ~~Describe \code{long2} here~~
}
  \item{lat2}{
%%     ~~Describe \code{lat2} here~~
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
function (long1, lat1, long2, lat2) 
{
    sqrt(((long2 - long1) * long.to.km((lat1 + lat2)/2))^2 + 
        ((lat2 - lat1) * lat.to.km((lat1 + lat2)/2))^2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
