\name{plotSexRatio.fn}
\alias{plotSexRatio.fn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots the sex ratio by length bins
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotSexRatio.fn(len, fn = median, circleSize = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{len}{
%%     ~~Describe \code{len} here~~
}
  \item{fn}{
%%     ~~Describe \code{fn} here~~
}
  \item{circleSize}{
%%     ~~Describe \code{circleSize} here~~
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
function (len, fn = median, circleSize = 0.1) 
{
    ratioF <- len$NumF/(len$NumF + len$NumM)
    yF <- lapply(split(ratioF, floor(len$Length)), fn)
    x <- names(split(ratioF, floor(len$Length)))
    nobs <- unlist(lapply(split(ratioF, floor(len$Length)), length))
    plot(x, yF, type = "l", col = "red")
    symbols(x, yF, circles = nobs, inches = circleSize, fg = "red", 
        bg = rgb(1, 0, 0, alpha = 0.5), add = T)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
