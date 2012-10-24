\name{plotBioStrata.fn}
\alias{plotBioStrata.fn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot biomass estimates by strata with CI's
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotBioStrata.fn(bio, CI = 0.95, scalar = 1e+06, gap = 0.03, ylab = "Biomass ('000 mt)", xlab = "Stratum", ylim = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bio}{
%%     ~~Describe \code{bio} here~~
}
  \item{CI}{
%%     ~~Describe \code{CI} here~~
}
  \item{scalar}{
%%     ~~Describe \code{scalar} here~~
}
  \item{gap}{
%%     ~~Describe \code{gap} here~~
}
  \item{ylab}{
%%     ~~Describe \code{ylab} here~~
}
  \item{xlab}{
%%     ~~Describe \code{xlab} here~~
}
  \item{ylim}{
%%     ~~Describe \code{ylim} here~~
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
function (bio, CI = 0.95, scalar = 1e+06, gap = 0.03, ylab = "Biomass ('000 mt)", 
    xlab = "Stratum", ylim = NULL, ...) 
{
    y <- as.numeric(as.character(bio$Value/scalar))
    x <- 1:nrow(bio)
    se <- as.numeric(as.character(bio$seLogB))
    logB <- log(bio$Value)
    ci <- exp(rbind(c(logB + qnorm(1 - (1 - CI)/2) * se), c(logB - 
        qnorm(1 - (1 - CI)/2) * se)))/scalar
    if (is.null(ylim)) {
        ylim <- c(0, 1.05 * max(ci))
    }
    gap <- gap * max(y)
    plot(x, y, ylab = ylab, xlab = xlab, ylim = ylim, xaxt = "n", 
        ...)
    axis(1, at = x, label = as.character(bio$Stratum), padj = 0, 
        las = 2, cex.axis = 0.5)
    segments(x, y + gap, x, ci[1, ])
    segments(x, y - gap, x, ci[2, ])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
