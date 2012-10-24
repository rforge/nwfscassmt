\name{plotFreqDataStrata.fn}
\alias{plotFreqDataStrata.fn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot age or length frequencies by strata as bubble plots
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotFreqDataStrata.fn(dat, inch = 0.15, ylab = "Bins", xlab = "", zero2NAs = T, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
%%     ~~Describe \code{dat} here~~
}
  \item{inch}{
%%     ~~Describe \code{inch} here~~
}
  \item{ylab}{
%%     ~~Describe \code{ylab} here~~
}
  \item{xlab}{
%%     ~~Describe \code{xlab} here~~
}
  \item{zero2NAs}{
%%     ~~Describe \code{zero2NAs} here~~
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
function (dat, inch = 0.15, ylab = "Bins", xlab = "", zero2NAs = T, 
    ...) 
{
    x <- 1:nrow(dat)
    xlabels <- as.character(dat$year)
    gender <- dat$gender[1]
    dat <- dat[, -c(1:6)]
    if (gender == 0) {
        dat <- dat[, -match("U0.1", names(dat))]
        dat <- dat[, -match("U0", names(dat))]
    }
    if (gender == 3) {
        dat <- dat[, -match("F0", names(dat))]
        dat <- dat[, -match("M0", names(dat))]
    }
    numLens <- ncol(dat)/2
    y <- as.numeric(substring(names(dat), 2))
    y <- y[1:numLens]
    if (zero2NAs) {
        dat[dat == 0] <- NA
    }
    if (gender == 0) {
        print(dat$year)
        z <- c(unlist(dat[, 1:numLens]), max(dat))
        symbols(c(rep(x, length(y)), 0), c(rep(y, each = length(x)), 
            0), circles = z, main = "Unsexed+Males+Females", 
            inches = inch, xlab = xlab, ylab = ylab, xlim = range(x), 
            xaxt = "n", ...)
        axis(1, at = x, label = xlabels, padj = 0, las = 2, cex.axis = 0.5)
    }
    if (gender == 3) {
        z <- c(unlist(dat[, 1:numLens]), max(dat))
        symbols(c(rep(x, length(y)), 0), c(rep(y, each = length(x)), 
            0), circles = z, main = "Female", inches = inch, 
            xlab = xlab, ylab = ylab, xlim = range(x), xaxt = "n", 
            ...)
        axis(1, at = x, label = xlabels, padj = 0, las = 2, cex.axis = 0.5)
        z <- c(unlist(dat[, (numLens + 1):ncol(dat)]), max(dat))
        symbols(c(rep(x, length(y)), 0), c(rep(y, each = length(x)), 
            0), circles = z, main = "Male", inches = inch, xlab = xlab, 
            ylab = ylab, xlim = range(x), xaxt = "n", ...)
        axis(1, at = x, label = xlabels, padj = 0, las = 2, cex.axis = 0.5)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
