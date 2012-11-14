\name{SimulatorFn}
\alias{SimulatorFn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SimulatorFn(Nreaders, M, SelexForm, ErrorParams, BiasParams, SelexParams, ReadsMat, RecCv = 0.6, RecAr1 = 0.8, Amax = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Nreaders}{
%%     ~~Describe \code{Nreaders} here~~
}
  \item{M}{
%%     ~~Describe \code{M} here~~
}
  \item{SelexForm}{
%%     ~~Describe \code{SelexForm} here~~
}
  \item{ErrorParams}{
%%     ~~Describe \code{ErrorParams} here~~
}
  \item{BiasParams}{
%%     ~~Describe \code{BiasParams} here~~
}
  \item{SelexParams}{
%%     ~~Describe \code{SelexParams} here~~
}
  \item{ReadsMat}{
%%     ~~Describe \code{ReadsMat} here~~
}
  \item{RecCv}{
%%     ~~Describe \code{RecCv} here~~
}
  \item{RecAr1}{
%%     ~~Describe \code{RecAr1} here~~
}
  \item{Amax}{
%%     ~~Describe \code{Amax} here~~
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
function (Nreaders, M, SelexForm, ErrorParams, BiasParams, SelexParams, 
    ReadsMat, RecCv = 0.6, RecAr1 = 0.8, Amax = 100) 
{
    RecDev = rnorm(Amax, mean = 0, sd = RecCv)
    for (i in 2:length(RecDev)) {
        RecDev[i] = RecDev[i] * sqrt(1 - RecAr1) + RecDev[i - 
            1] * sqrt(RecAr1)
    }
    AgeStruct = 1 * exp(-M * 1:Amax) * exp(RecDev - RecCv^2/2)
    if (SelexForm == "Logistic") {
        SelexAtAge = 1/(1 + exp((SelexParams[1] - 1:Amax) * SelexParams[2]))
    }
    else {
        stop("Selex not implemented")
    }
    Ages = sample(x = 1:Amax, size = sum(ReadsMat[, 1]), prob = AgeStruct * 
        SelexAtAge, replace = TRUE)
    AgeReads = vector()
    IndexI = 0
    for (RowI in 1:nrow(ReadsMat)) {
        for (OtolithI in 1:ReadsMat[RowI, 1]) {
            IndexI = IndexI + 1
            Row = vector()
            for (ReaderI in 1:Nreaders) {
                if (ReadsMat[RowI, ReaderI + 1] == 0) {
                  Row = c(Row, NA)
                }
                else {
                  Row = c(Row, round(Ages[IndexI] * BiasParams[ReaderI] + 
                    rnorm(1, mean = 0, sd = Ages[IndexI] * ErrorParams[ReaderI])))
                }
            }
            AgeReads = rbind(AgeReads, Row)
        }
    }
    return(AgeReads)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
