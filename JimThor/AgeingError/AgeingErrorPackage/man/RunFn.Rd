\name{RunFn}
\alias{RunFn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
RunFn(Data, SigOpt, BiasOpt, NDataSets, MinAge, MaxAge, RefAge, MinusAge, PlusAge, MaxSd, MaxExpectedAge, SaveFile, EffSampleSize = 0, Intern = TRUE, AdmbFile = NULL, JustWrite = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{
%%     ~~Describe \code{Data} here~~
}
  \item{SigOpt}{
%%     ~~Describe \code{SigOpt} here~~
}
  \item{BiasOpt}{
%%     ~~Describe \code{BiasOpt} here~~
}
  \item{NDataSets}{
%%     ~~Describe \code{NDataSets} here~~
}
  \item{MinAge}{
%%     ~~Describe \code{MinAge} here~~
}
  \item{MaxAge}{
%%     ~~Describe \code{MaxAge} here~~
}
  \item{RefAge}{
%%     ~~Describe \code{RefAge} here~~
}
  \item{MinusAge}{
%%     ~~Describe \code{MinusAge} here~~
}
  \item{PlusAge}{
%%     ~~Describe \code{PlusAge} here~~
}
  \item{MaxSd}{
%%     ~~Describe \code{MaxSd} here~~
}
  \item{MaxExpectedAge}{
%%     ~~Describe \code{MaxExpectedAge} here~~
}
  \item{SaveFile}{
%%     ~~Describe \code{SaveFile} here~~
}
  \item{EffSampleSize}{
%%     ~~Describe \code{EffSampleSize} here~~
}
  \item{Intern}{
%%     ~~Describe \code{Intern} here~~
}
  \item{AdmbFile}{
%%     ~~Describe \code{AdmbFile} here~~
}
  \item{JustWrite}{
%%     ~~Describe \code{JustWrite} here~~
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
function (Data, SigOpt, BiasOpt, NDataSets, MinAge, MaxAge, RefAge, 
    MinusAge, PlusAge, MaxSd, MaxExpectedAge, SaveFile, EffSampleSize = 0, 
    Intern = TRUE, AdmbFile = NULL, JustWrite = FALSE) 
{
    if (!is.null(AdmbFile)) 
        file.copy(from = paste(AdmbFile, "agemat.exe", sep = ""), 
            to = paste(SaveFile, "agemat.exe", sep = ""), overwrite = TRUE)
    write(c("# Maximum number of readers", ncol(Data) - 1), file = paste(SaveFile, 
        "agemat.dat", sep = ""))
    write(c("# Number of data sets", NDataSets), file = paste(SaveFile, 
        "agemat.dat", sep = ""), append = TRUE)
    write(c("# Number of points per data set", nrow(Data)), file = paste(SaveFile, 
        "agemat.dat", sep = ""), append = TRUE)
    write(c("# Readers per data set", ncol(Data) - 1), file = paste(SaveFile, 
        "agemat.dat", sep = ""), append = TRUE)
    write("# Which readers per data set", file = paste(SaveFile, 
        "agemat.dat", sep = ""), append = TRUE)
    write.table(rMx(1:(ncol(Data) - 1)), file = paste(SaveFile, 
        "agemat.dat", sep = ""), append = TRUE, row.names = FALSE, 
        col.names = FALSE)
    write(c("# Minimum age", MinAge), file = paste(SaveFile, 
        "agemat.dat", sep = ""), append = TRUE)
    write(c("# Maximum age", MaxAge), file = paste(SaveFile, 
        "agemat.dat", sep = ""), append = TRUE)
    write(c("# Reference age", RefAge), file = paste(SaveFile, 
        "agemat.dat", sep = ""), append = TRUE)
    write(c("# Minus groups", MinusAge), file = paste(SaveFile, 
        "agemat.dat", sep = ""), append = TRUE)
    write(c("# Plus groups", PlusAge), file = paste(SaveFile, 
        "agemat.dat", sep = ""), append = TRUE)
    write("# Option for bias", file = paste(SaveFile, "agemat.dat", 
        sep = ""), append = TRUE)
    write.table(rMx(BiasOpt), file = paste(SaveFile, "agemat.dat", 
        sep = ""), append = TRUE, row.names = FALSE, col.names = FALSE)
    write("# Option for standard deviation", file = paste(SaveFile, 
        "agemat.dat", sep = ""), append = TRUE)
    write.table(rMx(SigOpt), file = paste(SaveFile, "agemat.dat", 
        sep = ""), append = TRUE, row.names = FALSE, col.names = FALSE)
    write(c("# Option for effective sample size", EffSampleSize), 
        file = paste(SaveFile, "agemat.dat", sep = ""), append = TRUE)
    write(c("# Use Par File (1=Yes)", 0), file = paste(SaveFile, 
        "agemat.dat", sep = ""), append = TRUE)
    write("\n# Min, Max, Init, Phase for Bias", file = paste(SaveFile, 
        "agemat.dat", sep = ""), append = TRUE)
    for (BiasI in 1:length(BiasOpt)) {
        if (BiasOpt[BiasI] <= 0) {
        }
        if (BiasOpt[BiasI] == 1) 
            write.table(rMx(c(0.001, 3, 1, 2)), file = paste(SaveFile, 
                "agemat.dat", sep = ""), append = TRUE, col.names = FALSE, 
                row.names = FALSE)
        if (BiasOpt[BiasI] == 2) {
            write.table(rMx(c(0.001, 10, 1, 2)), file = paste(SaveFile, 
                "agemat.dat", sep = ""), append = TRUE, col.names = FALSE, 
                row.names = FALSE)
            write.table(rMx(c(-10, 1, 0.1, 2)), file = paste(SaveFile, 
                "agemat.dat", sep = ""), append = TRUE, col.names = FALSE, 
                row.names = FALSE)
            write.table(rMx(c(0.001, MaxAge * 2, MaxAge, 2)), 
                file = paste(SaveFile, "agemat.dat", sep = ""), 
                append = TRUE, col.names = FALSE, row.names = FALSE)
        }
    }
    write("\n# Min, Max, Init, Phase for Sigma", file = paste(SaveFile, 
        "agemat.dat", sep = ""), append = TRUE)
    for (SigI in 1:length(SigOpt)) {
        if (SigOpt[SigI] <= 0) {
        }
        if (SigOpt[SigI] == 1) 
            write.table(rMx(c(0.001, 3, 0.1, 2)), file = paste(SaveFile, 
                "agemat.dat", sep = ""), append = TRUE, col.names = FALSE, 
                row.names = FALSE)
        if (SigOpt[SigI] == 2) {
            write.table(rMx(c(0.001, 100, 1, 2)), file = paste(SaveFile, 
                "agemat.dat", sep = ""), append = TRUE, col.names = FALSE, 
                row.names = FALSE)
            write.table(rMx(c(-10, 1, 0.1, 2)), file = paste(SaveFile, 
                "agemat.dat", sep = ""), append = TRUE, col.names = FALSE, 
                row.names = FALSE)
            write.table(rMx(c(0.001, 100, 10, 2)), file = paste(SaveFile, 
                "agemat.dat", sep = ""), append = TRUE, col.names = FALSE, 
                row.names = FALSE)
        }
        if (SigOpt[SigI] == 3) {
            write.table(rMx(c(0.001, 3, 0.1, 2)), file = paste(SaveFile, 
                "agemat.dat", sep = ""), append = TRUE, col.names = FALSE, 
                row.names = FALSE)
            write.table(rMx(c(-10, 1, 0.1, 2)), file = paste(SaveFile, 
                "agemat.dat", sep = ""), append = TRUE, col.names = FALSE, 
                row.names = FALSE)
            write.table(rMx(c(0.001, 3, 0.1, 2)), file = paste(SaveFile, 
                "agemat.dat", sep = ""), append = TRUE, col.names = FALSE, 
                row.names = FALSE)
        }
    }
    write("\n# Min, Max, Phase for Probs", file = paste(SaveFile, 
        "agemat.dat", sep = ""), append = TRUE)
    write.table(rMx(c(-20, 20, 2)), file = paste(SaveFile, "agemat.dat", 
        sep = ""), append = TRUE, col.names = FALSE, row.names = FALSE)
    write("\n# Min, Max, Init, Phase for slopes", file = paste(SaveFile, 
        "agemat.dat", sep = ""), append = TRUE)
    for (DataSetI in 1:NDataSets) {
        if (MinAge > 1) 
            write.table(rMx(c(-10, 0, 0, 1)), file = paste(SaveFile, 
                "agemat.dat", sep = ""), append = TRUE, col.names = FALSE, 
                row.names = FALSE)
        write.table(rMx(c(-10, 0, 0, 1)), file = paste(SaveFile, 
            "agemat.dat", sep = ""), append = TRUE, col.names = FALSE, 
            row.names = FALSE)
    }
    write("\n# Data set", file = paste(SaveFile, "agemat.dat", 
        sep = ""), append = TRUE)
    write.table(Data, file = paste(SaveFile, "agemat.dat", sep = ""), 
        append = TRUE, row.names = FALSE, col.names = FALSE)
    write(c("# Test number", 123456), file = paste(SaveFile, 
        "agemat.dat", sep = ""), append = TRUE)
    if (JustWrite == FALSE) {
        setwd(SaveFile)
        Output = shell("agemat.exe -est", intern = Intern)
        Admb = scan(paste(SaveFile, "agemat.par", sep = ""), 
            comment.char = "#", quiet = TRUE)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
