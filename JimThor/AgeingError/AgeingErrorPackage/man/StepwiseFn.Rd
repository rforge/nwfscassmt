\name{StepwiseFn}
\alias{StepwiseFn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
StepwiseFn(SearchMat, Data, NDataSets, MinAge, MaxAge, RefAge, MaxSd, MaxExpectedAge, SaveFile, EffSampleSize = 0, Intern = TRUE, InformationCriterion = "AIC")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SearchMat}{
%%     ~~Describe \code{SearchMat} here~~
}
  \item{Data}{
%%     ~~Describe \code{Data} here~~
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
  \item{InformationCriterion}{
%%     ~~Describe \code{InformationCriterion} here~~
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
function (SearchMat, Data, NDataSets, MinAge, MaxAge, RefAge, 
    MaxSd, MaxExpectedAge, SaveFile, EffSampleSize = 0, Intern = TRUE, 
    InformationCriterion = "AIC") 
{
    Nages = MaxAge + 1
    Nreaders = ncol(Data) - 1
    ParamVecOpt = SearchMat[, 1]
    Stop = FALSE
    IcRecord = NULL
    StateRecord = NULL
    OuterIndex = 0
    while (Stop == FALSE) {
        OuterIndex = OuterIndex + 1
        Index = 0
        IcVec = NULL
        ParamMat = NULL
        Rep = NULL
        ParamVecOptPreviouslyEstimates = FALSE
        for (VarI in 1:nrow(SearchMat)) {
            for (ValueI in 1:length(na.omit(SearchMat[VarI, ]))) {
                ParamVecCurrent = ParamVecOpt
                ParamVecCurrent[VarI] = na.omit(SearchMat[VarI, 
                  ])[ValueI]
                if (all(ParamVecCurrent == ParamVecOpt) & ParamVecOptPreviouslyEstimates == 
                  FALSE | !all(ParamVecCurrent == ParamVecOpt)) {
                  if (all(ParamVecCurrent == ParamVecOpt)) 
                    ParamVecOptPreviouslyEstimates = TRUE
                  RunFile = paste(SaveFile, "Run\\", sep = "")
                  dir.create(RunFile)
                  file.copy(from = paste(SaveFile, "agemat.exe", 
                    sep = ""), to = paste(RunFile, "agemat.exe", 
                    sep = ""))
                  Index = Index + 1
                  print(paste("Loop=", OuterIndex, " Run=", Index, 
                    " StartTime=", date(), sep = ""))
                  SigOpt = ParamVecCurrent[1:Nreaders]
                  BiasOpt = ParamVecCurrent[Nreaders + 1:Nreaders]
                  MinusAge = ParamVecCurrent[2 * Nreaders + 1]
                  PlusAge = ParamVecCurrent[2 * Nreaders + 2]
                  RunFn(SigOpt = SigOpt, BiasOpt = BiasOpt, Data = Data, 
                    NDataSets = NDataSets, MinAge = MinAge, MaxAge = MaxAge, 
                    RefAge = RefAge, MinusAge = MinusAge, PlusAge = PlusAge, 
                    MaxSd = MaxSd, MaxExpectedAge = MaxExpectedAge, 
                    SaveFile = RunFile)
                  Df = as.numeric(scan(paste(RunFile, "agemat.par", 
                    sep = ""), comment.char = "\%", what = "character", 
                    quiet = TRUE)[6])
                  Nll = as.numeric(scan(paste(RunFile, "agemat.par", 
                    sep = ""), comment.char = "\%", what = "character", 
                    quiet = TRUE)[11])
                  n = sum(ifelse(Data[, -1] == -999, 0, 1))
                  Aic = 2 * Nll + 2 * Df
                  Aicc = Aic + 2 * Df * (Df + 1)/(n - Df - 1)
                  Bic = 2 * Nll + Df * log(n)
                  if (InformationCriterion == "AIC") 
                    IcVec = c(IcVec, Aic)
                  if (InformationCriterion == "AICc") 
                    IcVec = c(IcVec, Aicc)
                  if (InformationCriterion == "BIC") 
                    IcVec = c(IcVec, Bic)
                  ParamMat = rbind(ParamMat, ParamVecCurrent)
                  write.table(cbind(IcVec, ParamMat), paste(SaveFile, 
                    "Stepwise - Model loop ", OuterIndex, ".txt", 
                    sep = ""), sep = "\t", row.names = FALSE)
                  Rep[[Index]] = readLines(paste(RunFile, "agemat.rep", 
                    sep = ""))
                }
            }
        }
        IcRecord = rbind(IcRecord, IcVec)
        StateRecord = rbind(StateRecord, ParamVecOpt)
        capture.output(list(IcRecord = IcRecord, StateRecord = StateRecord), 
            file = paste(SaveFile, "Stepwise - Record.txt", sep = ""))
        Max = which.max(IcVec)
        if (all(ParamMat[Max, ] == ParamVecOpt)) 
            Stop = TRUE
        ParamVecOpt = ParamMat[Max, ]
        CurrentMinusAge = ParamVecOpt[length(ParamVecOpt) - 1]
        SearchMat[length(ParamVecOpt) - 1, 1:7] = c(CurrentMinusAge, 
            CurrentMinusAge - 10, CurrentMinusAge - 4, CurrentMinusAge - 
                1, CurrentMinusAge + 1, CurrentMinusAge + 4, 
            CurrentMinusAge + 10)
        SearchMat[length(ParamVecOpt) - 1, 1:7] = ifelse(SearchMat[length(ParamVecOpt) - 
            1, 1:7] < MinAge, NA, SearchMat[length(ParamVecOpt) - 
            1, 1:7])
        CurrentPlusAge = ParamVecOpt[length(ParamVecOpt)]
        SearchMat[length(ParamVecOpt), 1:7] = c(CurrentPlusAge, 
            CurrentPlusAge - 10, CurrentPlusAge - 4, CurrentPlusAge - 
                1, CurrentPlusAge + 1, CurrentPlusAge + 4, CurrentPlusAge + 
                10)
        SearchMat[length(ParamVecOpt), 1:7] = ifelse(SearchMat[length(ParamVecOpt), 
            1:7] > MaxAge, NA, SearchMat[length(ParamVecOpt), 
            1:7])
        writeLines(Rep[[Max]], con = paste(SaveFile, "agemat.rep", 
            sep = ""))
        PlotOutputFn(Data = Data, MaxAge = MaxAge, SaveFile = SaveFile, 
            PlotType = "JPG")
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
