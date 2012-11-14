\name{PlotOutputFn}
\alias{PlotOutputFn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PlotOutputFn(Data, MaxAge, SaveFile, PlotType = "PDF")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{
%%     ~~Describe \code{Data} here~~
}
  \item{MaxAge}{
%%     ~~Describe \code{MaxAge} here~~
}
  \item{SaveFile}{
%%     ~~Describe \code{SaveFile} here~~
}
  \item{PlotType}{
%%     ~~Describe \code{PlotType} here~~
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
function (Data, MaxAge, SaveFile, PlotType = "PDF") 
{
    Nreaders = ncol(Data) - 1
    Ages = Nages = MaxAge + 1
    Rep = scan(paste(SaveFile, "agemat.rep", sep = ""), comment.char = "\%", 
        what = "character", quiet = TRUE)
    Grep = grep("reader#", Rep)
    MisclassArray = array(NA, dim = c(Nreaders, Ages, Ages), 
        dimnames = list(paste("Reader", 1:Nreaders), paste("TrueAge", 
            0:MaxAge), paste("EstAge", 0:MaxAge)))
    for (i in 1:Nreaders) {
        MisclassArray[i, , ] = matrix(as.numeric(Rep[Grep[i] + 
            1 + 1:(Ages^2)]), ncol = Ages, byrow = TRUE)
    }
    Grep = grep("age-structure", Rep)
    AgeStruct = matrix(as.numeric(Rep[Grep[1] + 7 + 1:(2 * Ages)]), 
        ncol = 2, byrow = TRUE)
    Grep = grep("age", Rep)[3]
    Temp = Rep[Grep + 1:(5 * Nages * Nreaders)]
    ErrorAndBiasArray = array(as.numeric(Temp), dim = c(5, Nages, 
        Nreaders), dimnames = list(c("Reader", "True_Age", "CV", 
        "SD", "Expected_age"), paste("Age", 0:MaxAge), paste("Reader", 
        1:Nreaders)))
    AgeProbs = array(NA, dim = c(nrow(Data), Ages), dimnames = list(paste("Otolith", 
        1:nrow(Data)), paste("TrueAge", 0:MaxAge)))
    OtI = AgeI = ReadI = 1
    for (OtI in 1:nrow(Data)) {
        for (AgeI in 1:Ages) {
            AgeProbs[OtI, AgeI] = AgeStruct[AgeI, 2]
            for (ReadI in 1:Nreaders) {
                if (Data[OtI, ReadI + 1] != -999) {
                  AgeRead = Data[OtI, ReadI + 1]
                  AgeProbs[OtI, AgeI] = AgeStruct[AgeI, 2] * 
                    (MisclassArray[ReadI, AgeI, AgeRead + 1])^Data[OtI, 
                      1]
                }
            }
        }
    }
    TrueAge = apply(AgeProbs, MARGIN = 1, FUN = function(Vec) {
        order(Vec[-length(Vec)], decreasing = TRUE)[1]
    })
    Temp = ifelse(Data[, -1] == -999, NA, Data[, -1])
    Temp = tapply(ifelse(is.na(Temp), 0, 1), INDEX = Temp, FUN = sum)
    Prop = rep(0, MaxAge + 1)
    Prop[as.numeric(names(Temp)) + 1] = Temp
    cbind(0:MaxAge, Prop, round(AgeStruct[, 2] * sum(Prop), 1))
    cbind(0:MaxAge, round(Prop/sum(Prop), 3), round(AgeStruct[, 
        2], 3))
    Plot = function() {
        par(mar = c(3, 3, 2, 0), mgp = c(1.5, 0.25, 0), tck = -0.02, 
            oma = c(0, 0, 0, 0) + 0.1)
        plot(x = AgeStruct[, 1], y = AgeStruct[, 2], type = "s", 
            lwd = 2, xlab = "Age", ylab = "Prop", main = "Estimated=Black, Observed=Red")
        hist(ifelse(Data[, -1] == -999, NA, Data[, -1]), add = TRUE, 
            freq = FALSE, breaks = seq(0, MaxAge, by = 1), col = rgb(red = 1, 
                green = 0, blue = 0, alpha = 0.3))
    }
    if (PlotType == "PDF") {
        pdf(paste(SaveFile, "Estimated vs Observed Age Structure.pdf", 
            sep = ""), width = 6, height = 6)
        Plot()
        dev.off()
    }
    if (PlotType == "JPG") {
        jpeg(paste(SaveFile, "Estimated vs Observed Age Structure.jpg", 
            sep = ""), width = 6, height = 6, units = "in", res = 200)
        Plot()
        dev.off()
    }
    Ncol = ceiling(sqrt(Nreaders))
    Nrow = ceiling(Nreaders/Ncol)
    Plot = function() {
        par(mfrow = c(Nrow, Ncol), mar = c(3, 3, 2, 0), mgp = c(1.5, 
            0.25, 0), tck = -0.02, oma = c(0, 0, 5, 0) + 0.1)
        for (ReadI in 1:Nreaders) {
            Temp = cbind(TrueAge, Data[, ReadI + 1] + 0.5)
            Temp = Temp[which(Data[, ReadI + 1] != -999), ]
            plot(x = Temp[, 1], y = Temp[, 2], ylim = c(0, MaxAge), 
                xlim = c(0, MaxAge), col = rgb(red = 0, green = 0, 
                  blue = 0, alpha = 0.2), xlab = "Mode predicted age | parameters", 
                ylab = "Read age", lwd = 2, main = paste("Reader", 
                  ReadI), pch = 21, cex = 0.2)
            lines(x = c(0, MaxAge), y = c(0, MaxAge), lwd = 1, 
                lty = "dashed")
            lines(x = ErrorAndBiasArray["True_Age", , ReadI], 
                y = ErrorAndBiasArray["Expected_age", , ReadI], 
                type = "l", col = "red", lwd = 1)
            lines(x = ErrorAndBiasArray["True_Age", , ReadI], 
                y = ErrorAndBiasArray["SD", , ReadI], type = "l", 
                col = "blue", lwd = 1)
            lines(x = ErrorAndBiasArray["True_Age", , ReadI], 
                y = ErrorAndBiasArray["Expected_age", , ReadI] + 
                  2 * ErrorAndBiasArray["SD", , ReadI], type = "l", 
                col = "red", lwd = 1, lty = "dashed")
            lines(x = ErrorAndBiasArray["True_Age", , ReadI], 
                y = ErrorAndBiasArray["Expected_age", , ReadI] - 
                  2 * ErrorAndBiasArray["SD", , ReadI], type = "l", 
                col = "red", lwd = 1, lty = "dashed")
        }
        mtext(side = 3, outer = TRUE, text = "Reads(dot), Sd(blue), expected_read(red solid line),\n and 95\% CI for expected_read(red dotted line)", 
            line = 1)
    }
    if (PlotType == "PDF") {
        pdf(paste(SaveFile, "True vs Reads (by reader).pdf", 
            sep = ""), width = Ncol * 3, height = Nrow * 3)
        Plot()
        dev.off()
    }
    if (PlotType == "JPG") {
        jpeg(paste(SaveFile, "True vs Reads (by reader).jpg", 
            sep = ""), width = Ncol * 3, height = Nrow * 3, units = "in", 
            res = 200)
        Plot()
        dev.off()
    }
    Nll = as.numeric(scan(paste(SaveFile, "agemat.par", sep = ""), 
        comment.char = "\%", what = "character", quiet = TRUE)[11])
    Df = as.numeric(scan(paste(SaveFile, "agemat.par", sep = ""), 
        comment.char = "\%", what = "character", quiet = TRUE)[6])
    n = sum(ifelse(Data[, -1] == -999, 0, 1))
    Aic = 2 * Nll + 2 * Df
    Aicc = Aic + 2 * Df * (Df + 1)/(n - Df - 1)
    Bic = 2 * Nll + Df * log(n)
    Output = list(Aic = Aic, Aicc = Aicc, Bic = Bic)
    return(Output)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
