\name{NLTS-package}
\alias{NLTS}
\docType{package}
\title{
  Non-linear time series
}
\description{
  This package includes implementations of simplex and S-map algorithms for projecting time series data in the presence of non-linear (state-dependent) dynamics 
}
\details{
\tabular{ll}{
Package: \tab NLTS\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-11-15\cr
License: \tab GPL-2\cr
}
~~ An overview of how to use the package, including the most important functions ~~
}
\author{
James Thorson (JamesT.esq+nlts@gmail.com)
}
\references{
Sugihara, G., and May, R.M. 1990. Nonlinear forecasting as a way of distinguishing chaos from measurement error in time series. Nature 334: 734-741.
Sugihara, G. 1994. Nonlinear forecasting for the classification of natural time series. Philosophical Transactions of the Royal Society B: Biological Sciences 348: 477-495.
}
\keyword{ package }
\seealso{
}
\examples{
# Length of time series
Nobs = 100

# Generate tent function timeseries
Y = SimTentFn(Nobs=Nobs, S=1.75)
X = seq(0,1, length=1e4); plot(x=X, y=sapply(X, FUN=TentFn, S=1.75))
plot(y=Y[-1], x=Y[-Nobs])
plot(x=Y[-Nobs], y=Y[-1])

# Generate Lotka-Volterra timeseries
Y = SimPredPreyFn(Nobs=Nobs, Nt=100, A=0.4, B=0.5, C=0.2, E=1)$Y

# Simplex
Output = NltsFn(Y, PredInterval=2, Nembed=1, Method="Simplex")
EmbedFn(c(Y,NA,NA), PredInterval=2, Candidates=1:10)
NltsPred(c(Y,NA,NA), PredInterval=2, Nembed=2, PredNum=length(Y)+2, Method="Simplex")  # PredInterval is the number of years in the future (1 is next year, etc)

# S-map
Output = NltsFn(c(Y,NA,NA), PredInterval=2, Nembed=2, Method="Smap", Theta=1)
ThetaFn(c(Y,NA,NA), PredInterval=2, Nembed=2)
NltsPred(c(Y,NA,NA), PredInterval=2, Nembed=4, PredNum=length(Y)+1, Method="Smap", Theta=1)

}
