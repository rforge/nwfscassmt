\name{EmbedFn}
\alias{EmbedFn}
\title{
  Find optimal simplex 'embedding' dimension
}
\description{
  Function that returns the embedding dimention for simplex that maximizes cross-validation correlation
}
\usage{
  EmbedFn(Y, PredInterval, Candidates = 1:5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
  Time series
}
  \item{PredInterval}{
  Prediction interval, i.e. 1 represents predicting one observation ahead
}
  \item{Candidates}{
  Candidate embedding dimensions
}
}
\details{
}
\value{
  Returns an integer from the 'Candidates' which maximizes crossvalidation correlation
}
\references{
Sugihara, G., and May, R.M. 1990. Nonlinear forecasting as a way of distinguishing chaos from measurement error in time series. Nature 334: 734-741.
}
\author{
James Thorson
}
\note{
}

\seealso{
}
\examples{
# Time series length
Nobs = 100

# Generate Lotka-Volterra timeseries
Y = SimPredPreyFn(Nobs=Nobs, Nt=100, A=0.4, B=0.5, C=0.2, E=1)$Y

# Simplex
Output = NltsFn(Y, PredInterval=2, Nembed=1, Method="Simplex")
EmbedFn(c(Y,NA,NA), PredInterval=2, Candidates=1:10)

}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
