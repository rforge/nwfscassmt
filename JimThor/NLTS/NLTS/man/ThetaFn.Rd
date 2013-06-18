\name{ThetaFn}
\alias{ThetaFn}
\title{
  Optimize theta
}
\description{
  Estimate optimal degree of non-linearity for S-map
}
\usage{
ThetaFn(Y, PredInterval, ThetaSet = c(0, exp(seq(-3, 3))), Nembed)
}
\arguments{
  \item{Y}{
  Time series
}
  \item{PredInterval}{
  Prediction interval, i.e. 1 represents predicting one observation ahead
}
  \item{ThetaSet}{
  Vector containing candidate set of theta's to seach across
}
  \item{Nembed}{
  Embedding dimension
}
}
\details{
}
\value{
  Returns a list with the crossvalidation correlation for each candidate theta, as well as the theta that maximizes that correlation
}
\references{
Sugihara, G., and May, R.M. 1990. Nonlinear forecasting as a way of distinguishing chaos from measurement error in time series. Nature 334: 734-741.
Sugihara, G. 1994. Nonlinear forecasting for the classification of natural time series. Philosophical Transactions of the Royal Society B: Biological Sciences 348: 477-495.
}
\author{
James Thorson
}
\note{
}
\seealso{
}
\examples{
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
