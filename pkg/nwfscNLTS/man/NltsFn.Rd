\name{NltsFn}
\alias{NltsFn}
\title{
  Non-linear time series -- crossvalidation correlation
}
\description{
  Calculates correlation for crossvalidation 
}
\usage{
NltsFn(Y, PredInterval, Nembed, Theta = NA, Method)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
  Time series
}
  \item{PredInterval}{
  Prediction interval, i.e. 1 represents predicting one observation ahead
}
  \item{Nembed}{
  Embedding dimension (if using simplex)
}
  \item{Theta}{
  Magnitude of non-linearity (if using S-map)
}
  \item{Method}{
  Method to use, either 'Smap' or 'Simplex'
}
}
\details{
}
\value{
  Returns a table comparing crossvalidation predictions and observations, as well as the their correlation
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
