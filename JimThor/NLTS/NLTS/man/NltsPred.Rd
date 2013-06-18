\name{NltsPred}
\alias{NltsPred}
\title{
  Non-linear time series -- compute predictions
}
\description{
  Function to compute predictions for simplex and S-map
}
\usage{
NltsPred(Y, PredInterval, Nembed, PredNum, Method, Theta = NA)
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
  Embedding dimension
}
  \item{PredNum}{
  Index of observation to predict
}
  \item{Method}{
  Method to use, either 'Smap' or 'Simplex'
}
  \item{Theta}{
  Magnitude of non-linearity (if using S-map)
}
}
\details{
}
\value{
  Returns prediction of Y[PredNum]
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
