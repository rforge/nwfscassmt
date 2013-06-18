\name{SimTentFn}
\alias{SimTentFn}
\title{
  Simulate tent function
}
\description{
  This simulates time series data from a tent function, a simple functionm that generates 'nonlinear' dynamics 
}
\usage{
  SimTentFn(Nobs, S)
}
\arguments{
  \item{Nobs}{
  Time series length
}
  \item{S}{
  Tent height (producing non-linear dynamics given S>1.75)
}
}
\details{
}
\value{
  Generates a vector of observations
}
\references{
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

# Generate tent function timeseries
Y = SimTentFn(Nobs=Nobs, S=1.75)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
