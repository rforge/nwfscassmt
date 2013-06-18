\name{SimPredPreyFn}
\alias{SimPredPreyFn}
\title{
 Simulate predator-prey data 
}
\description{
  This function simulates using a fifference equation approximation to Lotka-Volterra dynamics
}
\usage{
SimPredPreyFn(Nobs, Nt = 100, A = 0.04, B = 0.05, C = 0.2, E = 1)
}
\arguments{
  \item{Nobs}{
  Length of time series data vectors
}
  \item{Nt}{
  Number of finite-difference approximations used between time series observations
}
  \item{A}{
  Prey intrinsic growth rate
}
  \item{B}{
  Prey mortality rate per capita due to per capita predation
}
  \item{C}{
  Predator assimilation rate
}
  \item{E}{
  Predator mortality rate
}
}
\details{
}
\value{
  List with observations for prey (X) and predator (Y)
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

# Generate Lotka-Volterra timeseries
Y = SimPredPreyFn(Nobs=Nobs, Nt=100, A=0.4, B=0.5, C=0.2, E=1)$Y
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
