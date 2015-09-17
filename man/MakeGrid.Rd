\name{MakeGrid}
\alias{MakeGrid}
\title{
Grid Construction
}
\description{
Computes a grid of points on the interval (0,1). This function is
useful for constructing the "alpha-grid" used in various r-value computations.
}
\usage{
MakeGrid(nunits, type = "log", ngrid = NULL, lower = 1/nunits, upper = 1 - lower)
}
\arguments{
  \item{nunits}{
The number of units in the data for which r-values
are to be calculated.
}
  \item{type}{
The type of grid; type can be set to \code{type="uniform"},
\code{type="log"}, or \code{type="log.symmetric"}.
}
  \item{ngrid}{
a number specifying the number of grid points
}
  \item{lower}{
the smallest grid point; must be greater than zero  
}
  \item{upper}{
the largest grid point; must be less than one
}
}
\details{
If \code{nunits} \eqn{\le 1000}, the default number of grid points is
equal to \code{nunits}. When \code{nunits} \eqn{> 1000}, the default 
number of grid points is determined by 
\deqn{1000 + 25*log(nunits - 1000)*(nunits - 1000)^{1/4}}.
}
\value{
A vector of grid points in (0,1).
}
\author{
Nicholas Henderson and Michael Newton
}
\seealso{
\code{\link{rvalues}}
}
\examples{
alpha.grid <- MakeGrid(1000,type="uniform",ngrid=200)

log.grid <- MakeGrid(40,type="log")
log.grid
hist(log.grid)
}
\keyword{ math }
