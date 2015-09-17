\name{npmixapply}
\alias{npmixapply}
\title{
Apply Functions over estimated unit-specific posterior
distributions
}
\description{
Using a nonparametric estimate of the mixing distribution, computes
a posterior quantity of interest for each unit.
}
\usage{
npmixapply(object, FUN, ...)
}
\arguments{
  \item{object}{ 
  an object of class "npmix"
}
  \item{FUN}{
  a user provided function
}
  \item{\dots}{
  optional arguments to FUN  
}
}
\details{
\code{object} is an object of class "npmix" containing 
a nonparametric estimate of the mixing distribution \eqn{F}
in the following two-level sampling model
\eqn{ X_i|\theta_i} ~ \eqn{p(x|\theta_i,\eta_i)}
and \eqn{\theta_i} ~ \eqn{F} for \eqn{i = 1,...,n}.
  
Using \code{npmixapply(object, f)}, then returns the 
posterior expectation of \eqn{f}:
\eqn{ E[ f(\theta_i) | X_i, \eta_i] }, for \eqn{i = 1,...,n}.  
}
\value{
a vector with length equal to \eqn{n}
}
\author{
Nicholas Henderson
}
\seealso{
\code{\link{npmle}}
}
\examples{
\dontrun{
data(hiv)
npobj <- npmle(hiv, family = gaussian, maxiter = 4)

### Compute unit-specific posterior means
pmean <- npmixapply(npobj, function(x) { x })

### Compute post. prob that \theta_i < .1
pp <- npmixapply(npobj, function(x) { x < .1})
}
}
\keyword{hstat}
