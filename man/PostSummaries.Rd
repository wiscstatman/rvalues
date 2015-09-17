\name{PostSummaries}
\alias{PostSummaries}
\title{
R-values from posterior summary quantities
}
\description{
Computes r-values assuming that, for each parameter of interest, the
user supplies a value for the posterior mean and the 
posterior standard deviation. The assumption here is
that the posterior distributions are Normal.
}
\usage{
PostSummaries(post.means, post.sds, hypers = NULL, qtheta = NULL, alpha.grid = NULL, 
               ngrid = NULL, smooth = 0)
}
\arguments{
  \item{post.means}{
a vector of posterior means
}
  \item{post.sds}{
a vector of posterior standard deviations
}
  \item{hypers}{
a list with two elements: mean and sd. These represent
the parameters in the (Normal) prior which was used
to generate the posterior means and sds. If hypers
is not supplied then one must supply the quantile
function qtheta.
}
  \item{qtheta}{
a function which returns the quantiles (for upper tail probs.) of 
theta. If this is not supplied,
the hyperparameter must be supplied.
}
  \item{alpha.grid}{
grid of values in (0,1); used for the discrete approximation
approach for computing r-values.
}
  \item{ngrid}{
number of grid points for alpha.grid; only relevant when \code{alpha.grid = NULL}
}
  \item{smooth}{
either \code{smooth="none"} or \code{smooth} takes 
a value between 0 and 10; this determines the level of smoothing applied to the
estimate of \eqn{\lambda(\alpha)}; if \code{smooth} is given a number, the
number is used as the \code{bass} argument in \code{\link{supsmu}}.
}
}
\value{
An object of class "rvals"
}
\author{
Nicholas Henderson and Michael Newton
}
\seealso{
\code{\link{rvalues}}
}
\examples{
n <- 500
theta <- rnorm(n)
sig_sq <- rgamma(n,shape=1,scale=1)
X <- theta + sqrt(sig_sq)*rnorm(n)

pm <- X/(sig_sq + 1)
psd <- sqrt(sig_sq/(sig_sq + 1))

rvs <- PostSummaries(post.means=pm,post.sds=psd,hypers=list(mean=0,sd=1))
hist(rvs$rvalues)
}
\keyword{ htest }
\keyword{ models }
