\name{OverlapCurve}
\alias{OverlapCurve}
\title{
Overlap Curve
}
\description{
Estimates the expected proportion of units in the top fraction
and those deemed to be in the top fraction by the r-value
procedure. If plot=TRUE, the curve is plotted before the
estimated function is returned.
}
\usage{
OverlapCurve(object, plot = TRUE, xlim, ylim, xlab, ylab, main, ...)
}
\arguments{
  \item{object}{
An object of class "rvals"
}
  \item{plot}{
logical. If TRUE, the estimated overlap curve is plotted.
}
  \item{xlim, ylim}{
x and y - axis limits for the plot
}
  \item{xlab,ylab}{
x and y - axis labels
}
  \item{main}{
the title of the plot
}
  \item{\dots}{
 additional arguments to \code{\link{plot.default}}
}
}
\details{
For parameters of interest \eqn{\theta_1,...,\theta_n} and corresponding
r-values \eqn{r_1,...,r_n}, the overlap at a particular value of \eqn{\alpha} 
is defined to be 
\deqn{ overlap(\alpha) = P(\theta_i \ge \theta_\alpha, r_i \le \alpha), }
where the threshold \eqn{\theta_\alpha} is the upper-\eqn{\alpha}th quantile of 
the distribution of the \eqn{\theta_i} (i.e., \eqn{P(\theta_i \ge \theta_\alpha) = \alpha}).
\code{OverlapCurve} estimates this overlap
for values of alpha across (0,1) and plots (if \code{plot=TRUE})
the resulting curve.
}
\value{
A function returning estimated overlap values.
}
\references{
Henderson, N.C. and Newton, M.A. (2014) Making the Cut: Improved Ranking and Selection 
for Large-Scale Inference. 
\url{http://arxiv.org/abs/1312.5776}
}
\author{
Nicholas Henderson and Michael Newton
}
\examples{
n <- 500
theta <- rnorm(n)
ses <- sqrt(rgamma(n,shape=1,scale=1))
XX <- theta + ses*rnorm(n)
dd <- cbind(XX,ses)

rvs <- rvalues(dd, family = gaussian)

OverlapCurve(rvs, cex.main = 1.5)
}
\keyword{ dplot }
\keyword{ htest }
