\name{FDRCurve}
\alias{FDRCurve}
\title{
FDR Curve
}
\description{
Estimates the expected proportion of misclassified units
when using a given r-value threshold. If plot=TRUE, the 
curve is plotted before the estimated function is returned.
}
\usage{
FDRCurve(object, q, threshold = 1, plot = TRUE, xlim, ylim, xlab, ylab, main, ...)
}
\arguments{
  \item{object}{
An object of class "rvals"
}
  \item{q}{
A value in between 0 and 1; the desired level of FDR control.
}
  \item{threshold}{
The r-value threshold.
}
  \item{plot}{
logical; if TRUE, the estimated FDR curve is plotted.
}
  \item{xlim,ylim}{
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
Consider parameters of interest \eqn{(\theta_1,...,\theta_n)} with an effect of size of interest 
\eqn{\tau}. That is, a unit is taken to be "null" if \eqn{\theta_i \le \tau} and
taken to be "non-null" if \eqn{\theta_i > \tau}.

For r-values \eqn{r_1,...,r_n} and a procedure which "rejects" units 
satisfying \eqn{r_i \le c}, the FDR is defined to be
\deqn{FDR(c) = P(\theta_i < \tau,r_i \le c)/P(r_i \le c).}
\code{FDRCurve} estimates \eqn{FDR(c)}
for values of \eqn{c} across (0,1) and plots (if \code{plot=TRUE})
the resulting curve.
}
\value{
A list with the following two components
  \item{fdrcurve}{A function which returns the estimated FDR for each r-value threshold.}
  \item{Rval.cut}{The largest r-value cutoff which still gives an estimated FDR less than q.}
}
\author{
Nicholas Henderson and Michael Newton
}
\seealso{
\code{\link{OverlapCurve}}
}
\examples{
n <- 500
theta <- rnorm(n)
ses <- sqrt(rgamma(n,shape=1,scale=1))
XX <- theta + ses*rnorm(n)
dd <- cbind(XX,ses)

rvs <- rvalues(dd, family = gaussian)

FDRCurve(rvs, q = .1, threshold = .3, cex.main = 1.5)
}
\keyword{ dplot }
\keyword{ hstat }